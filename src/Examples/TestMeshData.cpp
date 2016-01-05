#include <stdio.h>
#include <array>
#include <armadillo>

#include "MeshManipulations.h"
#include "MeshComponents.h"
#include "VTKExporter.h"
#include "PhaseEdge.h"
#include "Options.h"

voxel2tet::MeshData *createmesh(int n)
{
    // n is the number of vertices, i.e. n-1 is the number of boundary elements
    //
    arma::vec coords = arma::linspace(0.0, double(n)-1, double(n));
    arma::mat vertices (n*n, 3, arma::fill::zeros);

    int r=0;
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++) {
            vertices.row(r) = arma::rowvec({coords(j), coords(i), coords(0)});
            r++;
        }
    }

    r=0;
    arma::Mat<int> squares((n-1)*(n-1), 4);
    for (int i=0; i<n-1; i++){
        for (int j=0; j<n-1; j++) {
            squares.row(r) = arma::Row<int>({(i*n+j), (i*n+j+1), ((i+1)*n+j), ((i+1)*n+j+1)});
            r++;
        }
    }

    r=0;
    arma::Mat<int> triangles((n-1)*(n-1)*2, 3);
    for (unsigned int i=0; i<squares.n_rows; i++) {
        triangles.row(r) = arma::Row<int>({squares(i,0), squares(i,1), squares(i,3)});
        r++;
        triangles.row(r) = arma::Row<int>({squares(i,0), squares(i,3), squares(i,2)});
        r++;
    }

    //vertices.print();
    //squares.print();
    //triangles.print();

    voxel2tet::BoundingBoxType BoundingBox;
    BoundingBox.maxvalues[0] = double(n);
    BoundingBox.maxvalues[1] = double(n);
    BoundingBox.maxvalues[2] = double(n);
    BoundingBox.minvalues[0] = 0.0;
    BoundingBox.minvalues[1] = 0.0;
    BoundingBox.minvalues[2] = 0.0;

    voxel2tet::MeshData *Mesh = new voxel2tet::MeshData(BoundingBox);
    for (unsigned int i=0; i<triangles.n_rows; i++) {
        voxel2tet::TriangleType *t =
            Mesh->AddTriangle({vertices(triangles(i,0), 0), vertices(triangles(i,0), 1), vertices(triangles(i,0), 2)},
                         {vertices(triangles(i,1), 0), vertices(triangles(i,1), 1), vertices(triangles(i,1), 2)},
                         {vertices(triangles(i,2), 0), vertices(triangles(i,2), 1), vertices(triangles(i,2), 2)});
        t->InterfaceID = i;
    }
    return Mesh;
}

int main( int argc, char *argv[] ) {
    int vertexsidecount = 4;
    std::map <std::string, std::string> DefaultOptions;
    voxel2tet::Options *Options = new voxel2tet::Options(argc, argv, DefaultOptions);

    voxel2tet::MeshManipulations *mesh = static_cast <voxel2tet::MeshManipulations*> (createmesh(vertexsidecount)) ;

    voxel2tet::PhaseEdge West(Options), East(Options), North(Options), South(Options);

    for (voxel2tet::VertexType* v: mesh->Vertices) {
        std::array<double, 3> c=v->get_c();

        if (c[0] == 0) v->AddPhaseEdge(&West);
        if (c[1] == 0) v->AddPhaseEdge(&South);
        if (c[0] > double(vertexsidecount-1) - 1e-8) v->AddPhaseEdge(&East);
        if (c[1] > double(vertexsidecount-1) - 1e-8) v->AddPhaseEdge(&North);

    }

    for (voxel2tet::EdgeType* e: mesh->Edges) {
         if ((e->Vertices.at(0)->get_c(1) == 0.0) && (e->Vertices.at(1)->get_c(1) == 0.0)){
             South.EdgeSegments.push_back(e->Vertices);
         }
         if ((e->Vertices.at(0)->get_c(0) == 0.0) && (e->Vertices.at(1)->get_c(0) == 0.0)){
             West.EdgeSegments.push_back(e->Vertices);
         }

         if ((e->Vertices.at(0)->get_c(1) > double(vertexsidecount-1) - 1e-8) && (e->Vertices.at(1)->get_c(1) > double(vertexsidecount-1) - 1e-8)){
             North.EdgeSegments.push_back(e->Vertices);
         }

         if ((e->Vertices.at(0)->get_c(0) > double(vertexsidecount-1) - 1e-8) && (e->Vertices.at(1)->get_c(0) > double(vertexsidecount-1) - 1e-8)){
             East.EdgeSegments.push_back(e->Vertices);
         }
    }


    mesh->Vertices.at(11)->set_c(2.5, 0);

    mesh->ExportVTK("/tmp/test2D_0.vtp");
    int iter=1;
    bool EdgeCollapsed = true;
    while (EdgeCollapsed) {
        EdgeCollapsed=false;
        for (voxel2tet::EdgeType* e: mesh->Edges) {

            if (!mesh->CollapseEdge(e,0)) {
                EdgeCollapsed=mesh->CollapseEdge(e,1);
            } else {
                EdgeCollapsed=true;
            }


            if (EdgeCollapsed) {
                std::ostringstream FileName;
                FileName << "/tmp/test2D_" << iter << ".vtp";
                mesh->ExportVTK( FileName.str() );
                iter++;
                break;
            }
        }
    }



}
