#include <stdio.h>
#include <armadillo>

#include "MeshData.h"
#include "VTKExporter.h"

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

int main() {
    voxel2tet::MeshData *mesh = createmesh(10);
    voxel2tet::VTKExporter exporter(mesh);
    exporter.WriteData("/tmp/testXYZ.vtp");

}
