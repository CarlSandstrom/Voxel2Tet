#include "TetGenCaller.h"

namespace voxel2tet
{

TetGenCaller::TetGenCaller()
{

}

void TetGenCaller :: UpdateVertexMapping()
{
    this->VertexMapFromID.clear();
    this->VertexMapFromTetgen.clear();

    std::vector<VertexType *> VertexIDsFromTriangles;

    for (TriangleType *t: this->Mesh->Triangles) {
        for (VertexType *v: t->Vertices)
        VertexIDsFromTriangles.push_back(v);
    }

    std::sort(VertexIDsFromTriangles.begin(), VertexIDsFromTriangles.end());

    VertexIDsFromTriangles.erase( std::unique(VertexIDsFromTriangles.begin(), VertexIDsFromTriangles.end()), VertexIDsFromTriangles.end() );

    for (unsigned int i=0; i<VertexIDsFromTriangles.size(); i++) {
        this->VertexMapFromID[VertexIDsFromTriangles.at(i)] = i;
        this->VertexMapFromTetgen[i] = VertexIDsFromTriangles.at(i);
    }

}

void TetGenCaller :: CopyMeshFromSelf(tetgenio *in)
{

    this->UpdateVertexMapping();
    tetgenio::facet *f;
    tetgenio::polygon *p;

    in->initialize();
    in->firstnumber = 0;

    // Copy mesh vertices
    in->numberofpoints = this->VertexMapFromID.size();
    in->pointlist = new REAL[in->numberofpoints * 3];

    unsigned int cid=0; // Coordinate id. (0, 1, 2), (3, 4, 5), (6, 7, 8) ...

    std::map<VertexType *, int>::iterator iter;
    for (iter = this->VertexMapFromID.begin(); iter != this->VertexMapFromID.end(); iter++) {
         VertexType *v=iter->first;
         for (double coordvalue: v->get_c()) {
             in->pointlist[cid]=coordvalue;
             cid++;
         }
    }

    // Copy mesh triangles
    in->numberoffacets = this->Mesh->Triangles.size();
    in->facetlist = new tetgenio::facet[in->numberoffacets]; // List of triangles
    in->facetmarkerlist = new int[in->numberoffacets]; // List of tag for each triangle

    unsigned int tid=0;
    for (TriangleType *t: this->Mesh->Triangles) {

        f = &in->facetlist[tid];

        f->numberofpolygons = 1;
        f->holelist = (REAL *) NULL;
        f->numberofholes = 0;
        f->polygonlist = new tetgenio::polygon[f->numberofpolygons];

        p = &f->polygonlist[0];
        p->numberofvertices = 3;
        p->vertexlist = new int[p->numberofvertices];

        for (int i=0; i<3; i++) {
            VertexType *v = t->Vertices.at(i);
            int TetGenID = this->VertexMapFromID[v];
            p->vertexlist[i] = TetGenID;
        }

        in->facetmarkerlist[tid] = t->InterfaceID;

        tid++;

    }

}

void TetGenCaller :: CopyMeshToSelf(tetgenio *io)
{
    //io->pointlist
    for (int i=0; i<io->numberoftrifaces; i++) {
        int *triface[3];
        int marker = io->trifacemarkerlist[i];
        triface = &io->trifacelist[3*i];
        printf("(%u, %u, %u) : %u\n", triface[0], triface[1], triface[2], marker);
    }
}

void TetGenCaller :: EmbarrassingTestExample()
{
    tetgenio in, out;
    tetgenio::facet *f;
    tetgenio::polygon *p;
    int i;

    // All indices start from 1.
    in.firstnumber = 1;

    in.numberofpoints = 9;
    in.pointlist = new REAL[in.numberofpoints * 3];
    in.pointlist[0]  = 0;  // node 1.
    in.pointlist[1]  = 0;
    in.pointlist[2]  = 0;
    in.pointlist[3]  = 2;  // node 2.
    in.pointlist[4]  = 0;
    in.pointlist[5]  = 0;
    in.pointlist[6]  = 2;  // node 3.
    in.pointlist[7]  = 2;
    in.pointlist[8]  = 0;
    in.pointlist[9]  = 0;  // node 4.
    in.pointlist[10] = 2;
    in.pointlist[11] = 0;
    // Set node 5, 6, 7, 8.
    for (i = 4; i < 8; i++) {
      in.pointlist[i * 3]     = in.pointlist[(i - 4) * 3];
      in.pointlist[i * 3 + 1] = in.pointlist[(i - 4) * 3 + 1];
      in.pointlist[i * 3 + 2] = 12;
    }
    in.pointlist[24] = .1;
    in.pointlist[25] = .1;
    in.pointlist[26] = .1;

    in.numberoffacets = 7;
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];

    // Facet 1. The leftmost facet.
    f = &in.facetlist[0];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    p = &f->polygonlist[0];
    p->numberofvertices = 4;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0] = 1;
    p->vertexlist[1] = 2;
    p->vertexlist[2] = 3;
    p->vertexlist[3] = 4;

    // Facet 2. The rightmost facet.
    f = &in.facetlist[1];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    p = &f->polygonlist[0];
    p->numberofvertices = 4;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0] = 5;
    p->vertexlist[1] = 6;
    p->vertexlist[2] = 7;
    p->vertexlist[3] = 8;

    // Facet 3. The bottom facet.
    f = &in.facetlist[2];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    p = &f->polygonlist[0];
    p->numberofvertices = 4;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0] = 1;
    p->vertexlist[1] = 5;
    p->vertexlist[2] = 6;
    p->vertexlist[3] = 2;

    // Facet 4. The back facet.
    f = &in.facetlist[3];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    p = &f->polygonlist[0];
    p->numberofvertices = 4;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0] = 2;
    p->vertexlist[1] = 6;
    p->vertexlist[2] = 7;
    p->vertexlist[3] = 3;

    // Facet 5. The top facet.
    f = &in.facetlist[4];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    p = &f->polygonlist[0];
    p->numberofvertices = 4;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0] = 3;
    p->vertexlist[1] = 7;
    p->vertexlist[2] = 8;
    p->vertexlist[3] = 4;

    // Facet 6. The front facet.
    f = &in.facetlist[5];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    p = &f->polygonlist[0];
    p->numberofvertices = 3;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0] = 4;
    p->vertexlist[1] = 8;
    p->vertexlist[2] = 1;

    // Facet 7. The front facet.
    f = &in.facetlist[6];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    p = &f->polygonlist[0];
    p->numberofvertices = 3;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0] = 1;
    p->vertexlist[1] = 5;
    p->vertexlist[2] = 8;

    // Set 'in.facetmarkerlist'

    in.facetmarkerlist[0] = -1;
    in.facetmarkerlist[1] = -2;
    in.facetmarkerlist[2] = 0;
    in.facetmarkerlist[3] = 0;
    in.facetmarkerlist[4] = 0;
    in.facetmarkerlist[5] = 0;
    in.facetmarkerlist[6] = 0;

    // Output the PLC to files 'barin.node' and 'barin.poly'.
    in.save_nodes("/tmp/barin");
    in.save_poly("/tmp/barin");

    // Tetrahedralize the PLC. Switches are chosen to read a PLC (p),
    //   do quality mesh generation (q) with a specified quality bound
    //   (1.414), and apply a maximum volume constraint (a0.1).

    tetrahedralize("pq1.414a0.1", &in, &out);

    // Output mesh to files 'barout.node', 'barout.ele' and 'barout.face'.
    out.save_nodes("/tmp/barout");
    out.save_elements("/tmp/barout");
    out.save_faces("/tmp/barout");
}

void TetGenCaller :: Execute()
{

    tetgenio in, out;

    CopyMeshFromSelf(&in);

    in.save_nodes("plc");
    in.save_poly("plc");
    tetrahedralize("pd", &in, &out);
    tetrahedralize("p", &in, &out, NULL); //pq1.414a0.1

    CopyMeshToSelf(&out);

}

}
