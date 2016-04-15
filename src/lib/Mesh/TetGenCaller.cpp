#include "TetGenCaller.h"

namespace voxel2tet
{

TetGenCaller::TetGenCaller()
{

}

void TetGenCaller :: UpdateVertexMapping()
{
    this->VertexMap.clear();

}

void TetGenCaller :: CopyMeshFromSelf(tetgenio *io)
{

    this->UpdateVertexMapping();

    io->initialize();

    // Copy mesh vertices
    io->numberofpoints = this->Mesh->Vertices.size();
    io->pointlist = new REAL[io->numberofpoints * 3];

    int cid=0; // Coordinate id. (0, 1, 2), (3, 4, 5), (6, 7, 8) ...

    for (VertexType *v: this->Mesh->Vertices) {
        for (double coordvalue: v->get_c()) {
            io->pointlist[cid]=coordvalue;
            cid++;
        }
    }

    // Copy mesh triangles
    io->numberoffacets = this->Mesh->Triangles.size();
    io->facetlist = new tetgenio::facet[io->numberoffacets]; // List of triangles
    io->facetmarkerlist = new int[io->numberoffacets]; // List of tag for each triangle

    int tid=0;
    for (TriangleType *t: this->Mesh->Triangles) {

        tetgenio::facet *f = &io->facetlist[tid];

        f->numberofpolygons = 1;
        f->holelist = NULL;
        f->numberofholes = 0;
        f->polygonlist = new tetgenio::polygon[f->numberofpolygons];

        tetgenio::polygon *p = &f->polygonlist[0];
        p->numberofvertices = 3;
        p->vertexlist = new int[p->numberofvertices];

        for (int i=0; i<3; i++) {
            VertexType *v = t->Vertices.at(i);
            p->vertexlist[i] = v->ID;
        }

        //io->facetmarkerlist[tid] = t->InterfaceID;

        tid++;

    }



}

void TetGenCaller :: CopyMeshToSelf(tetgenio *io)
{

}

void TetGenCaller :: EmbarasingTestExample()
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

    in.numberoffacets = 6;
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
    p->numberofvertices = 4;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0] = 4;
    p->vertexlist[1] = 8;
    p->vertexlist[2] = 5;
    p->vertexlist[3] = 1;

    // Set 'in.facetmarkerlist'

    in.facetmarkerlist[0] = -1;
    in.facetmarkerlist[1] = -2;
    in.facetmarkerlist[2] = 0;
    in.facetmarkerlist[3] = 0;
    in.facetmarkerlist[4] = 0;
    in.facetmarkerlist[5] = 0;

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

    tetgenio *in, *out;
    in = new tetgenio;
    in->firstnumber = 1;

    CopyMeshFromSelf(in);

    in->save_edges("/tmp/test");
    in->save_nodes("/tmp/test");
    in->save_faces("/tmp/test");
    in->save_poly("/tmp/test");


    tetrahedralize("q", in, out);
    free (in);



}

}
