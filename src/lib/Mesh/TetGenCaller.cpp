#include "TetGenCaller.h"

namespace voxel2tet
{
TetGenCaller :: TetGenCaller()
{}

void TetGenCaller :: UpdateVertexMapping()
{
    this->VertexMapFromID.clear();
    this->VertexMapFromTetgen.clear();

    std :: vector< VertexType * >VertexIDsFromTriangles;

    for ( TriangleType *t : this->Mesh->Triangles ) {
        for ( VertexType *v : t->Vertices ) {
            VertexIDsFromTriangles.push_back(v);
        }
    }

    std :: sort( VertexIDsFromTriangles.begin(), VertexIDsFromTriangles.end() );

    VertexIDsFromTriangles.erase( std :: unique( VertexIDsFromTriangles.begin(), VertexIDsFromTriangles.end() ), VertexIDsFromTriangles.end() );

    for ( unsigned int i = 0; i < VertexIDsFromTriangles.size(); i++ ) {
        this->VertexMapFromID [ VertexIDsFromTriangles.at(i) ] = i;
        this->VertexMapFromTetgen [ i ] = VertexIDsFromTriangles.at(i);
    }
}

void TetGenCaller :: CopyMeshFromSelf(tetgenio *in)
{
    this->UpdateVertexMapping();
    tetgenio :: facet *f;
    tetgenio :: polygon *p;

    in->initialize();
    in->firstnumber = 0;

    // Copy mesh vertices
    in->numberofpoints = this->VertexMapFromID.size();
    in->pointlist = new REAL [ in->numberofpoints * 3 ];

    unsigned int cid = 0; // Coordinate id. (0, 1, 2), (3, 4, 5), (6, 7, 8) ...

    std :: map< VertexType *, int > :: iterator iter;
    for ( iter = this->VertexMapFromID.begin(); iter != this->VertexMapFromID.end(); iter++ ) {
        VertexType *v = iter->first;
        for ( double coordvalue : v->get_c() ) {
            in->pointlist [ cid ] = coordvalue;
            cid++;
        }
    }

    // Copy mesh triangles
    in->numberoffacets = this->Mesh->Triangles.size();
    in->facetlist = new tetgenio :: facet [ in->numberoffacets ]; // List of triangles
    in->facetmarkerlist = new int [ in->numberoffacets ]; // List of tag for each triangle

    unsigned int tid = 0;
    for ( TriangleType *t : this->Mesh->Triangles ) {
        f = & in->facetlist [ tid ];

        f->numberofpolygons = 1;
        f->holelist = ( REAL * ) NULL;
        f->numberofholes = 0;
        f->polygonlist = new tetgenio :: polygon [ f->numberofpolygons ];

        p = & f->polygonlist [ 0 ];
        p->numberofvertices = 3;
        p->vertexlist = new int [ p->numberofvertices ];

        for ( int i = 0; i < 3; i++ ) {
            VertexType *v = t->Vertices.at(i);
            int TetGenID = this->VertexMapFromID [ v ];
            p->vertexlist [ i ] = TetGenID;
        }

        in->facetmarkerlist [ tid ] = t->InterfaceID;

        tid++;
    }
}

MeshData *TetGenCaller :: CopyTetMesh(tetgenio *io)
{
    MeshData *NewMesh = new MeshData(this->Mesh->BoundingBox);

    // Add vertices
    for ( int i = 0; i < io->numberofpoints; i++ ) {
        double *c;
        c = & io->pointlist [ 3 * i ];
        NewMesh->VertexOctreeRoot->AddVertex(c [ 0 ], c [ 1 ], c [ 2 ]);
    }

    // Add triangles
    for ( int i = 0; i < io->numberoftrifaces; i++ ) {
        int *triface;
        int marker = io->trifacemarkerlist [ i ];
        triface = & io->trifacelist [ 3 * i ];
        TriangleType *t = NewMesh->AddTriangle({ triface [ 0 ], triface [ 1 ], triface [ 2 ] });
        t->InterfaceID = marker;
    }

    // Add tetrahedrons
    for ( int i = 0; i < io->numberoftetrahedra; i++ ) {
        int *tet = & io->tetrahedronlist [ 4 * i ];
        TetType *t = NewMesh->AddTetrahedron({ tet [ 0 ], tet [ 1 ], tet [ 2 ], tet [ 3 ] });

        if ( io->numberoftetrahedronattributes == 1 ) {
            int tetattr = io->tetrahedronattributelist [ i ];
            t->MaterialID = tetattr;
        } else {
            t->MaterialID = 0;
        }
    }

    return NewMesh;
}

MeshData *TetGenCaller :: Execute()
{
    STATUS("\n== Create tet-mesh with TetGen ===================================\n", 0);

    tetgenio in, out;

    CopyMeshFromSelf(& in);

    tetrahedralize( ( char * ) "pq2/10AV", & in, & out, NULL ); //pq1.414a0.1
    STATUS("==================================================================\n\n", 0);

    return CopyTetMesh(& out);
}

void TetGenCaller :: TestMesh()
{
    STATUS("\n== Test mesh with TetGen =========================================\n", 0);

    tetgenio in, out;

    CopyMeshFromSelf(& in);

    tetrahedralize( ( char * ) "pd", & in, & out );

    if (out.numberoftrifaces>0) {
        STATUS("Check failed\n", 0);

        for ( int i = 0; i < out.numberoftrifaces; i++ ) {
            int marker = out.trifacemarkerlist [ i ];
            STATUS ("\t%u\n", marker);
        }

        throw(0);
    }

    STATUS("==================================================================\n\n", 0);
}
}
