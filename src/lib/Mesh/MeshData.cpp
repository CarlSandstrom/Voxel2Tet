#include <vector>
#include "MeshData.h"
#include "MiscFunctions.h"
#include "TetGenExporter.h"
#include "OFFExporter.h"
#include "OOFEMExporter.h"
#include "AbaqusExporter.h"

namespace voxel2tet
{
MeshData::MeshData(BoundingBoxType BoundingBox)
{
    this->BoundingBox = BoundingBox;
    this->VertexOctreeRoot = new VertexOctreeNode(this->BoundingBox, &this->Vertices, 0);
    this->TriangleCounter = 0;
    this->EdgeCounter = 0;
}

MeshData::~MeshData()
{
    for (unsigned int i = 0; i < this->Edges.size(); i++) {
        delete this->Edges.at(i);
    }

    for (auto t : this->Triangles) {
        delete t;
    }

    for (auto v : this->Vertices) {
        delete v;
    }

    for (auto t : this->Tets) {
        delete t;
    }

    delete this->VertexOctreeRoot;
}

void MeshData::DoSanityCheck()
{
    // Check edge for duplicates and make sure that each edge conatins two unique vertices
    for (size_t i = 0; i < this->Edges.size(); i++) {
        EdgeType *e1 = this->Edges[i];
        if (e1->Vertices[0] == e1->Vertices[1]) {
            LOG ("Edge is a point\n", 0);
            throw (0);
        }
        for (size_t j = i + 1; j < this->Edges.size(); j++) {
            EdgeType *e2 = this->Edges[j];
            if (((e1->Vertices[0] == e2->Vertices[0]) & (e1->Vertices[1] == e2->Vertices[1])) |
                ((e1->Vertices[0] == e2->Vertices[1]) & (e1->Vertices[1] == e2->Vertices[0]))) {
                LOG("Duplicate edge! %u@%p (%u, %u) and %u@%p (%u, %u)\n", e1->ID, e1, e1->Vertices[0]->ID,
                    e1->Vertices[1]->ID, e2->ID, e2, e2->Vertices[0]->ID, e2->Vertices[1]->ID);
                throw (0);
            }
        }
    }

    // Check list of edges and triangles and ensure that each edge at least holds two triangles
/*    bool EdgeMissing = false;
    for (TriangleType *t: this->Triangles) {
        for (int i=0; i<3; i++) {
            EdgeType te;
            te.Vertices[0] = t->Vertices[i];
            te.Vertices[1] = ( i<2 ? t->Vertices[i+1] : t->Vertices[0]);

            bool EdgeFound = false;
            for (EdgeType *e: this->Edges) {
                if ( ((e->Vertices[0]==te.Vertices[0]) & (e->Vertices[1]==te.Vertices[1])) | ((e->Vertices[0]==te.Vertices[1]) & (e->Vertices[1]==te.Vertices[0]))) {

                    std :: vector< TriangleType * > ConnectedTriangles = e->GiveTriangles();
                    if (ConnectedTriangles.size()<2) {
                        LOG("To few triangles\n", 0);
                        throw (0);
                    }

                    EdgeFound = true;
                    break;
                }
            }

            if (!EdgeFound) {
                EdgeMissing = true;
                LOG("\tTriangle %u, edge %u\n", t->ID, i);
            }

        }

    }

    if (EdgeMissing) {
        throw (0);
    }*/

    /*
    // Check if vertices in triangles are the same as vertices in edges
    std::vector<VertexType *> VerticesTriangles;
    std::vector<VertexType *> VerticesEdges;

    for (TriangleType *t: this->Triangles) {
        for (int i=0; i<3; i++) {
            VerticesTriangles.push_back( t->Vertices[i] );
        }
    }

    for (EdgeType *e: this->Edges) {
        for (int i=0; i<2; i++) {
            VerticesEdges.push_back(e->Vertices[i]);
        }
    }

    std :: sort( VerticesTriangles.begin(), VerticesTriangles.end() );
    VerticesTriangles.erase( std :: unique( VerticesTriangles.begin(), VerticesTriangles.end() ), VerticesTriangles.end() );

    std :: sort( VerticesEdges.begin(), VerticesEdges.end() );
    VerticesEdges.erase( std :: unique( VerticesEdges.begin(), VerticesEdges.end() ), VerticesEdges.end() );

    if (VerticesEdges.size() != VerticesTriangles.size()) {
        STATUS("Set of vertices in edges and set of vertices in triangles differ\n", 0);

        std::vector<VertexType *> InTriNotEdges;
        std::set_difference(VerticesTriangles.begin(), VerticesTriangles.end(), VerticesEdges.begin(), VerticesEdges.end(), std::inserter(InTriNotEdges, InTriNotEdges.end()));
        LOG("Vertices in triangles not on any edge:\n", 0);
        for (VertexType *v: InTriNotEdges) {LOG("\t%u\n", v->ID);}

        std::vector<VertexType *> OnEdgesNotTri;
        std::set_difference(VerticesEdges.begin(), VerticesEdges.end(), VerticesTriangles.begin(), VerticesTriangles.end(), std::inserter(OnEdgesNotTri, OnEdgesNotTri.end()));
        LOG("Vertices on edge not in any triangle:\n", 0);
        for (VertexType *v: OnEdgesNotTri) {LOG("\t%u\n", v->ID);}

        throw(0);
    } */

/*
    // Does the list of triangles vertices correspond to the list of vertices triangles?
    for ( TriangleType *t : this->Triangles ) {
        bool TriangleFound = false;
        for ( VertexType *v : t->Vertices ) {
            for ( TriangleType *vt : v->Triangles ) {
                if ( vt == t ) {
                    TriangleFound = true;
                }
            }
        }
        if ( !TriangleFound ) {
            STATUS("Triangles and vertices does not match!\n", 0);
            throw 0;
        }
        if ( t->ID > this->TriangleCounter ) {
            STATUS("Invalid triangle ID!\n", 0);
            throw 0;
        }
    }

    // Does the list of edge vertices correspond to the list of vertices triangles?
    for ( EdgeType *e : this->Edges ) {
        bool EdgeFound = false;
        for ( VertexType *v : e->Vertices ) {
            for ( EdgeType *ve : v->Edges ) {
                if ( ve == e ) {
                    EdgeFound = true;
                }
            }
        }
        if ( !EdgeFound ) {
            STATUS("Edges and vertices does not match!\n", 0);
            throw 0;
        }
    }

    // Check the same, but for vertices
    for ( VertexType *v : this->Vertices ) {
        bool VertexFoundInTriangle = ( v->Triangles.size() == 0 ) ? true : false;

        for ( TriangleType *t : v->Triangles ) {
            for ( VertexType *vt : t->Vertices ) {
                if ( vt == v ) {
                    VertexFoundInTriangle = true;
                }
            }
        }

        if ( !VertexFoundInTriangle ) {
            STATUS("Triangle does not contain vertex while vertex contains triangle\n", 0);
            throw 0;
        }

        bool VertexFoundInEdge = ( v->Edges.size() == 0 ) ? true : false;

        for ( EdgeType *e : v->Edges ) {
            for ( VertexType *ve : e->Vertices ) {
                if ( ve == v ) {
                    VertexFoundInEdge = true;
                }
            }
        }

        if ( !VertexFoundInEdge ) {
            STATUS("Vertex not found in edge list while edge is found in vertex list of edges\n", 0);
            throw 0;
        }
    }


    // Are all triangles unique?
    for ( TriangleType *t1 : this->Triangles ) {
        for ( TriangleType *t2 : this->Triangles ) {
            if ( t1 != t2 ) {
                bool ispermutation = std :: is_permutation( t1->Vertices.begin(), t1->Vertices.end(), t2->Vertices.begin() );
                if ( ispermutation ) {
                    STATUS("\nDuplicate triangles\n", 0);
                    throw 0;
                }
            }
        }
    }*/

    LOG("DoSanityCheck finished\n", 0);

}

void MeshData::ExportSurface(std::string FileName, Exporter_FileTypes FileType, int VolumeID)
{
    STATUS("Export surface to %s\n", FileName.c_str());
    Exporter *exporter;
    switch (FileType) {
        case FT_OFF: {
            exporter = new OFFExporter(&this->Triangles, &this->Vertices, &this->Edges, &this->Tets);
            break;
        }
        case FT_Poly: {
            exporter = new TetGenExporter(&this->Triangles, &this->Vertices, &this->Edges, &this->Tets);
            break;
        }
        case FT_VTK: {
            exporter = new VTKExporter(&this->Triangles, &this->Vertices, &this->Edges, &this->Tets);
            break;
        }
        case FT_SIMPLE: {
            exporter = new SimpleExporter(&this->Triangles, &this->Vertices, &this->Edges, &this->Tets);
            break;
        }
        case FT_STL: {
            exporter = new STLExporter(&this->Triangles, &this->Vertices, &this->Edges, &this->Tets);
            break;
        }
        default: {
            LOG("File type not supported!\n", 0);
            throw (0);
        }
    }
    exporter->WriteSurfaceData(FileName, VolumeID);
    free(exporter);
}

void MeshData::ExportVolume(std::string FileName, Exporter_FileTypes FileType)
{
    STATUS("Export volume to %s\n", FileName.c_str());
    Exporter *exporter;
    switch (FileType) {
        case FT_OFF: {
            exporter = new OFFExporter(&this->Triangles, &this->Vertices, &this->Edges, &this->Tets);
            break;
        }
        case FT_Poly: {
            exporter = new TetGenExporter(&this->Triangles, &this->Vertices, &this->Edges, &this->Tets);
            break;
        }
        case FT_VTK: {
            exporter = new VTKExporter(&this->Triangles, &this->Vertices, &this->Edges, &this->Tets);
            break;
        }
        case FT_OOFEM: {
            exporter = new OOFEMExporter(&this->Triangles, &this->Vertices, &this->Edges, &this->Tets);
            break;
        }
        case FT_ABAQUS: {
            exporter = new AbaqusExporter(&this->Triangles, &this->Vertices, &this->Edges, &this->Tets, false);
            break;
        }
        case FT_ABAQUSPHON: {
            exporter = new AbaqusExporter(&this->Triangles, &this->Vertices, &this->Edges, &this->Tets, true);
            break;
        }
        default:
            STATUS("No export filter found for %s\n", FileName.c_str());
            return;
    }
    exporter->WriteVolumeData(FileName);
    free(exporter);
}

EdgeType *MeshData::AddEdge(std::array<int, 2> VertexIDs)
{
    EdgeType *NewEdge = new EdgeType;
    for (unsigned int i : {0, 1}) {
        NewEdge->Vertices.at(i) = this->Vertices.at(VertexIDs.at(i));
    }

    return AddEdge(NewEdge);
}

EdgeType *MeshData::AddEdge(EdgeType *e)
{
    LOG("Create edge from Vertex IDs %u and %u: %p\n", e->Vertices[0]->ID, e->Vertices[1]->ID, e);

    // Check if edge already exists (only check edges connected to the first vertex)
    VertexType *ThisVertex = e->Vertices[0];
    VertexType *OtherVertex = e->Vertices[1];

    std::vector<EdgeType *> Edges = e->Vertices.at(0)->Edges;

    for (auto Edge : Edges) {
        if (((Edge->Vertices[0] == ThisVertex) & (Edge->Vertices[1] == OtherVertex)) |
            ((Edge->Vertices[1] == ThisVertex) & (Edge->Vertices[0] == OtherVertex))) {
            LOG("Edge %u@%p already exists\n", Edge->ID, Edge);
            return Edge;
        }
    }

    if (e->Vertices[0] == e->Vertices[1]) {
        STATUS("Vertices are the same\n", 0);
        throw (0);
    }

    e->ID = this->EdgeCounter;
    this->EdgeCounter++;
    for (VertexType *v : e->Vertices) {
        v->AddEdge(e);
    }
    LOG("Add edge %u@%p from vertices (%u, %u)\n", e->ID, e, e->Vertices[0]->ID, e->Vertices[1]->ID);
    this->Edges.push_back(e);
    return e;
}

void MeshData::RemoveEdge(EdgeType *e)
{
    LOG("Remove edge %u@%p (%u, %u)\n", e->ID, e, e->Vertices[0]->ID, e->Vertices[1]->ID);
    for (VertexType *v : e->Vertices) {
        v->RemoveEdge(e);
    }
    this->Edges.erase(std::remove(this->Edges.begin(), this->Edges.end(), e), this->Edges.end());
    delete e;
    e = NULL;
}

void MeshData::RemoveTriangle(TriangleType *t)
{
    LOG("Remove triangle %u\n", t->ID);
    for (VertexType *v : t->Vertices) {
        v->RemoveTriangle(t);
    }
    this->Triangles.erase(std::remove(this->Triangles.begin(), this->Triangles.end(), t), this->Triangles.end());
    delete t;
    t = NULL;
}

TriangleType *MeshData::AddTriangle(std::array<double, 3> v0, std::array<double, 3> v1, std::array<double, 3> v2)
{
    // Insert vertices and create a triangle using the indices returned
    int VertexIDs[3];

    VertexIDs[0] = this->VertexOctreeRoot->AddVertex(v0[0], v0[1], v0[2]);
    VertexIDs[1] = this->VertexOctreeRoot->AddVertex(v1[0], v1[1], v1[2]);
    VertexIDs[2] = this->VertexOctreeRoot->AddVertex(v2[0], v2[1], v2[2]);

    return this->AddTriangle({VertexIDs[0], VertexIDs[1], VertexIDs[2]});
}

TriangleType *MeshData::AddTriangle(std::array<int, 3> VertexIDs)
{
    TriangleType *NewTriangle = new TriangleType;
    LOG("Create triangle %p from vertices (%u, %u, %u)@(%p, %p, %p)\n", NewTriangle, VertexIDs.at(0),
        VertexIDs.at(1), VertexIDs.at(2),
        this->Vertices.at(VertexIDs.at(0)), this->Vertices.at(VertexIDs.at(1)), this->Vertices.at(VertexIDs.at(2)));

    for (int i = 0; i < 3; i++) {
        if (i < 2) {
            this->AddEdge({{VertexIDs.at(i), VertexIDs.at(i + 1)}});
        } else {
            this->AddEdge({{VertexIDs.at(i), VertexIDs.at(0)}});
        }
    }

    // Update vertices and triangles with references to each other
    for (int i = 0; i < 3; i++) {
        this->Vertices.at(VertexIDs.at(i))->AddTriangle(NewTriangle);
        NewTriangle->Vertices[i] = this->Vertices.at(VertexIDs.at(i));
    }

    return AddTriangle(NewTriangle);
}

TriangleType *MeshData::AddTriangle(TriangleType *NewTriangle)
{
#if SANITYCHECK == 1
    int i = 0;
    for ( TriangleType *t : this->Triangles ) {
        bool permutation = std :: is_permutation( t->Vertices.begin(), t->Vertices.end(), NewTriangle->Vertices.begin() );
        if ( permutation ) {
            STATUS("\nTriangle already exist. Existing ID = %u (index %u in list)!\n", t->ID, i); //TODO: Add a logging command for errors
            return t;

            throw 0;
        }
        i++;
    }
#endif

    NewTriangle->UpdateNormal();
    NewTriangle->ID = TriangleCounter;
    LOG("Add triangle %u to set\n", NewTriangle->ID);

    for (VertexType *v : NewTriangle->Vertices) {
        v->AddTriangle(NewTriangle);
    }

    TriangleCounter++;
    this->Triangles.push_back(NewTriangle);
    return NewTriangle;
}

TetType *MeshData::AddTetrahedron(std::array<int, 4> VertexIDs)
{
    TetType *NewTet = new TetType;
    NewTet->Vertices = {{this->Vertices.at(VertexIDs[0]), this->Vertices.at(VertexIDs[1]), this->Vertices.at(
            VertexIDs[2]), this->Vertices.at(VertexIDs[3])}};
    return this->AddTetrahedron(NewTet);
}

TetType *MeshData::AddTetrahedron(TetType *NewTet)
{
    // TODO: Ensure that all elements in Tets are unique
    NewTet->ID = this->Tets.size();
    this->Tets.push_back(NewTet);
    return NewTet;
}

void MeshData::RemoveTetrahedron(TetType *t)
{
    throw (0); // TODO: Complete this
}

std::vector<TriangleType *> MeshData::GetTrianglesAround(std::array<double, 3> c, double r)
{
    std::vector<VertexType *> NearVertices = this->VertexOctreeRoot->GiveVerticesWithinSphere(c[0], c[1], c[2], r);
    std::sort(NearVertices.begin(), NearVertices.end(), SortByID<VertexType *>);
    NearVertices.erase(std::unique(NearVertices.begin(), NearVertices.end()), NearVertices.end());

    std::vector<TriangleType *> NearTriangles;
    for (VertexType *v : NearVertices) {
        for (TriangleType *t : v->Triangles) {
            NearTriangles.push_back(t);
        }
    }
    std::sort(NearTriangles.begin(), NearTriangles.end(), SortByID<TriangleType *>);
    NearTriangles.erase(std::unique(NearTriangles.begin(), NearTriangles.end()), NearTriangles.end());

    return NearTriangles;
}

bool MeshData::CheckSameOrientation(TriangleType *t1, TriangleType *t2)
{
    for (int i = 0; i < 3; i++) {
        int ThisIndext1 = i;
        int NextIndext1 = i == 2 ? 0 : i + 1;
        int PrevIndext1 = i == 0 ? 2 : i - 1;
        for (int j = 0; j < 3; j++) {
            int ThisIndext2 = j;
            int NextIndext2 = j == 2 ? 0 : j + 1;
            int PrevIndext2 = j == 0 ? 2 : j - 1;
            if (t1->Vertices[ThisIndext1] == t2->Vertices[ThisIndext2]) {
                if ((t1->Vertices[NextIndext1] == t2->Vertices[PrevIndext2]) |
                    (t1->Vertices[PrevIndext1] == t2->Vertices[NextIndext2])) {
                    return true;
                } else {
                    return false;
                }
            }
        }
    }
    return false;
}

FC_MESH MeshData::CheckTrianglePenetration(TriangleType *t1, TriangleType *t2)
{
    int sharedvertices;
    FC_MESH Result = CheckTrianglePenetration(t1->Vertices, t2->Vertices, sharedvertices);
    if (Result == FC_DUPLICATETRIANGLE) {
        // It is naturally ok if a newly created, not yet used, triangle (negative ID) shared vertices with an existing triangle
        if (t1->ID == -t2->ID) {
            Result = FC_OK;
        }
    }
    return Result;
}

FC_MESH MeshData::CheckTrianglePenetration(std::array<VertexType *, 3> t1, std::array<VertexType *, 3> t2,
                                           int &sharedvertices)
{

    //LOG("Check for penetration of triangles %u [%u, %u, %u] and %u [%u, %u, %u]\n", t1[0]->ID, t1[1]->ID, t1[2]->ID, t2[0]->ID, t2[1]->ID, t2[2]->ID);

    sharedvertices = 0;
    std::vector<VertexType *> SharedVerticesList;
    for (VertexType *v1 : t1) {
        for (VertexType *v2 : t2) {
            if (v1 == v2) {
                SharedVerticesList.push_back(v1);
                sharedvertices++;
            }
        }
    }

    if (sharedvertices == 0) {
        // I am a bit unsure if this is correct. Naturally two triangles with one shared vertex can intersect, but will we ever have that situation?
        // Anyway, best would be to check if edges not members of the other triangle penetrates the surface since there is a "bug" in the used algorithm
        // that gives a "false" true if vertices are shared.

        double V0[3] = {
                t1[0]->get_c(0), t1[0]->get_c(1), t1[0]->get_c(2)
        };
        double V1[3] = {
                t1[1]->get_c(0), t1[1]->get_c(1), t1[1]->get_c(2)
        };
        double V2[3] = {
                t1[2]->get_c(0), t1[2]->get_c(1), t1[2]->get_c(2)
        };

        double U0[3] = {
                t2[0]->get_c(0), t2[0]->get_c(1), t2[0]->get_c(2)
        };
        double U1[3] = {
                t2[1]->get_c(0), t2[1]->get_c(1), t2[1]->get_c(2)
        };
        double U2[3] = {
                t2[2]->get_c(0), t2[2]->get_c(1), t2[2]->get_c(2)
        };

        int intersects = tri_tri_intersect(V0, V1, V2, U0, U1, U2);

        if (intersects == 1) {
            return FC_TRIANGLESINTERSECT;
        }
    } else if (sharedvertices == 1) {
        VertexType *SharedVertex = SharedVerticesList[0];
        for (int i = 0; i < 2; i++) { // Edge index
            std::array<VertexType *, 3> OtherTriangle = (i == 0) ? t2 : t1;
            std::array<VertexType *, 3> ThisTriangle = (i == 0) ? t1 : t2;  // Triangle being checked if penetrated

            std::array<VertexType *, 2> FreeEdge;
            int k = 0;
            for (int j = 0; j < 3; j++) {
                if (OtherTriangle[j] != SharedVertex) {
                    FreeEdge[k] = OtherTriangle[j];
                    k++;
                }
            }

            // Check if the free edge penetrates ThisTriangle
            double V0[3] = {
                    ThisTriangle[0]->get_c(0), ThisTriangle[0]->get_c(1), ThisTriangle[0]->get_c(2)
            };
            double V1[3] = {
                    ThisTriangle[1]->get_c(0), ThisTriangle[1]->get_c(1), ThisTriangle[1]->get_c(2)
            };
            double V2[3] = {
                    ThisTriangle[2]->get_c(0), ThisTriangle[2]->get_c(1), ThisTriangle[2]->get_c(2)
            };

            double O[3] = {
                    FreeEdge[0]->get_c(0), FreeEdge[0]->get_c(1), FreeEdge[0]->get_c(2)
            };
            double D[3] = {
                    FreeEdge[1]->get_c(0) - O[0], FreeEdge[1]->get_c(1) - O[1], FreeEdge[1]->get_c(2) - O[2]
            };
            float alpha;

            if (triangle_ray_intersection(V0, V1, V2, O, D, &alpha)) {
                if ((alpha >= 0.0) & (alpha <= 1.0)) {
                    return FC_TRIANGLESINTERSECT;
                }
            } else {
                // Check if points are located inside triangle (ray could lie in the same plane as triangle)
                for (VertexType *v: FreeEdge) {
                    double P[3] = {v->get_c(0), v->get_c(1), v->get_c(2)};
                    if (point_in_tri(V0, V1, V2, P)) {
                        return FC_TRIANGLESINTERSECT;
                    }
                }
            }
        }
    } else if (sharedvertices == 2) {
        // If two triangles share an edge, check if each remaining (not shared) vertex is located within the other triangle
        std::vector<VertexType *> UniqueVertices;
        for (std::array<VertexType *, 3> list : {t1, t2}) {
            for (VertexType *v : list) {
                if ((v != SharedVerticesList[0]) & (v != SharedVerticesList[1])) {
                    UniqueVertices.push_back(v);
                }
            }
        }

        double s0[3] = {
                SharedVerticesList[0]->get_c(0), SharedVerticesList[0]->get_c(1), SharedVerticesList[0]->get_c(2)
        };
        double s1[3] = {
                SharedVerticesList[1]->get_c(0), SharedVerticesList[1]->get_c(1), SharedVerticesList[1]->get_c(2)
        };

        double u0[3] = {
                UniqueVertices[0]->get_c(0), UniqueVertices[0]->get_c(1), UniqueVertices[0]->get_c(2)
        };
        double u1[3] = {
                UniqueVertices[1]->get_c(0), UniqueVertices[1]->get_c(1), UniqueVertices[1]->get_c(2)
        };

        if (tri_tri_intersect_shared_edge(s0, s1, u0, u1)) {
            return FC_TRIANGLESINTERSECT;
        }
    } else if (sharedvertices == 3) {
        return FC_DUPLICATETRIANGLE;
    }
    return FC_OK;
}
}
