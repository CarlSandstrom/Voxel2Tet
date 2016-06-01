#include <algorithm>
#include <vector>
#include <iterator>
#include <time.h>
#include <omp.h>

#include "Voxel2Tet.h"
#include "Importer.h"
#include "Dream3DDataReader.h"
#include "CallbackImporter.h"
#include "TetGenCaller.h"
#include "Smoother.h"

namespace voxel2tet
{

Voxel2Tet::Voxel2Tet(Options *Opt)
{
    this->Opt = Opt;

    // Set defult options
    this->Opt->AddDefaultMap("spring_c_factor", "0.5");
    this->Opt->AddDefaultMap("edge_spring_c_factor", "0.5");

    LOG("Starting Voxel2Tet\n", 0);
}

Voxel2Tet::~Voxel2Tet()
{
    for (unsigned int i=0; i<PhaseEdges.size(); i++) {
        delete this->PhaseEdges.at(i);
    }
    
    for (auto s: this->Surfaces) {
        delete s;
    }
    
    delete this->Mesh;
    
    //delete this->Imp;
    
    //for (auto p: this->PhaseEdges) delete p;
    
}

void Voxel2Tet::LoadCallback(cbMaterialIDByCoordinate MaterialIDByCoordinate, std::array<double,3> origin, std::array<double,3> spacing, std::array<int,3> dimensions)
{
    STATUS ("Setup callback functions\n", 0);
    CallbackImporter *DataReader = new CallbackImporter(MaterialIDByCoordinate, origin, spacing, dimensions);
    this->Imp = DataReader;
    FinalizeLoad();
}

void Voxel2Tet::LoadFile(std::string Filename)
{
    STATUS ("Load file %s\n", Filename.c_str());
    Dream3DDataReader *DataReader = new Dream3DDataReader(this->Opt->GiveStringValue("DataContainer"), this->Opt->GiveStringValue("MaterialId") );
    DataReader->LoadFile(Filename);
    this->Imp = DataReader;
    FinalizeLoad();
    
}

void Voxel2Tet :: FinalizeLoad()
{
    
    double cellspace[3];
    int dim[3];
    this->Imp->GiveSpacing(cellspace);
    this->Imp->GiveDimensions(dim);
    
    STATUS("\tVoxel dimensions are %f * %f * %f\n", cellspace[0],cellspace[1],cellspace[2]);
    STATUS("\tNumber of voxels:%u\n", dim[0]*dim[1]*dim[2]);

    auto_c = false;
    
    if (this->Opt->has_key("spring_alpha")) {
        spring_alpha = Opt->GiveDoubleValue("spring_alpha");
    } else {
        spring_alpha = 4;
    }
    
    if (this->Opt->has_key("spring_c")) {
        spring_c = Opt->GiveDoubleValue("spring_c");
    } else {
        spring_c = Compute_c(cellspace[0]*Opt->GiveDoubleValue("spring_c_factor"), spring_alpha);
    }
    
    if (this->Opt->has_key("edge_spring_alpha")) {
        edgespring_alpha = Opt->GiveDoubleValue("edge_spring_alpha");
    } else {
        edgespring_alpha = 4;
    }
    
    if (this->Opt->has_key("edge_spring_c")) {
        edgespring_c = Opt->GiveDoubleValue("edge_spring_c");
    } else {
        edgespring_c = Compute_c(cellspace[0]*Opt->GiveDoubleValue("edge_spring_c_factor"), edgespring_alpha);
    }
    
    STATUS("Using spring_alpha=%f, spring_c=%f\n", spring_alpha, spring_c);
    STATUS("Using edgespring_alpha=%f, edgespring_c=%f\n", edgespring_alpha, edgespring_c);
    
    BoundingBoxType bb;
    double spacing[3];
    
    bb = this->Imp->GiveBoundingBox();
    this->Imp->GiveSpacing(spacing);
    
    for (int i=0; i<3; i++) {
        bb.maxvalues[i] = bb.maxvalues[i] + spacing[i];
        bb.minvalues[i] = bb.minvalues[i] - spacing[i];
    }
    Mesh = new MeshManipulations(bb);
}

void Voxel2Tet::LoadData()
{
    STATUS ("Load data\n", 0);
    
    if (this->Opt->has_key("i")) {
        LoadFile(this->Opt->GiveStringValue("i"));
    }
    
}

void Voxel2Tet :: UpdateSurfaces()
{
    for (Surface *s: this->Surfaces) {
        s->Triangles.clear();
    }
    
    for (TriangleType *t: this->Mesh->Triangles) {
        this->Surfaces.at(t->InterfaceID)->Triangles.push_back(t);
    }
}

void Voxel2Tet::ExportSurface(std::string FileName, Exporter_FileTypes FileType)
{
    this->Mesh->ExportSurface(FileName, FileType);
}

void Voxel2Tet :: Tetrahedralize()
{
    TetGenCaller Generator;
    Generator.Mesh = this->Mesh;
    Generator.TestMesh();

    MeshData *Mesh = Generator.Execute();

    MeshManipulations *NewMesh;
    NewMesh = static_cast <MeshManipulations*> (Mesh);

    // Update tetrahedrons with material information

    // Create mapping
    std::map<int, int> Tetgen2Self; // first: Tetgen material ID, second: Phase
    std::map<int, int>::iterator MapIter;

    for (TetType *t: NewMesh->Tets) {

        // Check if tetgen ID is already handled
        MapIter = Tetgen2Self.find(t->MaterialID);

        if (MapIter == Tetgen2Self.end()) { // Material not yet mapped
            std::array<double, 3> cm = t->GiveCenterOfMass();
            Volume *v = this->FindVolumeContainingPoint(cm);
            Tetgen2Self[t->MaterialID] = v->Phase;
        }

        t->MaterialID = Tetgen2Self[t->MaterialID];

    }

    NewMesh->ExportVolume("/tmp/TetVolume0.vtu", FT_VTK);
    NewMesh->ExportVolume("/tmp/TetVolume1.vtu", FT_VTK);
    NewMesh->ExportVolume("/tmp/TetVolume.in", FT_OOFEM);

}

void Voxel2Tet :: ExportVolume(std::string FileName, Exporter_FileTypes FileType)
{
    this->Mesh->ExportVolume(FileName, FileType);
}

void Voxel2Tet::FindSurfaces()
{
    
    STATUS("Find surfaces\n", 0);
    
    int dim[3];
    double spacing[3], origin[3];
    this->Imp->GiveDimensions(dim);
    this->Imp->GiveSpacing(spacing);
    this->Imp->GiveOrigin(origin);

    STATUS("Total volume: %f\n", dim[0]*spacing[0] * dim[1]*spacing[1] * dim[2]*spacing[2]);
    
    std::vector <double> signs = {1, -1};
    
    
    for (int i=0; i<dim[0]; i++) {
        for (int j=0; j<dim[1]; j++) {
            for (int k=0; k<dim[2]; k++) {
                // Specifies directions in which we will look for different materials
                std::vector <std::vector <double>> testdirections = {{1,0,0}, {0,1,0}, {0,0,1}};
                // If we the adjacent material is of other type, we will create a square by varying the coordinates marked '1' in vdirections
                std::vector <std::vector <double>> vdirections = {{0,1,1}, {1,0,1}, {1,1,0}};
                // vindex is the indices of the coordinates
                std::vector <std::vector <int>> vindex = {{1,2},{0,2},{0,1}};
                
                // If we are on a boundary, we need to check what is outside of that boundary
                if (i==0) {
                    testdirections.push_back({-1,0,0});
                    vdirections.push_back({0,1,1});
                    vindex.push_back({1,2});
                }
                if (j==0) {
                    testdirections.push_back({0,-1,0});
                    vdirections.push_back({1,0,1});
                    vindex.push_back({0,2});
                }
                if (k==0) {
                    testdirections.push_back({0,0,-1});
                    vdirections.push_back({1,1,0});
                    vindex.push_back({0,1});
                }
                
                int ThisPhase = this->Imp->GiveMaterialIDByIndex(i, j, k);
                int NeighboringPhase;
                bool SamePhase;
                
                // Check material in each direction
                
                for (unsigned int m=0; m<testdirections.size(); m++) {
                    int testi = testdirections.at(m).at(0) + i;
                    int testj = testdirections.at(m).at(1) + j;
                    int testk = testdirections.at(m).at(2) + k;
                    
                    // If comparing inside the domain, simply compare
                    if ( (testi>=0) & (testj>=0) & (testk>=0) & (testi<dim[0]) & (testj<dim[1]) & (testk<dim[2])) {
                        NeighboringPhase = this->Imp->GiveMaterialIDByIndex(testi, testj, testk);
                        SamePhase = (ThisPhase == NeighboringPhase);
                    } else {
                        // If we are comparing with the outside, take into account that a we might have void (i.e. 0) in both voxels
                        if (ThisPhase!=0) {
                            SamePhase = false;
                            NeighboringPhase = this->Imp->GiveMaterialIDByIndex(testi, testj, testk);
                        } else {
                            // Void-to-void connection
                            SamePhase = true;
                            NeighboringPhase = ThisPhase;
                        }
                    }
                    
                    if (!SamePhase) {   // Add surface. I.e. a square separating the two voxels
                        
                        // Compute centre off square
                        double c[3];
                        c[0] = (double(i) + double(testdirections.at(m).at(0))/2.0) * spacing[0] + origin[0] + spacing[0]/2.0;
                        c[1] = (double(j) + double(testdirections.at(m).at(1))/2.0) * spacing[1] + origin[1] + spacing[1]/2.0;
                        c[2] = (double(k) + double(testdirections.at(m).at(2))/2.0) * spacing[2] + origin[2] + spacing[2]/2.0;
                        
                        // Compute coordinate of corner point
                        double delta[3];
                        delta[0] = spacing[0]*vdirections.at(m)[0] / 2.0;
                        delta[1] = spacing[1]*vdirections.at(m)[1] / 2.0;
                        delta[2] = spacing[2]*vdirections.at(m)[2] / 2.0;
                        
                        std::vector <int> VoxelIDs;
                        
                        
                        for (auto s1: signs) {
                            for (auto s2: signs) {
                                double newvertex[3];
                                newvertex[0] = c[0];
                                newvertex[1] = c[1];
                                newvertex[2] = c[2];
                                
                                newvertex[vindex.at(m).at(0)] = newvertex[vindex.at(m).at(0)] + s1*delta[vindex.at(m)[0]];
                                newvertex[vindex.at(m).at(1)] = newvertex[vindex.at(m).at(1)] + s2*delta[vindex.at(m)[1]];
                                
                                int id = Mesh->VertexOctreeRoot->AddVertex(newvertex[0], newvertex[1], newvertex[2]);
                                LOG ("Corner (id=%u) at (%f, %f, %f)\n", id, newvertex[0], newvertex[1], newvertex[2]);
                                VoxelIDs.push_back(id);
                            }
                        }
                        AddSurfaceSquare(VoxelIDs, {ThisPhase, NeighboringPhase}, NeighboringPhase);
                        
                    }
                }
            }
        }
    }

    this->UpdateSurfaces();


    STATUS("Find volumes\n", 0);
    for (Surface *s: this->Surfaces) {
        for (int p: s->Phases) {
            if (p>=0) {
                bool PhaseFound = false;
                for (Volume *v: this->Volumes) {
                    if (v->Phase==p) {
                        v->Surfaces.push_back(s);
                        PhaseFound = true;
                        break;
                    }
                }
                if (!PhaseFound) {
                    Volume *NewVolume = new Volume(p);
                    this->Volumes.push_back(NewVolume);

                    NewVolume->Surfaces.push_back(s);
                }
            }
        }
    }
}


void Voxel2Tet :: FindEdges()
{
    
    STATUS ("Find edges\n\tIdentify vertices shared by surfaces...\n", 0);
    
    // Sort vertex vectors on all surfaces
    for(auto surface: this->Surfaces) {
        std::sort(surface->Vertices.begin(), surface->Vertices.end());
    }
    
    std::vector <VertexType*> EdgeVertices;
    
    // Find all vertices that are shared among the surfaces. Since this is in 3D,
    // one vertex can be shared by several surfaces while not being an endpoint
    // of the edge. Thus, first find all shared vertices and then trace along x, y, z
    // to find the edge
    
    for (unsigned int s1=0; s1<this->Surfaces.size(); s1++) {
        Surface *surface1 = this->Surfaces.at(s1);
        std::vector <VertexType*> SurfaceEdgeVertices;
        for (unsigned int s2=0; s2<this->Surfaces.size(); s2++) {
            Surface *surface2 = this->Surfaces.at(s2);
            
            for (int i=0; i<2; i++) {
                if (s1==s2) break;
                bool SharedPhaseFound=false;
                std::vector <VertexType*> SharedVertices;
                int mat1 = surface1->Phases[i];
                for (int j=0; j<2; j++) {
                    int mat2 = surface2->Phases[j];
                    if (mat1==mat2) {
                        
                        LOG("Surfaces %p and %p shares phase %i\n", surface1, surface2, mat1);
                        
                        SharedPhaseFound = true;
                        
                        std::set_intersection(surface1->Vertices.begin(), surface1->Vertices.end(),
                                              surface2->Vertices.begin(), surface2->Vertices.end(),
                                              back_inserter(SharedVertices));
                        
                        // If both sets of vertices has an intersection, that intersection is an edge
                        if (SharedVertices.size()>0) {
                            LOG ("Surfaces intersects\n",0);
                            for (auto Vertex: SharedVertices) {
                                EdgeVertices.push_back(Vertex);
                                SurfaceEdgeVertices.push_back(Vertex);
                            }
                        }
                        break;
                    }
                }
                if (SharedPhaseFound) break;
            }
            
        }
        std::sort(SurfaceEdgeVertices.begin(), SurfaceEdgeVertices.end());
        SurfaceEdgeVertices.erase( std::unique(SurfaceEdgeVertices.begin(), SurfaceEdgeVertices.end()), SurfaceEdgeVertices.end());
        
        for (auto v: SurfaceEdgeVertices) {
            v->Fixed = {true, true, true};
        }
        
    }
    
    STATUS ("\tTrace edges...\n", 0);
    
    std::sort(EdgeVertices.begin(), EdgeVertices.end());
    EdgeVertices.erase( std::unique(EdgeVertices.begin(), EdgeVertices.end()), EdgeVertices.end());
    
    // Trace edges by for each vertex in EdgeVertices:
    //  * Check if the neightbour is in EdgeVertices as well
    //    * If so, add the edge segment to edge separating the materials surounding the midpoint of edge segment
    
    // Used to check which material surround the midpoint. See FindSurfaces for explanation.
    std::vector <std::vector <double>> testdirections = {{1,0,0}, {0,1,0}, {0,0,1}};
    std::vector <std::vector <double>> vdirections = {{0,1,1}, {1,0,1}, {1,1,0}};
    std::vector <std::vector <int>> vindex = {{1,2},{0,2},{0,1}};
    std::vector <double> signs = {1, -1};
    
    double spacing[3];
    this->Imp->GiveSpacing(spacing);
    
    for (auto v: EdgeVertices) {
        for (int i=0; i<3; i++) {
            // Find neighbour
            double c[3];
            for (int j=0; j<3; j++) c[j] = v->get_c(j) + testdirections.at(i).at(j)*spacing[j];
            VertexType *Neighbour = this->Mesh->VertexOctreeRoot->FindVertexByCoords(c[0], c[1], c[2]);
            
            if (std::find(EdgeVertices.begin(), EdgeVertices.end(), Neighbour) != EdgeVertices.end()) {
                LOG ("Found Neightbour %p for %p\n", Neighbour, v);
                double cm[3];
                for (int j=0; j<3; j++) cm[j]=v->get_c(j) + testdirections.at(i).at(j)*spacing[j]/2;
                
                // Check phases surrounding cm
                double delta[3];
                delta[0] = spacing[0]*vdirections.at(i)[0] / 2.0;
                delta[1] = spacing[1]*vdirections.at(i)[1] / 2.0;
                delta[2] = spacing[2]*vdirections.at(i)[2] / 2.0;
                std::vector <int> Phases;
                
                for (auto s1: signs) {
                    for (auto s2: signs) {
                        double testpoint[3];
                        testpoint[0] = cm[0];
                        testpoint[1] = cm[1];
                        testpoint[2] = cm[2];
                        
                        testpoint[vindex.at(i).at(0)] = testpoint[vindex.at(i).at(0)] + s1*delta[vindex.at(i)[0]];
                        testpoint[vindex.at(i).at(1)] = testpoint[vindex.at(i).at(1)] + s2*delta[vindex.at(i)[1]];
                        
                        int matid = this->Imp->GiveMaterialIDByCoordinate(testpoint[0], testpoint[1], testpoint[2]);
                        Phases.push_back(matid);
                        
                    }
                }
                
                // If 3 or 4 phases surrounds cm, this is a PhaseEdge
                std::sort(Phases.begin(), Phases.end());
                Phases.erase( std::unique(Phases.begin(), Phases.end()), Phases.end());
                
                if (Phases.size()>=3) {
                    LOG ("PhaseEdge found. Add it to the list and edd PhaseEdge to Vertex\n", 0);
                    AddPhaseEdge({v, Neighbour}, Phases);
                }
                
            }
        }
    }
    
    STATUS ("\tSort and fix non-connected edges...\n", 0);
    LOG("Fix and sort edges\n", 0);
    
    // Ensure that only have inner connected edges. I.e. max two vertices not connected to any other vertex on the edge
    unsigned int i=0;
    while (i<this->PhaseEdges.size()) {
        std::vector <PhaseEdge*> *FixedEdges = new std::vector <PhaseEdge*>;
        
        this->PhaseEdges.at(i)->SortAndFixBrokenEdge(FixedEdges);
        
        // Erase current PhaseEdge and replace it with the ones in FixedEdges
        // TODO: Should really delete this...
        //delete this->PhaseEdges.at(i);
        this->PhaseEdges.erase(this->PhaseEdges.begin() + i);
        this->PhaseEdges.insert(this->PhaseEdges.begin() + i, FixedEdges->begin(), FixedEdges->end());
        
        i=i+FixedEdges->size();
        
        // Cleanup
        /*for (auto fe: *FixedEdges) {
            delete fe;
        }
        delete FixedEdges;*/
        
    }
    
    STATUS("\tSplit phase edges at shared points\n",0);
    // Split phase edges at shared points
    
    // Add all vertices to a list from unique lists of vertices of each PhaseEdge
    std::vector <VertexType *> VertexList;
    for (auto e: this->PhaseEdges) {
        std::vector <VertexType *> FlatList = e->GetFlatListOfVertices();
        VertexList.insert(VertexList.end(), FlatList.begin(), FlatList.end());
    }
    
    // Count occurences of each vertex
    std::map <VertexType*, int> Counter;
    
    for (auto v: VertexList) {
        if (Counter.find(v) == Counter.end()) {
            Counter[v] = 1;
        } else {
            Counter[v] = Counter[v] + 1;
        }
    }
    
    // Add vertices occuring more than once to a list
    std::vector <VertexType *> SharedVertices;
    for (auto v: Counter) {
        if (v.second > 1) {
            SharedVertices.push_back(v.first);
        }
    }
    
    // Split edges at SharedVertices
    for (auto v: SharedVertices) {
        unsigned int i=0;
        
        //for (auto p: this->PhaseEdges) {
        while (i<this->PhaseEdges.size()) {
            PhaseEdge *p = PhaseEdges.at(i);
            
            std::vector<PhaseEdge*> SplitEdges;
            if (p->SplitAtVertex(v, &SplitEdges)) {
                
                this->PhaseEdges.insert(this->PhaseEdges.end(), SplitEdges.begin(), SplitEdges.end());
                this->PhaseEdges.erase(this->PhaseEdges.begin()+i);
                
                // Free memory from old PhaseEdge
                delete p;
                
            }
            i++;
        }
        
    }
    
    
    // Add PhaseEdges to surfaces
    // TODO: Performance can be increased by making sure that the vertices in Surfaces and PhaseEdges are sorted. Also, save the sorted list in PhaseEdges
    int x=0;
    for (Surface* s: this->Surfaces) {
        std::vector <VertexType *> SurfaceVertices;
        for (VertexType *v: s->Vertices) {
            SurfaceVertices.push_back(v);
        }
        std::sort (SurfaceVertices.begin(), SurfaceVertices.end());
        int SurfacePhases[2];
        SurfacePhases[0] = s->Phases[0];
        SurfacePhases[1] = s->Phases[1];
        std::sort(SurfacePhases, SurfacePhases+2);
        
        for (PhaseEdge *p: this->PhaseEdges) {
            
            std::vector<int> PhaseEdgePhases = p->Phases;
            std::sort(PhaseEdgePhases.begin(), PhaseEdgePhases.end());
            
            // If the surface phases are a subset of the edges phases, they might be connected
            if (std::includes(PhaseEdgePhases.begin(), PhaseEdgePhases.end(), SurfacePhases, SurfacePhases+2)) {
                std::vector <VertexType *> PhaseEdgeVertices = p->GetFlatListOfVertices();
                std::sort(PhaseEdgeVertices.begin(), PhaseEdgeVertices.end());
                
                if (std::includes(s->Vertices.begin(), s->Vertices.end(),
                                  PhaseEdgeVertices.begin(), PhaseEdgeVertices.end())) {
                    s->PhaseEdges.push_back(p);
                }
            }
            x++;
            
        }
    }
    
}

void Voxel2Tet :: SmoothEdgesIndividually(double c, double alpha, double charlength, bool Automatic_c)
{
    STATUS("Smooth edges (individually)\n", 0);
    for (unsigned int i=0; i<this->PhaseEdges.size(); i++) {
        LOG ("Smooth edge %i\n", i);
        PhaseEdge *e = this->PhaseEdges.at(i);
        e->Smooth(this->Mesh, c, alpha, charlength, Automatic_c);
    }
}

void Voxel2Tet :: SmoothEdgesSimultaneously(double c, double alpha, double charlength, bool Automatic_c)
{
    STATUS("Smooth edges (simultaneously)\n", 0);
    
    // We want to sort several lists according to one list. (http://stackoverflow.com/questions/1723066/c-stl-custom-sorting-one-vector-based-on-contents-of-another)
    struct VertexConnectivity {
        VertexType *v;
        std::vector<VertexType *> Connections;
        std::array<bool,3> FixedDirections;
    };
    
    struct by_vertexptr {
        bool operator()(VertexConnectivity *a, VertexConnectivity *b) { return (a->v < b->v); }
    };
    
    std::vector<VertexConnectivity *> VertexConnections;
    
    // Collect all vertices on phase edges
    int j=0;
    for (PhaseEdge *p: this->PhaseEdges) {
        std::vector<std::vector<VertexType *>> Connections;
        std::vector<std::array<bool,3>> FixedDirectionsList;
        
        std::vector<VertexType *> EdgeVertices = p->GetFlatListOfVertices();
        p->GiveTopologyLists(&Connections, &FixedDirectionsList);
        
        for (unsigned int i=0; i<EdgeVertices.size(); i++) {
            VertexType *v = EdgeVertices.at(i);
            VertexConnectivity *vc = new VertexConnectivity;
            vc->v = v;
            vc->Connections.insert(vc->Connections.begin(), Connections.at(i).begin(), Connections.at(i).end());
            vc->FixedDirections = FixedDirectionsList.at(i);
            VertexConnections.push_back(vc);
        }
        j++;
    }
    
    // Sort list by VertexType address and merge information of duplicate vertices
    std::sort(VertexConnections.begin(), VertexConnections.end(), by_vertexptr());
    
    unsigned int i=0;
    while (i<VertexConnections.size()) {
        // Compare element i to element i+1. If the elements points to the same vertex, merge connections and remove element i+1
        if (i<(VertexConnections.size()-1)) {
            if (VertexConnections.at(i)->v == VertexConnections.at(i+1)->v) {
                VertexConnectivity *thisvc = VertexConnections.at(i);
                VertexConnectivity *nextvc = VertexConnections.at(i+1);
                
                thisvc->Connections.insert(thisvc->Connections.end(), nextvc->Connections.begin(), nextvc->Connections.end());
                
                // Uniqueify connections
                std::sort(thisvc->Connections.begin(), thisvc->Connections.end());
                std::vector<VertexType *>::iterator it;
                it = std::unique(thisvc->Connections.begin(), thisvc->Connections.end());
                if (it!=thisvc->Connections.end()) {
                    thisvc->Connections.erase(it);
                }
                
                // Erase next connection
                VertexConnections.erase(VertexConnections.begin()+i+1);
            } else {
                i++;
            }
        } else i++;
        
    }
    
    // Build connectivity and FixedDirectionsLists
    std::vector<VertexType *> VertexList;
    std::vector<std::vector<VertexType *>> Connections;
    std::vector<bool> FixedDirectionsList;
    
    for (unsigned int i=0; i<VertexConnections.size(); i++) {
        VertexType *v = VertexConnections.at(i)->v;
        
        VertexList.push_back(v);
        Connections.push_back(VertexConnections.at(i)->Connections);
        
        // Determine which directions are locked TODO: This should be done elsewhere
        for (int j=0; j<3; j++) {
            if ( (v->get_c(j) > (this->Imp->GiveBoundingBox().maxvalues[j]-eps)) | (v->get_c(j) < (this->Imp->GiveBoundingBox().minvalues[j]+eps))) {
                v->Fixed[j]=true;
            } else {
                v->Fixed[j]=false;
            }
        }
        FixedDirectionsList.push_back(false);
    }
    
    for (VertexConnectivity *v: VertexConnections) {
        delete v;
    }
    
    SpringSmooth(VertexList, FixedDirectionsList, Connections, c, alpha, charlength, Automatic_c, this->Mesh);
    
}

void Voxel2Tet :: SmoothSurfaces(double c, double alpha, double charlength, bool Automatic_c)
{
    STATUS("Smooth surfaces\n", 0);

    int i=1;
    for (Surface *s: this->Surfaces) {
        STATUS("Smoothing surface %i (%i)\n", i, this->Surfaces.size());
        s->Smooth(this->Mesh, c, alpha, charlength, Automatic_c);
        i++;
    }
}

void Voxel2Tet :: SmoothAllAtOnce(double c, double alpha, double charlength, bool Automatic_c)
{
    STATUS("Smooth complete structure\n", 0);
    
    std::vector<std::vector<VertexType *>> Connections;
    std::vector<bool> FixedDirectionsList;
    
    // Create connections matrix
    for (unsigned int i=0; i<this->Mesh->Vertices.size(); i++) {
        VertexType *ThisVertex = this->Mesh->Vertices.at(i);
        
        // Find connected vertices
        std::vector <VertexType*> NeighbouringVertices = ThisVertex->FetchNeighbouringVertices();
        std::sort (NeighbouringVertices.begin(), NeighbouringVertices.end());
        std::vector <VertexType*> ConnectedVertices;
        
        // Create list of indices of connected vertices
        std::set_intersection(NeighbouringVertices.begin(), NeighbouringVertices.end(),
                              this->Mesh->Vertices.begin(), this->Mesh->Vertices.end(), back_inserter(ConnectedVertices));
        Connections.push_back(NeighbouringVertices);
        
        std::array<bool,3> FixedDirections;
        // Lock vertices on boundary surfaces such that they only move in the plane
        FixedDirections = {false, false, false};
        for (int j=0; j<3; j++) {
            if (ThisVertex->get_c(j)>=(this->Imp->GiveBoundingBox().maxvalues[j])-eps) {
                FixedDirections[j] = true;
            }
            if (ThisVertex->get_c(j)<=(this->Imp->GiveBoundingBox().minvalues[j])+eps) {
                FixedDirections[j] = true;
            }
        }
        
        FixedDirectionsList.push_back(false);
        
    }
    
    SpringSmooth(this->Mesh->Vertices, FixedDirectionsList, Connections, c, alpha, charlength, Automatic_c);
}

PhaseEdge* Voxel2Tet :: AddPhaseEdge(std::vector<VertexType*> EdgeSegment, std::vector<int> Phases)
{
    // Ensure that Phases argument is unique
    std::sort(Phases.begin(), Phases.end());
    PhaseEdge* ThisPhaseEdge=NULL;
    Phases.erase( std::unique(Phases.begin(), Phases.end()), Phases.end());
    
    // Find PhaseEdge
    for (auto pe: this->PhaseEdges) {
        std::vector <int> PhaseDiff;
        std::set_difference(Phases.begin(), Phases.end(),
                            pe->Phases.begin(), pe->Phases.end(),
                            back_inserter(PhaseDiff));
        if (PhaseDiff.size()==0) {
            ThisPhaseEdge = pe;
            break;
        }
    }
    
    // If PhaseEdge does not exists, create it
    if (ThisPhaseEdge==NULL) {
        ThisPhaseEdge = new PhaseEdge(this->Opt);
        ThisPhaseEdge->Phases = Phases;
        this->PhaseEdges.push_back(ThisPhaseEdge);
    }
    
    for (VertexType* v: EdgeSegment) {
        v->AddPhaseEdge(ThisPhaseEdge);
    }
    
    ThisPhaseEdge->AddPhaseEdgeSegment(EdgeSegment.at(0), EdgeSegment.at(1));
    //ThisPhaseEdge->EdgeSegments.push_back({EdgeSegment.at(0), EdgeSegment.at(1)});
    return ThisPhaseEdge;
    
}

void Voxel2Tet :: AddSurfaceSquare(std::vector<int> VertexIDs, std::vector<int> phases, int normalphase)
{
    // Check is surface exists
    Surface *ThisSurface = NULL;
    int SurfaceID;
    for (unsigned int i=0; i<this->Surfaces.size(); i++) {
        if ( ( (this->Surfaces.at(i)->Phases[0] == phases.at(0) ) & ( this->Surfaces.at(i)->Phases[1] == phases.at(1) ) ) |
             ( (this->Surfaces.at(i)->Phases[0] == phases.at(1) ) & ( this->Surfaces.at(i)->Phases[1] == phases.at(0) ) ) ) {
            ThisSurface = this->Surfaces.at(i);
            SurfaceID = i;
            break;
        }
    }
    
    // If not, create it and add it to the list
    if (ThisSurface==NULL) {
        ThisSurface = new Surface(phases.at(0), phases.at(1), this->Opt);
        this->Surfaces.push_back(ThisSurface);
        SurfaceID = this->Surfaces.size()-1;
    }
    
    // Create square (i.e. two triangles)
    TriangleType *triangle0, *triangle1;
    triangle0 = Mesh->AddTriangle({VertexIDs.at(0), VertexIDs.at(1), VertexIDs.at(2)});
    triangle1 = Mesh->AddTriangle({VertexIDs.at(2), VertexIDs.at(1), VertexIDs.at(3)});

    triangle0->InterfaceID = triangle1->InterfaceID = SurfaceID;

    // Check phases - This is kind of a dirty hack. Should not be neccessary to check the phase on the positive side once again.

    std::array<double, 3> cm = triangle1->GiveCenterOfMass();
    std::array<double, 3> normal = triangle0->GiveNormal();
    double spacing[3];
    this->Imp->GiveSpacing(spacing);

    for (int i=0; i<3; i++) {
        normal[i]=normal[i]*spacing[i]*.5;
    }

    int pPhase;
    pPhase = this->Imp->GiveMaterialIDByCoordinate(cm[0]+normal[0], cm[1]+normal[1], cm[2]+normal[2]);

    triangle0->PosNormalMatID = triangle1->PosNormalMatID = pPhase;

    if (phases[0] == pPhase) {
        triangle0->NegNormalMatID = triangle1->NegNormalMatID = phases[1];
    } else {
        triangle0->NegNormalMatID = triangle1->NegNormalMatID = phases[0];
    }

    // Update surface
    ThisSurface->AddTriangle(triangle0);
    ThisSurface->AddTriangle(triangle1);
    for (int i=0; i<4; i++) {
        ThisSurface->AddVertex(this->Mesh->Vertices.at(VertexIDs.at(i)));
    }
    
}

double Voxel2Tet::GetListOfVolumes(std::vector<double> &VolumeList, std::vector<int> &PhaseList)
{

    VolumeList.clear();
    PhaseList.clear();

    double TotalVolume = 0.0;
    for (size_t i=0; i<Volumes.size(); i++) {
        double v = this->Volumes.at(i)->ComputeVolume();
        VolumeList.push_back(v);
        PhaseList.push_back(this->Volumes.at(i)->Phase);
        TotalVolume =+ v;
    }
    return TotalVolume;

}

Volume * Voxel2Tet::FindVolumeContainingPoint (std::array<double, 3> P)
{
    for (Volume *v: this->Volumes) {
        if (v->IsPointInside(P)) {
            return v;
        }
    }
    return NULL;
}

void Voxel2Tet::Process()
{
    STATUS ("Proccess content\n", 0);
    
#ifdef OPENMP
    STATUS ("\tUsing OpenMP and %u threads\n", omp_get_max_threads());
#endif
    TetGenCaller Generator;
    Generator.Mesh = this->Mesh;

    int outputindex = 0;
    
    this->FindSurfaces();
    
    std::vector<std::vector<double>> PhaseVolumes;
    std::vector<double> CurrentVolumes;
    std::vector<int> PhaseList;

    double TotalVolume = GetListOfVolumes(CurrentVolumes, PhaseList);
    LOG("Total volume: %f\n", TotalVolume);
    PhaseVolumes.push_back(CurrentVolumes);
    
    this->Mesh->ExportSurface(strfmt("/tmp/Voxeltest%u.vtp", outputindex++) , FT_VTK );
        
    this->FindEdges();

    clock_t t1 = clock();

    this->SmoothEdgesSimultaneously(edgespring_c, edgespring_alpha, 0, false);
    Generator.TestMesh();

    GetListOfVolumes(CurrentVolumes, PhaseList);
    PhaseVolumes.push_back(CurrentVolumes);

    this->Mesh->ExportSurface(strfmt("/tmp/Voxeltest%u.vtp", outputindex++), FT_VTK);
    this->SmoothSurfaces(spring_c, spring_alpha, 0, false);
    clock_t t2 = clock();
    Generator.TestMesh();

    GetListOfVolumes(CurrentVolumes, PhaseList);
    PhaseVolumes.push_back(CurrentVolumes);

    STATUS("Smoothing completed in %fs\n", float(t2-t1)/(double)CLOCKS_PER_SEC);
    
    this->Mesh->ExportSurface(strfmt("/tmp/Voxeltest%u.vtp", outputindex++), FT_VTK);
    this->Mesh->FlipAll();

    this->UpdateSurfaces();

    GetListOfVolumes(CurrentVolumes, PhaseList);
    PhaseVolumes.push_back(CurrentVolumes);

    this->Mesh->ExportSurface(strfmt("/tmp/Voxeltest%u.vtp", outputindex++), FT_VTK);
    this->Mesh->CoarsenMeshImproved();

    this->UpdateSurfaces();

    GetListOfVolumes(CurrentVolumes, PhaseList);
    PhaseVolumes.push_back(CurrentVolumes);

    this->Mesh->ExportSurface(strfmt("/tmp/Voxeltest%u.vtp", outputindex++), FT_VTK);
    this->Mesh->FlipAll();
    this->UpdateSurfaces();

    GetListOfVolumes(CurrentVolumes, PhaseList);
    PhaseVolumes.push_back(CurrentVolumes);

    this->Mesh->ExportSurface(strfmt("/tmp/Voxeltest%u.vtp", outputindex++), FT_VTK);

    for (int p: PhaseList) {
        printf("%u\t\t", p);
    }
    
    printf("\n");

    for (size_t i=0; i<PhaseVolumes.size(); i++) {
        for (size_t j=0; j<PhaseVolumes.at(i).size(); j++) {
            printf("%f\t", PhaseVolumes.at(i).at(j));
        }
        printf("\n");
    }
#if EXPORT_MESH_COARSENING == 1
    dooutputlogmesh(*this->Mesh, (char*) "/tmp/finalcoarsening%u.vtp", 0);
#endif
    
}

}
