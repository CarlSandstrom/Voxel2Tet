#include "PhaseEdge.h"
#include <array>

namespace voxel2tet
{

PhaseEdge::PhaseEdge(Options *Opt)
{
    this->Opt = Opt;
}

std::vector<VertexType*> PhaseEdge :: GetFlatListOfVertices()
{
    // TODO: If we assume that the lists are ordered, this can be optimized
    std::vector<VertexType*> FlatList;

    for (auto e: this->EdgeSegments) {
        FlatList.push_back(e.at(0));
        FlatList.push_back(e.at(1));
    }

    // Uniquify
    FlatList.erase( std::unique(FlatList.begin(), FlatList.end()), FlatList.end());
    return FlatList;
}

void PhaseEdge :: SplitAtVertex(VertexType *Vertex, std::vector<PhaseEdge*> *SplitEdges)
{
    int EdgeStart = 0;
    unsigned int i=0;

    for (; i<this->EdgeSegments.size(); i++) {

        std::vector<VertexType*> e;
        e = EdgeSegments.at(i);

        if ( (e.at(1)==Vertex) & (e.at(0)==Vertex)) {
            PhaseEdge *NewPhaseEdge = new PhaseEdge(this->Opt);
            NewPhaseEdge->Phases = this->Phases;

            std::copy(this->EdgeSegments.begin()+EdgeStart, this->EdgeSegments.begin()+i, std::back_inserter( NewPhaseEdge->EdgeSegments ));

            SplitEdges->push_back(NewPhaseEdge);
            EdgeStart=i;
        }
    }

    // Copy remaining part.
    PhaseEdge *NewPhaseEdge = new PhaseEdge(this->Opt);
    NewPhaseEdge->Phases = this->Phases;

    std::copy(this->EdgeSegments.begin()+EdgeStart, this->EdgeSegments.begin()+i, std::back_inserter( NewPhaseEdge->EdgeSegments ));
    SplitEdges->push_back(NewPhaseEdge);

}

void PhaseEdge :: SortAndFixBrokenEdge(std::vector<PhaseEdge*> *FixedEdges)
{

    while (this->EdgeSegments.size()>0) {
        std::vector<VertexType*> ThisLink = this->EdgeSegments.at(0);
        this->EdgeSegments.erase(this->EdgeSegments.begin(), this->EdgeSegments.begin()+1);

        LOG("Find connections for link (%p, %p)\n", ThisLink.at(0), ThisLink.at(1));

        PhaseEdge *NewPhaseEdge=new PhaseEdge(this->Opt);
        NewPhaseEdge->Phases = this->Phases;
        NewPhaseEdge->EdgeSegments.push_back(ThisLink);
        FixedEdges->push_back(NewPhaseEdge);

        for (int i=0; i<2; i++) {
            VertexType* VertexToFind = ThisLink.at(i);

            bool LastConnectionFound = false;
            while (!LastConnectionFound) {

                std::vector<VertexType*> NextLink;
                for (unsigned int j=0; j<this->EdgeSegments.size(); j++) {
                    if ((this->EdgeSegments.at(j).at(0)==VertexToFind) | (this->EdgeSegments.at(j).at(1)==VertexToFind)) {
                        LOG("Next link found (%p, %p)\n", this->EdgeSegments.at(j).at(0), this->EdgeSegments.at(j).at(1));

                        NextLink = this->EdgeSegments.at(j);
                        this->EdgeSegments.erase(this->EdgeSegments.begin()+j, this->EdgeSegments.begin()+j+1);
                        break;
                    }
                }

                if (NextLink.size()>0) {
                    std::vector<VertexType*> *NewEdgeSegment;
                    VertexType *NextLastVertex;

                    if (NextLink.at(0)==VertexToFind) {
                        NextLastVertex = NextLink.at(1);
                    } else {
                        NextLastVertex = NextLink.at(0);
                    }

                    if (i==0) {
                        NewEdgeSegment = new std::vector<VertexType*> {NextLastVertex, VertexToFind};
                        NewPhaseEdge->EdgeSegments.insert(NewPhaseEdge->EdgeSegments.begin(), *NewEdgeSegment);
                    } else {
                        NewEdgeSegment = new std::vector<VertexType*> {VertexToFind, NextLastVertex};
                        NewPhaseEdge->EdgeSegments.push_back(*NewEdgeSegment);
                    }
                    VertexToFind = NextLastVertex;

                } else {
                    LOG("No more connections found for i=%u\n", i);
                    LastConnectionFound = true;
                }
            }

        }
        if (this->EdgeSegments.size()>0) {
            LOG("Items left, split\n",0);
        } else {
            LOG("No items left. \n",0);
        }
    }
}

bool PhaseEdge :: IsClosed()
{
    return this->EdgeSegments.at(0).at(0) == this->EdgeSegments.at(this->EdgeSegments.size()-1).at(1);
}

void PhaseEdge :: GiveTopologyLists(std::vector<std::vector<VertexType *>> *Connections, std::vector<std::array<bool,3>> *FixedDirectionsList)
{
    Connections->clear();
    FixedDirectionsList->clear();

    // Push end vertices to FixedVertices.
    if (this->EdgeSegments.at(0).at(0) != this->EdgeSegments.at( this->EdgeSegments.size()-1 ).at(1)) {
        this->FixedVertices.push_back( this->EdgeSegments.at(0).at(0) );
        this->FixedVertices.push_back( this->EdgeSegments.at( this->EdgeSegments.size()-1 ).at(1) );
    }

    std::vector<VertexType *> FlatList = this->GetFlatListOfVertices();
    bool closed = this->IsClosed();

    // Build connection matrix
    for (unsigned int i=0; i<FlatList.size(); i++) {

        std::vector<VertexType *> MyConnections;

        unsigned int previndex = i-1;
        unsigned int nextindex = i+1;

        if (closed) {
            if (i==0) {
                previndex = FlatList.size()-1;
                nextindex = i+1;
            } else if (i==FlatList.size()-1) {
                previndex = i-1;
                nextindex = 0;
            }
        }

        if (previndex!=-1) MyConnections.push_back(FlatList.at(previndex));
        if (nextindex!=FlatList.size()) MyConnections.push_back(FlatList.at(nextindex));

        // Ignore edges not aligned with the principal axes


        Connections->push_back(MyConnections);

        std::array<bool,3> FixedDirections;


        if (std::find(this->FixedVertices.begin(), this->FixedVertices.end(), FlatList.at(i))==this->FixedVertices.end()  ) {//if ((i!=0) & (i!=(FlatList.size()-1))) {
            FixedDirections = {false, false, false};
        } else {
            FixedDirections = {true, true, true};
        }

        FixedDirectionsList->push_back(FixedDirections);
    }

}

void PhaseEdge :: Smooth(MeshData *Mesh)
{
    if (this->EdgeSegments.size()==1) return;
    double K = this->Opt->GiveDoubleValue("spring_const");

    std::vector<VertexType *> FlatList = this->GetFlatListOfVertices();

    std::vector<std::vector<VertexType *>> Connections;
    std::vector<std::array<bool,3>> FixedDirectionsList;

    this->GiveTopologyLists(&Connections, &FixedDirectionsList);

    SpringSmooth(FlatList, FixedDirectionsList, Connections, K, Mesh);

}

}
