#include "PhaseEdge.h"
#include <array>

namespace voxel2tet
{

PhaseEdge::PhaseEdge()
{

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
            PhaseEdge *NewPhaseEdge = new PhaseEdge;
            NewPhaseEdge->Phases = this->Phases;

            std::copy(this->EdgeSegments.begin()+EdgeStart, this->EdgeSegments.begin()+i, std::back_inserter( NewPhaseEdge->EdgeSegments ));

            SplitEdges->push_back(NewPhaseEdge);
            EdgeStart=i;
        }
    }

    // Copy remaining part.
    PhaseEdge *NewPhaseEdge = new PhaseEdge;
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

        PhaseEdge *NewPhaseEdge=new PhaseEdge();
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

void PhaseEdge :: Smooth()
{
    if (this->EdgeSegments.size()==1) return;

    // Push end vertices to FixedVertices.
    this->FixedVertices.push_back( this->EdgeSegments.at(0).at(0) );
    this->FixedVertices.push_back( this->EdgeSegments.at( this->EdgeSegments.size()-1 ).at(1) );

    std::vector<std::array<double, 3>> CurrentPositions;
    std::vector<std::array<double, 3>> PreviousPositions;
    std::vector<VertexType *> FlatList = this->GetFlatListOfVertices();

    for (unsigned int i=0; i<FlatList.size(); i++) {
        LOG("C = (%f, %f, %f)\n", FlatList.at(i)->c[0], FlatList.at(i)->c[1], FlatList.at(i)->c[2]);
        std::array<double, 3> cp;// = {FlatList.at(i)->c[0], FlatList.at(i)->c[1], FlatList.at(i)->c[2]};
        std::array<double, 3> pp;
        for (int j=0; j<3; j++) {
            cp.at(j) = pp.at(j) = FlatList.at(i)->c[j];
        }
        CurrentPositions.push_back(cp);
        PreviousPositions.push_back(pp);
    }

}

}
