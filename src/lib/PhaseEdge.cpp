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

void PhaseEdge :: Smooth()
{
    if (this->EdgeSegments.size()==1) return;
    double K = this->Opt->GiveDoubleValue("spring_const");

    // Push end vertices to FixedVertices.
    if (this->EdgeSegments.at(0).at(0) != this->EdgeSegments.at( this->EdgeSegments.size()-1 ).at(1)) {
        this->FixedVertices.push_back( this->EdgeSegments.at(0).at(0) );
        this->FixedVertices.push_back( this->EdgeSegments.at( this->EdgeSegments.size()-1 ).at(1) );
    }

    std::vector<std::array<double, 3>> CurrentPositions;
    std::vector<std::array<double, 3>> PreviousPositions;
    std::vector<VertexType *> FlatList = this->GetFlatListOfVertices();

    std::vector<int> FixedVerticesIndices = FindSubsetIndices(FlatList, this->FixedVertices);

    for (unsigned int i=0; i<FlatList.size(); i++) {
        std::array<double, 3> cp;
        std::array<double, 3> pp;
        for (int j=0; j<3; j++) {
            cp.at(j) = pp.at(j) = FlatList.at(i)->c[j];
        }
        CurrentPositions.push_back(cp);
        PreviousPositions.push_back(pp);
    }

    int itercount = 0;
    while (itercount < 100) {

        for (unsigned int i=0; i<FlatList.size(); i++) {
            // Move only if the vertex is not in the FixedVerticesIndecis list
            int previndex = i-1;
            int nextindex = i+1;
            if ((previndex == -1)) previndex = FlatList.size()-1;
            if ((nextindex == FlatList.size())) nextindex = 0;
            if (std::find(FixedVerticesIndices.begin(), FixedVerticesIndices.end(), i)==FixedVerticesIndices.end()) {
                // Laplace
                for (int j=0; j<3; j++) {
                    CurrentPositions.at(i)[j] = (PreviousPositions.at(previndex)[j] + PreviousPositions.at(nextindex)[j]) / 2.0;
                }

                // Pull back

                std::array<double, 3> delta, unitdelta;
                for (int j=0; j<3; j++) delta[j]=CurrentPositions.at(i)[j]-FlatList.at(i)->c[j];
                double d0 = sqrt( pow(delta[0], 2) + pow(delta[1], 2) + pow(delta[2], 2) );

                if (d0>1e-8) {
                    for (int j=0; j<3; j++) unitdelta[j] = delta[j]/d0;

                    double F = d0*1.0;
                    double d = d0;

                    double change=1e8;
                    while (change>1e-8) {
                        double NewDelta = F*1.0/exp(pow(d , 2) / K);
                        change = fabs(d-NewDelta);
                        d = NewDelta;
                    }

                    for (int j=0; j<3; j++) CurrentPositions.at(i)[j] = FlatList.at(i)->c[j] + unitdelta[j]*d;
                }
            }
        }




        // Update previous positions
        for (unsigned int j=0; j<PreviousPositions.size(); j++) {
            PreviousPositions.at(j)=CurrentPositions.at(j);
        }
        itercount ++;

    }

    //Update vertices
    for (unsigned int i=0; i<FlatList.size(); i++) {
        for (int j=0; j<3; j++) {
            FlatList.at(i)->c[j] = CurrentPositions.at(i)[j];
        }
    }

}

}
