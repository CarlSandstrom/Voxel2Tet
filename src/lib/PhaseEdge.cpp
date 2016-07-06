#include "PhaseEdge.h"

namespace voxel2tet
{
PhaseEdge :: PhaseEdge(Options *Opt, SpringSmoother *EdgeSmoother)
{
    this->Opt = Opt;
    this->EdgeSmoother = EdgeSmoother;
}

std :: vector< VertexType * >PhaseEdge :: GetFlatListOfVertices()
{
    // TODO: If we assume that the lists are ordered, this can be optimized
    std :: vector< VertexType * >FlatList;

    for ( auto e : this->EdgeSegments ) {
        FlatList.push_back( e.at(0) );
        FlatList.push_back( e.at(1) );
    }

    // Uniquify
    FlatList.erase( std :: unique( FlatList.begin(), FlatList.end() ), FlatList.end() );
    return FlatList;
}

std :: vector< VertexType * >PhaseEdge :: GiveVerticesConnectedToVertex(VertexType *v)
{
    std :: vector< VertexType * >ResultsList;
    for ( std :: array< VertexType *, 2 >e : this->EdgeSegments ) {
        if ( e [ 0 ] == v ) {
            ResultsList.push_back(e [ 1 ]);
        }
        if ( e [ 1 ] == v ) {
            ResultsList.push_back(e [ 0 ]);
        }
    }
    return ResultsList;
}

void PhaseEdge :: SortAndFixBrokenEdge(std :: vector< PhaseEdge * > *FixedEdges)
{
    // For each EdgeSegment find the two connecting segments until all segments has been processed. If links
    // to the end vertices of a chain of segments cannot be found while there are segments in the list, a new
    // list of segments is created.

    FixedEdges->clear();

    while ( this->EdgeSegments.size() > 0 ) {
        std :: array< VertexType *, 2 >ThisLink = this->EdgeSegments.at(0);
        this->EdgeSegments.erase(this->EdgeSegments.begin(), this->EdgeSegments.begin() + 1);

        LOG("Find connections for vertices (%u, %u)\n", ThisLink.at(0)->ID, ThisLink.at(1)->ID);

        PhaseEdge *NewPhaseEdge = new PhaseEdge(this->Opt, this->EdgeSmoother);
        NewPhaseEdge->Phases = this->Phases;
        NewPhaseEdge->EdgeSegments.push_back(ThisLink);
        FixedEdges->push_back(NewPhaseEdge);

        // Find next link, i.e. the link containing the i:th vertex of the current segment
        for ( int i = 0; i < 2; i++ ) {
            VertexType *VertexToFind = ThisLink.at(i);

            bool LastConnectionFound = false;
            while ( !LastConnectionFound ) {
                // Remove self from list of PhaseEdges
                LOG("Remove PhaseEdge object from list for vertex %u\n", VertexToFind->ID);
                std :: vector< PhaseEdge * > :: iterator ErasePhaseEdges = std :: remove(VertexToFind->PhaseEdges.begin(), VertexToFind->PhaseEdges.end(), this);
                VertexToFind->PhaseEdges.erase( ErasePhaseEdges, VertexToFind->PhaseEdges.end() );
                VertexToFind->PhaseEdges.push_back(NewPhaseEdge);

                std :: array< VertexType *, 2 >NextLink;
                bool NextLinkFound = false;
                for ( unsigned int j = 0; j < this->EdgeSegments.size(); j++ ) {
                    if ( ( this->EdgeSegments.at(j).at(0) == VertexToFind ) | ( this->EdgeSegments.at(j).at(1) == VertexToFind ) ) {
                        LOG( "Next link found (%p, %p)\n", this->EdgeSegments.at(j).at(0), this->EdgeSegments.at(j).at(1) );

                        NextLink = this->EdgeSegments.at(j);
                        NextLinkFound = true;
                        this->EdgeSegments.erase(this->EdgeSegments.begin() + j, this->EdgeSegments.begin() + j + 1);
                        break;
                    }
                }

                if ( NextLinkFound ) {
                    std :: array< VertexType *, 2 > *NewEdgeSegment;
                    VertexType *NextLastVertex;

                    if ( NextLink.at(0) == VertexToFind ) {
                        NextLastVertex = NextLink.at(1);
                    } else {
                        NextLastVertex = NextLink.at(0);
                    }

                    if ( i == 0 ) {
                        NewEdgeSegment = new std :: array< VertexType *, 2 > { {
                                                                                   NextLastVertex, VertexToFind
                                                                               } };
                        NewPhaseEdge->EdgeSegments.insert(NewPhaseEdge->EdgeSegments.begin(), * NewEdgeSegment);
                    } else {
                        NewEdgeSegment = new std :: array< VertexType *, 2 > { {
                                                                                   VertexToFind, NextLastVertex
                                                                               } };
                        NewPhaseEdge->EdgeSegments.push_back(* NewEdgeSegment);
                    }
                    VertexToFind = NextLastVertex;
                } else {
                    LOG("No more connections found for i=%u\n", i);
                    LastConnectionFound = true;
                }
            }

            // If no EdgeSegments are left, exit
            if ( this->EdgeSegments.size() == 0 ) {
                break;
            }
        }
        if ( this->EdgeSegments.size() > 0 ) {
            LOG("Items left, split\n", 0);
        } else {
            LOG("No items left. \n", 0);
        }
    }
}

bool PhaseEdge :: IsClosed()
{
    return this->EdgeSegments.at(0).at(0) == this->EdgeSegments.at(this->EdgeSegments.size() - 1).at(1);
}

void PhaseEdge :: GiveTopologyLists(std :: vector< std :: vector< VertexType * > > *Connections, std :: vector< std :: array< bool, 3 > > *FixedDirectionsList) // TODO: Should this really supply the FixedDirectionList? Is it used?
{
    Connections->clear();
    FixedDirectionsList->clear();

    // Push end vertices to FixedVertices.
    if ( this->EdgeSegments.at(0).at(0) != this->EdgeSegments.at(this->EdgeSegments.size() - 1).at(1) ) {
        this->FixedVertices.push_back( this->EdgeSegments.at(0).at(0) );
        this->FixedVertices.push_back( this->EdgeSegments.at(this->EdgeSegments.size() - 1).at(1) );
    }

    std :: vector< VertexType * >FlatList = this->GetFlatListOfVertices();
    bool closed = this->IsClosed();

    // Build connection matrix
    for ( unsigned int i = 0; i < FlatList.size(); i++ ) {
        std :: vector< VertexType * >MyConnections;

        signed int previndex = i - 1;
        unsigned int nextindex = i + 1;

        if ( closed ) {
            if ( i == 0 ) {
                previndex = FlatList.size() - 1;
                nextindex = i + 1;
            } else if ( i == FlatList.size() - 1 ) {
                previndex = i - 1;
                nextindex = 0;
            }
        }

        if ( previndex != -1 ) {
            MyConnections.push_back( FlatList.at(previndex) );
        }
        if ( nextindex != FlatList.size() ) {
            MyConnections.push_back( FlatList.at(nextindex) );
        }

        // Ignore edges not aligned with the principal axes
        Connections->push_back(MyConnections);

        std :: array< bool, 3 >FixedDirections;

        // Set start and end points to fixed
        if ( std :: find( this->FixedVertices.begin(), this->FixedVertices.end(), FlatList.at(i) ) == this->FixedVertices.end() ) {
            FixedDirections = { { false, false, false } };
        } else {
            FixedDirections = { { true, true, true } };
        }

        FixedDirectionsList->push_back(FixedDirections);
    }
}

void PhaseEdge :: Smooth(MeshData *Mesh)
{
    if ( this->EdgeSegments.size() == 1 ) {
        return;
    }

    std :: vector< VertexType * >FlatList = this->GetFlatListOfVertices();

    std :: vector< std :: vector< VertexType * > >Connections;
    std :: vector< std :: array< bool, 3 > >FixedDirectionsList;
    std :: vector< bool >FixedList;

    this->GiveTopologyLists(& Connections, & FixedDirectionsList);

    for ( size_t i = 0; i < FlatList.size(); i++ ) {
        FixedList.push_back(false);
    }

    this->EdgeSmoother->Smooth(FlatList, FixedList, Connections, Mesh);
}

void PhaseEdge :: AddPhaseEdgeSegment(VertexType *v1, VertexType *v2)
{
    EdgeSegments.push_back({ { v1, v2 } });
}
}
