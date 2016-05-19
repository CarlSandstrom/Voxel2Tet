#include <vector>
#include <math.h>
#include <sstream>

#include "MeshComponents.h"
#include "MeshData.h"
#include "TetGenCaller.h"

namespace voxel2tet
{

void SpringSmooth(std::vector<VertexType*> Vertices, std::vector<bool> FixedDirections,
                  std::vector<std::vector<VertexType*>> Connections, double K, MeshData *Mesh=NULL);


}
