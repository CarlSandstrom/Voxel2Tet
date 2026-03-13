#include "TetPipeline.h"

namespace voxel2tet
{
TetPipeline::TetPipeline(TetrahedralizeDelegate TetrahedralizeStep)
{
    this->TetrahedralizeStep = TetrahedralizeStep;
}

void TetPipeline::Tetrahedralize() const
{
    this->TetrahedralizeStep();
}
}
