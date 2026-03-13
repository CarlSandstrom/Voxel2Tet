#ifndef TETPIPELINE_H
#define TETPIPELINE_H

#include <functional>

namespace voxel2tet
{
class TetPipeline
{
public:
    typedef std::function<void(void)> TetrahedralizeDelegate;

    explicit TetPipeline(TetrahedralizeDelegate TetrahedralizeStep);

    void Tetrahedralize() const;

private:
    TetrahedralizeDelegate TetrahedralizeStep;
};
}

#endif // TETPIPELINE_H
