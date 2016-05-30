#ifndef VOLUME_H
#define VOLUME_H

#include <vector>

#include "Surface.h"
#include "MeshComponents.h"
#include "MiscFunctions.h"
#include "TriTriIntersect.h"

namespace voxel2tet
{

class Volume
{
public:
    Volume();

    Volume(int Phase);

    std::vector<Surface*> Surfaces;

    std::vector<TriangleType*> GiveTriangles();

    std::vector<VertexType*> GiveVertices();

    /**
     * @brief IsPointInside determines if point P lies withing the volume
     *
     * Uses the algorithm described in http://www.academia.edu/841865/An_efficient_point_classification_algorithm_for_triangle_meshes
     *
     * @param P [in] Point coodrinates
     * @return True if point is within the volume, false if outside.
     */
    bool IsPointInside(std::array<double, 3> P);

    double ComputeVolume();

    int Phase;
};

}
#endif // VOLUME_H
