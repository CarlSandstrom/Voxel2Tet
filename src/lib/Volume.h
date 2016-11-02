#ifndef VOLUME_H
#define VOLUME_H

#include <vector>

#include "Surface.h"
#include "MeshComponents.h"
#include "MiscFunctions.h"
#include "TriTriIntersect.h"

namespace voxel2tet
{

/**
 * @brief The Volume class holds information for a volume produces using Surface objects. A Volume contains only one material.
 */
class Volume
{
public:
    /**
     * @brief Constructor
     */
    Volume();

    /**
     * @brief Constructor
     * @param Phase The material ID of the material contained within the volume
     */
    Volume(int Phase);

    /**
     * @brief List of Surface objects that describe the boundary of the volume.
     */
    std::vector<Surface *> Surfaces;

    /**
     * @brief GiveTriangles returns all triangles on the surface of the volume
     * @return Vector of triangles
     */
    std::vector<TriangleType *> GiveTriangles();

    /**
     * @brief GiveVertices returns all vertices on the surface of the volume
     * @return Vector of vertices
     */
    std::vector<VertexType *> GiveVertices();

    /**
     * @brief IsPointInside determines if point P lies withing the volume
     *
     * Uses the algorithm described in http://www.academia.edu/841865/An_efficient_point_classification_algorithm_for_triangle_meshes
     *
     * @param P [in] Point coodrinates
     * @return True if point is within the volume, false if outside.
     */
    bool IsPointInside(std::array<double, 3> P);

    /**
     * @brief Computes the volume of the Volume object.
     *
     * The simple algorithm uses the signed volume of the tetrahedtral formed by adding origo to the set of vertices
     * for each triangle. The sign of the contribution is the decided by the orientation of the tetrahedron.
     *
     * @return Volume
     */
    double ComputeVolume();

    /**
     * @brief Material ID of material in volume
     */
    int Phase;
};
}
#endif // VOLUME_H
