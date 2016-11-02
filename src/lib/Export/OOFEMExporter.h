#ifndef OOFEMEXPORTER_H
#define OOFEMEXPORTER_H

#include "Exporter.h"

namespace voxel2tet
{

/**
 * @brief The OOFEMExporter class provides functionality for exporting to OOFEMs .in format.
 *
 * It exports the volume and supplies several sets for the user.
 *
 *  - Set 1: Set of all tetrahedrons
 *  - Set 2: Set of all nodes on the outer boundary
 *  - Sets 3-8: Sets of nodes on each outer boundary. E.g. one set for all nodes on the X-most surface etc.
 *  - Sets 9-12: Sets of element facets on each outer boundary.
 *
 */
class OOFEMExporter : public Exporter
{
public:
    /**
     * @copydoc Exporter::Exporter
     */
    OOFEMExporter(std::vector<TriangleType *> *Triangles, std::vector<VertexType *> *Vertices,
                  std::vector<EdgeType *> *Edges, std::vector<TetType *> *Tets);

    virtual void WriteSurfaceData(std::string Filename)
    {}

    virtual void WriteVolumeData(std::string Filename);
};
}

#endif // OOFEMEXPORTER_H
