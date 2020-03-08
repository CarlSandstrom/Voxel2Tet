#ifndef TETGETEXPORTER_H
#define TETGETEXPORTER_H

#include <string>
#include <array>

#include "Exporter.h"
#include "MeshComponents.h"
#include "MiscFunctions.h"
#include "Exporter.h"

namespace voxel2tet
{

/**
 * @brief The TetGenExporter class provides functionality to export to TetGens .poly format
 */
class TetGenExporter : public Exporter
{
public:
    /**
     * @copydoc Exporter::Exporter
     */
    TetGenExporter(std::vector<TriangleType *> *Triangles, std::vector<VertexType *> *Vertices,
                   std::vector<EdgeType *> *Edges, std::vector<TetType *> *Tets);

    virtual void WriteSurfaceData(std::string Filename, int VolumeID=-1);

    virtual void WriteVolumeData(std::string Filename)
    {}
};
}
#endif // TETGETEXPORTER_H
