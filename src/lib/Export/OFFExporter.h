#ifndef OFFEXPORTER_H
#define OFFEXPORTER_H

#include <string>
#include <array>

#include "Exporter.h"
#include "MeshComponents.h"
#include "MiscFunctions.h"
#include "Exporter.h"

namespace voxel2tet
{

/**
 * @brief The OFFExporter class provides functionality to export to the .OFF file format used by e.g. TetGen.
 */
class OFFExporter : public Exporter
{
public:
    /**
     * @copydoc Exporter::Exporter
     */
    OFFExporter(std :: vector< TriangleType * > *Triangles, std :: vector< VertexType * > *Vertices, std :: vector< EdgeType * > *Edges, std :: vector< TetType * > *Tets);
    virtual void WriteSurfaceData(std :: string Filename);
    virtual void WriteVolumeData(std :: string Filename) {}
};
}

#endif // OFFEXPORTER_H
