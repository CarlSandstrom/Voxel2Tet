#ifndef STLEXPORTER_H
#define STLEXPORTER_H

#include <iostream>
#include <iomanip>

#include "Exporter.h"
#include "MeshComponents.h"
#include "MiscFunctions.h"
#include "Exporter.h"

namespace voxel2tet
{
/**
 * @brief The STLExporter class provides functionality to export to the .STL file format
 */
class STLExporter : public Exporter
{
public:
    /**
     * @copydoc Exporter::Exporter
     */
    STLExporter(std::vector<TriangleType *> *Triangles, std::vector<VertexType *> *Vertices,
                std::vector<EdgeType *> *Edges, std::vector<TetType *> *Tets);

    virtual void WriteSurfaceData(std::string Filename, int VolumeID);

    virtual void WriteVolumeData(std::string Filename)
    {}
};
}

#endif // STLEXPORTER_H
