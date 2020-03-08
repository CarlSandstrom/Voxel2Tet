#ifndef SIMPLEEXPORTER_H
#define SIMPLEEXPORTER_H

#include "Exporter.h"

namespace voxel2tet
{

/**
 * @brief The SimpleExporter class provides functionality to export to a simple ASCII file.
 *
 * The basic format consists of three blocks. The first block holds the vertices, the second the edges and the third the triangles.
 * Each block starts vith an integer telling the number of comming lines for that block.
 *
 */
class SimpleExporter : public Exporter
{
public:

    /**
     * @copydoc Exporter::Exporter
     */
    SimpleExporter(std::vector<TriangleType *> *Triangles, std::vector<VertexType *> *Vertices,
                   std::vector<EdgeType *> *Edges, std::vector<TetType *> *Tets);

    virtual void WriteSurfaceData(std::string Filename, int VolumeID);

    virtual void WriteVolumeData(std::string Filename)
    {}

};

}
#endif // SIMPLEEXPORTER_H
