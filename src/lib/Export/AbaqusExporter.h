#ifndef ABAQUSEXPORTER_H
#define ABAQUSEXPORTER_H

#include "Exporter.h"

namespace voxel2tet
{

/**
 * Exporter class for Abaqus.
 */
class AbaqusExporter : public Exporter
{
private:
    bool UsePhon;
public:
    /**
     * Constructor of Abaqus export class
     * @param Triangles List of triangles to export
     * @param Vertices List of vertices to export
     * @param Edges List of edges to export
     * @param Tets List of tetrahedrons to export
     * @param UsePhon If true, special sets for Phon are exported.
     * @return
     */
    AbaqusExporter(std::vector<TriangleType *> *Triangles, std::vector<VertexType *> *Vertices,
                   std::vector<EdgeType *> *Edges, std::vector<TetType *> *Tets, bool UsePhon = false);

    virtual void WriteSurfaceData(std::string Filename, int VolumeID)
    {}

    virtual void WriteVolumeData(std::string Filename);
};

}

#endif // ABAQUSEXPORTER_H
