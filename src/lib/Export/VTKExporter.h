#ifndef VTKEXPORTER_H
#define VTKEXPORTER_H

#include <string>
#include <map>

#include "Exporter.h"
#include "MeshComponents.h"
#include "MiscFunctions.h"

#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkTetra.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>

namespace voxel2tet
{

class VTKExporter : public Exporter
{
private:
    std::map <VertexType*, int> VertexMap;

    vtkPoints *SetupVertices();
    vtkCellArray *SetupTriangles();
    vtkCellArray *SetupTetrahedrons();
    vtkIntArray *SetupInterfaceIDs();
    vtkIntArray *SetupTriangleIDs();
    vtkIntArray *SetupTrianglePosNormal();
    vtkIntArray *SetupTetIDs();
    vtkIntArray *SetupVertexIDs();
    vtkIntArray *SetupVertexTags();
public:
    VTKExporter();
    VTKExporter(std :: vector <TriangleType*> *Triangles, std :: vector <VertexType*> *Vertices, std :: vector <EdgeType*> *Edges, std::vector<TetType *> *Tets);
    void WriteSurfaceData(std::string Filename);
    void WriteVolumeData(std::string Filename);
};

}

#endif // VTKEXPORTER_H
