#ifndef VTKEXPORTER_H
#define VTKEXPORTER_H

#include <string>
#include <map>

#include "Exporter.h"
#include "MeshComponents.h"
#include "MiscFunctions.h"

#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>


namespace voxel2tet
{

class VTKExporter : public Exporter
{
private:
    std::map <VertexType*, int> VertexMap;

    vtkPoints* SetupVertices();
    vtkCellArray* SetupTriangles();
    vtkUnsignedCharArray *SetupInterfaceIDs();
public:
    VTKExporter();
    VTKExporter(std :: vector <TriangleType*> *Triangles, std :: vector <VertexType*> *Vertices, std :: vector <EdgeType*> *Edges);
    void WriteData(std::string Filename);
};

}

#endif // VTKEXPORTER_H
