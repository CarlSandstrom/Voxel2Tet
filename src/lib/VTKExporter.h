#ifndef VTKEXPORTER_H
#define VTKEXPORTER_H

#include <string>
#include <map>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>

#include "Exporter.h"
#include "MeshComponents.h"
#include "MiscFunctions.h"


namespace voxel2tet
{

//class MeshData;

class VTKExporter : public Exporter
{
private:
    std::map <VertexType*, int> VertexMap;

    vtkPoints* SetupVertices();
    vtkCellArray* SetupTriangles();
    vtkUnsignedCharArray *SetupInterfaceIDs();
public:
    VTKExporter();
    VTKExporter(MeshData *mesh);
    void WriteData(std::string Filename);
};

}

#endif // VTKEXPORTER_H
