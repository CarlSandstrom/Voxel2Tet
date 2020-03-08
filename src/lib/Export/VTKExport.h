#ifndef VTKEXPORT_H
#define VTKEXPORT_H

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
#include <vtkFloatArray.h>

namespace voxel2tet
{

/**
 * @brief The VTKExporter class provides functionality to export to ParaViews binary .vtu format for volumes and .vtp format for surfaces.This is used
 * mainly for visualization.
 */
class VTKExporter : public Exporter
{
private:
    std::map<VertexType *, int> VertexMap;

    vtkSmartPointer<vtkPoints> SetupVertices();

    vtkSmartPointer<vtkCellArray> SetupTriangles();

    vtkSmartPointer<vtkCellArray> SetupTetrahedrons();

    vtkSmartPointer<vtkIntArray> SetupVertexField(std::string Name, int VertexType::*FieldPtr );

    vtkSmartPointer<vtkFloatArray> SetupVertexField(std::string Name, double ( VertexType::*FieldPtr ));

    vtkSmartPointer<vtkIntArray> SetupTriangleField(std::string Name, int TriangleType::*FieldPtr);

    vtkSmartPointer<vtkIntArray> SetupTetField(std::string Name, int TetType::*FieldPtr);

public:
    /**
     * @copydoc Exporter::Exporter
     */
    VTKExporter(std::vector<TriangleType *> *Triangles, std::vector<VertexType *> *Vertices,
                std::vector<EdgeType *> *Edges, std::vector<TetType *> *Tets);

    void WriteSurfaceData(std::string Filename, int VolumeID=-1);

    void WriteVolumeData(std::string Filename);
};
}

#endif // VTKEXPORT_H
