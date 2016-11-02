#include "VTKExporter.h"

namespace voxel2tet
{

VTKExporter::VTKExporter(std::vector<TriangleType *> *Triangles, std::vector<VertexType *> *Vertices,
                         std::vector<EdgeType *> *Edges, std::vector<TetType *> *Tets) :
        Exporter(Triangles, Vertices, Edges, Tets)
{}

vtkSmartPointer<vtkPoints> VTKExporter::SetupVertices()
{
    vtkSmartPointer<vtkPoints> Points = vtkPoints::New();

    for (unsigned int i = 0; i < this->Vertices->size(); i++) {
        VertexType *v = this->Vertices->at(i);
        this->VertexMap[v] = i;
        Points->InsertNextPoint(v->get_c(0), v->get_c(1), v->get_c(2));
    }
    return Points;
}

vtkSmartPointer<vtkCellArray> VTKExporter::SetupTriangles()
{
    vtkSmartPointer<vtkCellArray> Cells = vtkCellArray::New();
    for (unsigned int i = 0; i < this->Triangles->size(); i++) {
        TriangleType *t = this->Triangles->at(i);
        vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
        for (int j = 0; j < 3; j++) {
            VertexType *v = t->Vertices[j];
            int VertexId = this->VertexMap[v];
            triangle->GetPointIds()->SetId(j, VertexId);
        }
        Cells->InsertNextCell(triangle);
    }
    return Cells;
}

vtkSmartPointer<vtkCellArray> VTKExporter::SetupTetrahedrons()
{
    vtkSmartPointer<vtkCellArray> Cells = vtkCellArray::New();
    for (unsigned int i = 0; i < this->Tets->size(); i++) {
        TetType *t = this->Tets->at(i);
        vtkSmartPointer<vtkTetra> tet = vtkSmartPointer<vtkTetra>::New();
        for (int j = 0; j < 4; j++) {
            VertexType *v = t->Vertices[j];
            tet->GetPointIds()->SetId(j, this->VertexMap[v]);
        }
        Cells->InsertNextCell(tet);
    }
    return Cells;
}

vtkSmartPointer<vtkIntArray> VTKExporter::SetupTriangleField(std::string Name, int TriangleType::*FieldPtr)
{
    vtkSmartPointer<vtkIntArray> Triangles = vtkIntArray::New();
    Triangles->SetNumberOfComponents(1);
    Triangles->SetName(Name.c_str());
    for (unsigned int i = 0; i < this->Triangles->size(); i++) {
        TriangleType *v = this->Triangles->at(i);
        int ThisValue = v->*FieldPtr;
        Triangles->InsertNextValue(ThisValue);
    }
    return Triangles;
}

vtkSmartPointer<vtkIntArray> VTKExporter::SetupTetField(std::string Name, int TetType::*FieldPtr)
{
    vtkSmartPointer<vtkIntArray> Triangles = vtkIntArray::New();
    Triangles->SetNumberOfComponents(1);
    Triangles->SetName(Name.c_str());
    for (unsigned int i = 0; i < this->Tets->size(); i++) {
        TetType *t = this->Tets->at(i);
        int ThisValue = t->*FieldPtr;
        Triangles->InsertNextValue(ThisValue);
    }
    return Triangles;
}

vtkSmartPointer<vtkIntArray> VTKExporter::SetupVertexField(std::string Name, int VertexType::*FieldPtr)
{
    vtkSmartPointer<vtkIntArray> VertexIDs = vtkIntArray::New();
    VertexIDs->SetNumberOfComponents(1);
    VertexIDs->SetName(Name.c_str());
    for (unsigned int i = 0; i < this->Vertices->size(); i++) {
        VertexType *v = this->Vertices->at(i);
        int ThisID = v->*FieldPtr;
        VertexIDs->InsertNextValue(ThisID);
    }
    return VertexIDs;
}

vtkSmartPointer<vtkFloatArray> VTKExporter::SetupVertexField(std::string Name, double ( VertexType::*FieldPtr ))
{
    vtkSmartPointer<vtkFloatArray> VertexIDs = vtkFloatArray::New();
    VertexIDs->SetNumberOfComponents(1);
    VertexIDs->SetName(Name.c_str());
    for (unsigned int i = 0; i < this->Vertices->size(); i++) {
        VertexType *v = this->Vertices->at(i);
        VertexIDs->InsertNextValue(v->*FieldPtr);
    }
    return VertexIDs;
}

void VTKExporter::WriteSurfaceData(std::string Filename)
{
    vtkSmartPointer<vtkPoints> Points = this->SetupVertices();
    vtkSmartPointer<vtkCellArray> TriangleArrays = this->SetupTriangles();

    vtkSmartPointer<vtkPolyData> PolyData = vtkSmartPointer<vtkPolyData>::New();
    PolyData->SetPoints(Points);

    // ID
    vtkSmartPointer<vtkIntArray> VertexID = SetupVertexField("Vertex ID", &VertexType::ID);
    PolyData->GetPointData()->AddArray(VertexID);

    // tag
    vtkSmartPointer<vtkIntArray> VertexTag = SetupVertexField("Vertex tag", &VertexType::tag);
    PolyData->GetPointData()->AddArray(VertexTag);

    // Error
    vtkSmartPointer<vtkFloatArray> VertexError = SetupVertexField("Coarsen error", &VertexType::error);
    PolyData->GetPointData()->AddArray(VertexError);

    // Triangles
    PolyData->SetPolys(TriangleArrays);

    // Interface ID
    vtkSmartPointer<vtkIntArray> InterfaceID = SetupTriangleField("Interface ID", &TriangleType::InterfaceID);
    PolyData->GetCellData()->AddArray(InterfaceID);

    // Triangle ID
    vtkSmartPointer<vtkIntArray> TriangleID = SetupTriangleField("Triangle ID", &TriangleType::ID);
    PolyData->GetCellData()->AddArray(TriangleID);

    // Positive normal phase
    vtkSmartPointer<vtkIntArray> TrianglePosNormalPhase = SetupTriangleField("PosNormalMatID",
                                                                             &TriangleType::PosNormalMatID);
    PolyData->GetCellData()->AddArray(TrianglePosNormalPhase);

    // Negative normal phase
    vtkSmartPointer<vtkIntArray> TriangleNegNormalPhase = SetupTriangleField("NegNormalMatID",
                                                                             &TriangleType::NegNormalMatID);
    PolyData->GetCellData()->AddArray(TriangleNegNormalPhase);

    vtkSmartPointer<vtkXMLPolyDataWriter> XMLWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    XMLWriter->SetFileName(Filename.c_str());
#if VTK_MAJOR_VERSION <= 5
    XMLWriter->SetInput(PolyData);
#else
    XMLWriter->SetInputData(PolyData);
#endif
    XMLWriter->Write();
}

void VTKExporter::WriteVolumeData(std::string Filename)
{
    vtkSmartPointer<vtkPoints> Points = this->SetupVertices();
    vtkSmartPointer<vtkCellArray> TetArrays = this->SetupTetrahedrons();

    vtkSmartPointer<vtkUnstructuredGrid> UnstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    UnstructuredGrid->SetPoints(Points);

    UnstructuredGrid->SetCells(VTK_TETRA, TetArrays);

    vtkSmartPointer<vtkIntArray> TetID = SetupTetField("Tet ID", &TetType::ID);
    UnstructuredGrid->GetCellData()->AddArray(TetID);

    vtkSmartPointer<vtkIntArray> MatID = SetupTetField("Mat ID", &TetType::MaterialID);
    UnstructuredGrid->GetCellData()->AddArray(MatID);

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> XMLWriter = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    XMLWriter->SetFileName(Filename.c_str());
#if VTK_MAJOR_VERSION <= 5
    XMLWriter->SetInput(UnstructuredGrid);
#else
    XMLWriter->SetInputData(UnstructuredGrid);
#endif
    XMLWriter->SetDataModeToAscii();
    XMLWriter->Write();
}
}
