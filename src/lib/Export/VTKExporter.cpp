#include "VTKExporter.h"

namespace voxel2tet
{

VTKExporter::VTKExporter()
{

}

VTKExporter::VTKExporter(std :: vector <TriangleType*> *Triangles, std :: vector <VertexType*> *Vertices, std :: vector <EdgeType*> *Edges, std :: vector <TetType*> *Tets) :
    Exporter (Triangles, Vertices, Edges, Tets)
{
    //LOG("Create VTKExporter object\n",0);
}

vtkPoints* VTKExporter :: SetupVertices()
{
    vtkPoints *Points = vtkPoints::New();

    for (unsigned int i=0; i<this->Vertices->size(); i++) {
        VertexType *v=this->Vertices->at(i);
        this->VertexMap[v]=i;
        Points->InsertNextPoint(v->get_c(0), v->get_c(1), v->get_c(2) );
    }
    return Points;
}

vtkCellArray *VTKExporter::SetupTriangles()
{
    vtkCellArray *Cells = vtkCellArray::New();
    for (unsigned int i=0; i<this->Triangles->size(); i++) {
        TriangleType *t=this->Triangles->at(i);
        vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
        for (int j=0; j<3; j++) {
            VertexType *v=t->Vertices[j];
            int VertexId = this->VertexMap[v];
            triangle->GetPointIds()->SetId ( j, VertexId );
        }
        Cells->InsertNextCell(triangle);
    }
    return Cells;
}

vtkCellArray *VTKExporter :: SetupTetrahedrons()
{
    vtkCellArray *Cells = vtkCellArray::New();
    for (unsigned int i=0; i<this->Tets->size(); i++) {
        TetType *t = this->Tets->at(i);
        vtkSmartPointer<vtkTetra> tet = vtkSmartPointer<vtkTetra>::New();
        for (int j=0; j<4; j++) {
            VertexType *v=t->Vertices[j];
            tet->GetPointIds()->SetId (j,  this->VertexMap[v]);
        }
        Cells->InsertNextCell(tet);
    }
    return Cells;
}

vtkIntArray *VTKExporter :: SetupInterfaceIDs()
{
    vtkIntArray *InterfaceIDs = vtkIntArray::New();
    InterfaceIDs->SetNumberOfComponents(1);
    InterfaceIDs->SetName("Surface ID");
    for (unsigned int i=0; i<this->Triangles->size(); i++) {
        TriangleType *t=this->Triangles->at(i);
        InterfaceIDs->InsertNextTuple1(t->InterfaceID);
    }
    return InterfaceIDs;
}


vtkIntArray *VTKExporter :: SetupTriangleIDs()
{
    vtkIntArray *TriangleIDs = vtkIntArray::New();
    TriangleIDs->SetNumberOfComponents(1);
    TriangleIDs->SetName("ID");
    for (unsigned int i=0; i<this->Triangles->size(); i++) {
        TriangleType *t=this->Triangles->at(i);
        TriangleIDs->InsertNextTuple1(t->ID);
    }
    return TriangleIDs;
}

vtkIntArray *VTKExporter :: SetupTrianglePosNormal()
{
    vtkIntArray *TrianglePosNormal = vtkIntArray::New();
    TrianglePosNormal->SetNumberOfComponents(1);
    TrianglePosNormal->SetName("PosNormalPhaseID");
    for (unsigned int i=0; i<this->Triangles->size(); i++) {
        TriangleType *t=this->Triangles->at(i);
        TrianglePosNormal->InsertNextTuple1(t->PosNormalMatID);
    }
    return TrianglePosNormal;
}

vtkIntArray *VTKExporter :: SetupTetIDs()
{
    vtkIntArray *TetIDs = vtkIntArray::New();
    TetIDs->SetNumberOfComponents(1);
    TetIDs->SetName("ID");
    for (unsigned int i=0; i<this->Tets->size(); i++) {
        TetType *t=this->Tets->at(i);
        TetIDs->InsertNextTuple1(t->ID);
    }
    return TetIDs;
}

vtkIntArray *VTKExporter :: SetupTetPhaseID()
{
    vtkIntArray *TetPhaseIDs = vtkIntArray::New();
    TetPhaseIDs->SetNumberOfComponents(1);
    TetPhaseIDs->SetName("PhaseID");
    for (unsigned int i=0; i<this->Tets->size(); i++) {
        TetType *t=this->Tets->at(i);
        TetPhaseIDs->InsertNextTuple1(t->MaterialID);
    }
    return TetPhaseIDs;
}

vtkIntArray *VTKExporter :: SetupVertexIDs()
{
    vtkIntArray *VertexIDs = vtkIntArray::New();
    VertexIDs->SetNumberOfComponents(1);
    VertexIDs->SetName("Vertex ID");
    for (unsigned int i=0; i<this->Vertices->size(); i++) {
        VertexType *v=this->Vertices->at(i);
        VertexIDs->InsertNextValue(v->ID);
    }
    return VertexIDs;
}

vtkIntArray *VTKExporter :: SetupVertexTags()
{
    vtkIntArray *VertexIDs = vtkIntArray::New();
    VertexIDs->SetName("Tag");
    for (unsigned int i=0; i<this->Vertices->size(); i++) {
        VertexType *v=this->Vertices->at(i);
        VertexIDs->InsertNextValue(v->tag);
    }
    return VertexIDs;
}

void VTKExporter :: WriteSurfaceData(std::string Filename)
{

    vtkSmartPointer<vtkPoints> Points;
    vtkSmartPointer<vtkCellArray> TriangleArrays;
    Points.TakeReference(this->SetupVertices());
    TriangleArrays.TakeReference(this->SetupTriangles());

    vtkSmartPointer<vtkPolyData> PolyData = vtkSmartPointer<vtkPolyData>::New();
    PolyData->SetPoints (Points);

    // ID
    vtkSmartPointer<vtkIntArray> VertexID;
    VertexID.TakeReference(this->SetupVertexIDs());

    vtkPointData *pd;
    pd = PolyData->GetPointData();
    pd->SetScalars(VertexID);
    PolyData->Modified();

    // tag
/*    vtkSmartPointer<vtkIntArray> VertexTag;
    VertexTag.TakeReference(this->SetupVertexTags());

    vtkPointData *pdtag;
    pdtag = PolyData->GetPointData();
    pdtag->SetScalars(VertexTag);
    PolyData->Modified();*/

    // Triangles
    PolyData->SetPolys( TriangleArrays );

    // Interface ID
    vtkSmartPointer<vtkIntArray> InterfaceID;
    InterfaceID.TakeReference(this->SetupInterfaceIDs());
    PolyData->GetCellData()->AddArray(InterfaceID);

    PolyData->Modified();

    // Triangle ID
    vtkSmartPointer<vtkIntArray> TriangleID;
    TriangleID.TakeReference(this->SetupTriangleIDs());
    PolyData->GetCellData()->AddArray(TriangleID);

    PolyData->Modified();

    // Positive normal phase
    vtkSmartPointer<vtkIntArray> TrianglePosNormalPhase;
    TrianglePosNormalPhase.TakeReference(this->SetupTrianglePosNormal());
    PolyData->GetCellData()->AddArray(TrianglePosNormalPhase);

    PolyData->Modified();


    vtkSmartPointer<vtkXMLPolyDataWriter> XMLWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    XMLWriter->SetFileName(Filename.c_str());
#if VTK_MAJOR_VERSION <= 5
    XMLWriter->SetInput(PolyData);
#else
    XMLWriter->SetInputData(PolyData);
#endif
    XMLWriter->Write();

}

void VTKExporter :: WriteVolumeData(std::string Filename)
{

    vtkSmartPointer<vtkPoints> Points;
    vtkSmartPointer<vtkCellArray> TetArrays;
    Points.TakeReference(this->SetupVertices());
    TetArrays.TakeReference(this->SetupTetrahedrons());

    vtkSmartPointer<vtkUnstructuredGrid> UnstructuredGrid= vtkSmartPointer<vtkUnstructuredGrid>::New();
    UnstructuredGrid->SetPoints (Points);
    UnstructuredGrid->SetCells(VTK_TETRA, TetArrays );

    vtkSmartPointer<vtkIntArray> RegionID;
    RegionID.TakeReference(this->SetupTetIDs());
    UnstructuredGrid->GetCellData()->AddArray(RegionID);

    vtkSmartPointer<vtkIntArray> PhaseID;
    PhaseID.TakeReference(this->SetupTetPhaseID());
    UnstructuredGrid->GetCellData()->AddArray(PhaseID);

    UnstructuredGrid->Modified();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> XMLWriter = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    XMLWriter->SetFileName(Filename.c_str());
#if VTK_MAJOR_VERSION <= 5
    XMLWriter->SetInput(UnstructuredGrid);
#else
    XMLWriter->SetInputData(UnstructuredGrid);
#endif
    XMLWriter->Write();

}

}

