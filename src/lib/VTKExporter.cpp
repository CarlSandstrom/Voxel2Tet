#include "VTKExporter.h"

namespace voxel2tet
{

VTKExporter::VTKExporter()
{

}

VTKExporter::VTKExporter(MeshData *mesh) : Exporter (mesh)
{
    log("Create VTKExporter object\n",0);
}

vtkPoints* VTKExporter :: SetupVertices()
{
    vtkPoints *Points = vtkPoints::New();

    for (unsigned int i=0; i<this->Vertices->size(); i++) {
        VertexType *v=this->Vertices->at(i);
        this->VertexMap[v]=i;
        Points->InsertNextPoint(v->c[0], v->c[1], v->c[2] );
    }
    return Points;
}

vtkCellArray *VTKExporter::SetupTriangles()
{
    vtkCellArray *Cells = vtkCellArray::New();
    for (unsigned int i=0; i<this->Triangles->size(); i++) {
        log("Setup triangle %u@(%p)\n", i, this->Triangles->at(i));
        TriangleType *t=this->Triangles->at(i);
        vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
        for (int j=0; j<3; j++) {
            VertexType *v=t->Vertices[j];
            int VertexId = this->VertexMap[v];
            log("Add vertex %u@(%p) to triangle\n", VertexId, this->Vertices->at(VertexId));
            triangle->GetPointIds()->SetId ( j, VertexId );
        }
        Cells->InsertNextCell(triangle);
    }
    return Cells;
}

vtkUnsignedCharArray *VTKExporter :: SetupInterfaceIDs()

{
    vtkUnsignedCharArray *InterfaceIDs = vtkUnsignedCharArray::New();
    InterfaceIDs->SetNumberOfComponents(1);
    InterfaceIDs->SetName("Surface ID");
    for (unsigned int i=0; i<this->Triangles->size(); i++) {
        TriangleType *t=this->Triangles->at(i);
        InterfaceIDs->InsertNextTuple1(t->InterfaceID);
    }
    return InterfaceIDs;
}

void VTKExporter :: WriteData(std::string Filename)
{

    vtkSmartPointer<vtkPoints> Points;
    vtkSmartPointer<vtkCellArray> CellArrays;
    Points.TakeReference(this->SetupVertices());
    CellArrays.TakeReference(this->SetupTriangles());

    vtkSmartPointer<vtkPolyData> PolyData = vtkSmartPointer<vtkPolyData>::New();
    PolyData->SetPoints (Points);
    PolyData->SetPolys( CellArrays );

    vtkSmartPointer<vtkUnsignedCharArray> InterfaceID;
    InterfaceID.TakeReference(this->SetupInterfaceIDs());
    PolyData->GetCellData()->SetScalars(InterfaceID);
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


}
