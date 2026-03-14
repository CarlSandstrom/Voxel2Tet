#include "TetPipeline.h"

#include "Importer.h"
#include "MeshManipulations.h"
#include "TetGenCaller.h"
#include "Volume.h"

#include <cstdlib>
#include <map>

namespace voxel2tet
{
TetPipeline::TetPipeline(Importer *Imp, std::unique_ptr<MeshManipulations> *Mesh, std::vector<Volume *> *Volumes,
                         double eps)
{
    this->Imp = Imp;
    this->Mesh = Mesh;
    this->Volumes = Volumes;
    this->eps = eps;
}

void TetPipeline::Tetrahedralize() const
{
    TetGenCaller Generator;
    Generator.Mesh = this->Mesh->get();
#if TEST_MESH_BETWEEN_STEPS_TETGEN == 1
    Generator.TestMesh();
#endif

    MeshData *TetMesh = Generator.Execute();
    MeshManipulations *NewMesh = static_cast<MeshManipulations *>(TetMesh);

    std::map<int, int> Tetgen2Self;

    for (TetType *t : NewMesh->Tets) {
        if (Tetgen2Self.find(t->MaterialID) == Tetgen2Self.end()) {
            std::array<double, 3> cm = t->GiveCenterOfMass();
            Volume *v = this->FindVolumeContainingPoint(cm);
            Tetgen2Self[t->MaterialID] = v->Phase;
        }

        t->MaterialID = Tetgen2Self[t->MaterialID];
    }

    this->Mesh->reset(NewMesh);

    this->UpdateFixed();
}

Volume *TetPipeline::FindVolumeContainingPoint(std::array<double, 3> P) const
{
    for (Volume *v : *this->Volumes) {
        if (v->IsPointInside(P)) {
            return v;
        }
    }

    return nullptr;
}

void TetPipeline::UpdateFixed() const
{
    for (VertexType *v : this->Mesh->get()->Vertices) {
        for (int j = 0; j < 3; j++) {
            if ((v->get_c(j) > (this->Imp->GiveBoundingBox().maxvalues[j] - this->eps)) |
                    (v->get_c(j) < (this->Imp->GiveBoundingBox().minvalues[j] + this->eps))) {
                v->Fixed[j] = true;
            } else {
                v->Fixed[j] = false;
            }
        }

        if (v->IsFixedVertex()) {
            for (int i = 0; i < 3; i++) {
                v->Fixed[i] = true;
            }
        }
    }
}
}
