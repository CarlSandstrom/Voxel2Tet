#ifdef OPENMP
#include <omp.h>
#endif
#include"Smoother.h"

namespace voxel2tet
{

double Compute_c(double l, double alpha)
{
    double c=l; // Initial guess
    double R=exp(pow(l/c, alpha))-1-l;
    double err = fabs(R);

    while (err > 1e-8) {
        double tangent = -exp(pow(l/c, alpha))*alpha*pow(l/c, alpha-1)*l*pow(c,-2);
        double deltac = -1/tangent*R;
        c=c+deltac;
        R=exp(pow(l/c, alpha))-1-l;
        err = fabs(R);
    }
    return c;
}

arma::vec ComputeOutOfBalance(std::vector<std::array<double, 3> > ConnectionCoords, arma::vec xc, arma::vec x0, double alpha, double c)
{

    arma::vec F = {0., 0., 0.};

    // Compute nonlinear part of force
    arma::vec n0;

    double d0 = arma::norm(x0-xc);
    if (d0<1e-8) {
        n0 = {0.,0.,0.};
    } else {
        n0 = (x0-xc)/d0;
    }

    F = (exp(pow(d0/c, alpha))-1) * n0;

    // Compute linear part of force
    for (unsigned int i=0; i<ConnectionCoords.size(); i++) {
        arma::vec xi = {ConnectionCoords.at(i)[0], ConnectionCoords.at(i)[1], ConnectionCoords.at(i)[2]};
        arma::vec nj;
        double dj = arma::norm(xi-xc);
        if (dj<1e-8) {
            nj = {0., 0., 0.,};
        } else {
            nj = (xi-xc)/dj;
        }
        F = F + dj*nj;
    }

    return F;

}

arma::mat  ComputeNumericalTangent(std::vector<std::array<double, 3> > ConnectionCoords, arma::vec xc, arma::vec x0, double alpha, double c)
{

    double eps=1e-10;
    arma::mat Tangent=arma::zeros<arma::mat>(3,3);

    arma::vec Fval=ComputeOutOfBalance(ConnectionCoords, xc, x0, alpha, c);

    for (int i=0; i<3; i++) {
        arma::vec xi = xc;
        xi[i] = xi[i] + eps;
        arma::vec Fvali=ComputeOutOfBalance(ConnectionCoords, xi, x0, alpha, c);
        arma::vec dF = (Fvali-Fval)/eps;
        Tangent.col(i) = dF;
    }

    return Tangent;

}

arma::mat ComputeAnalyticalTangent(std::vector<std::array<double, 3> > ConnectionCoords, arma::vec xc, arma::vec x0, double alpha, double c)
{
    arma::mat Tangent=arma::zeros<arma::mat>(3,3);

    arma::vec a0 = x0-xc;
    double d0 = arma::norm(a0);
    arma::vec n0;

    // Non-linear part
    // We run into numerical trouble if d0=0. However, in the case of d0=0, everyting nonlinear is zero...
    arma::vec Dexp;
    arma::mat Dn0;
    if (d0 < 1e-8) {
        n0 = {0.,0.,0.};
        Dexp = {0.,0.,0.};
        Dn0 = -arma::zeros<arma::mat>(3,3);

    } else {
        n0 = a0/d0;
        Dexp = -std::exp(std::pow(d0/c, alpha))*alpha/c*pow(d0/c, alpha-1)*n0;
        Dn0 = -arma::eye<arma::mat>(3,3)/d0 + a0*a0.t()/pow(d0,3);
    }
    Tangent = Tangent + n0*Dexp.t() + Dn0 * ( std::exp(std::pow(d0/c, alpha)) - 1);


    // Linear part
    Tangent = Tangent - arma::eye<arma::mat>(3,3)*ConnectionCoords.size();

    return Tangent;
}

arma::mat ComputeAnalyticalTangentGlobal(std::vector<std::array<double, 3> > ConnectionCoords, arma::vec xc, arma::vec x0, double alpha, double c)
{

    arma::mat Tangent=arma::zeros<arma::mat>(3,ConnectionCoords.size()*3+3);
    arma::mat TangentSelf=arma::zeros<arma::mat>(3,3);

    arma::vec a0 = x0-xc;
    double d0 = arma::norm(a0);
    arma::vec n0;

    // delta x_i part

    // Non-linear part
    // We run into numerical trouble if d0=0. However, in the case of d0=0, everyting nonlinear is zero...
    arma::vec Dexp;
    arma::mat Dn0;
    if (d0 < 1e-8) {
        n0 = {0.,0.,0.};
        Dexp = {0.,0.,0.};
        Dn0 = -arma::zeros<arma::mat>(3,3);

    } else {
        n0 = a0/d0;
        Dexp = -std::exp(std::pow(d0/c, alpha))*alpha/c*pow(d0/c, alpha-1)*n0;
        Dn0 = -arma::eye<arma::mat>(3,3)/d0 + a0*a0.t()/pow(d0,3);
    }
    TangentSelf = n0*Dexp.t() + Dn0 * ( std::exp(std::pow(d0/c, alpha)) - 1);

    // Linear part
    TangentSelf = TangentSelf - arma::eye<arma::mat>(3,3)*ConnectionCoords.size();

    // Assemble self part to tangent
    Tangent(arma::span(0,2),arma::span(0,2)) = TangentSelf;

    // delta x_j part
    for (unsigned int i=0; i<ConnectionCoords.size(); i++) {
        for (int j=0; j<3; j++) {
            Tangent(j, 3*(i+1)+j) = 1.0;// /ConnectionCoords.size();
        }
    }
    return Tangent;
}

void SpringSmoothGlobal(std::vector<VertexType*> Vertices, std::vector<bool> Fixed, std::vector<std::vector<VertexType*>> Connections, double c, double alpha, double charlength, bool Automatic_c, voxel2tet::MeshData *Mesh)
{
    int MAX_ITER_COUNT=1000;

    //Fixed.at(0)=true;

    if (Automatic_c) {
        // c = Compute_c(charlength);
    }

    // Create vectors for current and previous positions
    std::vector<std::array<double, 3>> CurrentPositions;
    std::vector<std::array<int, 3>> dofids;

    int ndof = 0;

    for (unsigned int i=0; i<Vertices.size(); i++) {
        std::array<double, 3> cp;
        std::array<int, 3> dofid = {{-1, -1, -1}};
        for (int j=0; j<3; j++) {
            cp.at(j) = Vertices.at(i)->get_c(j);
            if (!Fixed[i]) {
                if (!Vertices.at(i)->Fixed[j]) {
                    dofid[j]=ndof;
                    ndof++;
                }
            }
        }
        CurrentPositions.push_back(cp);
        dofids.push_back(dofid);
    }

    // Create vector-vector for accessing neightbours
    std::vector <std::vector <int>> ConnectionVertexIndex;
    for (std::vector<VertexType*> n: Connections) {
        ConnectionVertexIndex.push_back({});
        for (VertexType *v: n) {
            int VertexIndex = std::distance(Vertices.begin(), std::find(Vertices.begin(), Vertices.end(), v));
            ConnectionVertexIndex.at(ConnectionVertexIndex.size()-1).push_back(VertexIndex);
        }
    }

#if EXPORT_SMOOTHING_ANIMATION == 1
    std::ostringstream FileName;
    if (Mesh!=NULL) {
        FileName << "/tmp/Smoothing" << 0 << ".vtp";
        Mesh->ExportSurface( FileName.str(), FT_VTK );
    }
#endif

    int itercount = 0;
    double err = 1e8;

    arma::sp_mat Kff(ndof, ndof);
    arma::vec Rf = arma::ones<arma::vec>(ndof);

    while ((itercount < MAX_ITER_COUNT) & (err>1e-6)) {

        Kff.zeros();

        for (size_t i=0; i<Vertices.size(); i++) {
            if (!Fixed[i]) {

                std::vector<std::array<double, 3>> ConnectionCoords;
                std::vector<VertexType*> MyConnections = Connections.at(i);

                for (unsigned k=0; k<MyConnections.size(); k++) {
                    ConnectionCoords.push_back({CurrentPositions.at(ConnectionVertexIndex.at(i).at(k))[0],
                                                CurrentPositions.at(ConnectionVertexIndex.at(i).at(k))[1],
                                                CurrentPositions.at(ConnectionVertexIndex.at(i).at(k))[2]});
                }

                arma::vec xc = {CurrentPositions.at(i)[0], CurrentPositions.at(i)[1], CurrentPositions.at(i)[2]};
                arma::vec x0 = {Vertices.at(i)->get_c(0), Vertices.at(i)->get_c(1), Vertices.at(i)->get_c(2)};

                // Assemble out-of-balance vector
                arma::vec R = ComputeOutOfBalance(ConnectionCoords, xc, x0, alpha, c);

                for (int j=0; j<3; j++) {
                    if (dofids[i][j]!=-1) { // Not stationary
                        Rf[dofids[i][j]] = R[j];
                    }
                }

                // Assemble to tangent
                std::vector <int> ConnectionIndices = ConnectionVertexIndex.at(i);
                arma::mat T = ComputeAnalyticalTangent(ConnectionCoords, xc, x0, alpha, c);

                std::vector <int> Rows = {dofids[i][0], dofids[i][1], dofids[i][2]};
                std::vector <int> Cols = {dofids[i][0], dofids[i][1], dofids[i][2]};

                for (int index: ConnectionIndices) {
                    for (int j=0; j<3; j++) {
                        Cols.push_back(dofids[index][j]);
                    }
                }

                for (int j=0; j<3; j++) {
                    if (Rows[j]!=-1) {
                        for (size_t k=0; k<Cols.size(); k++) {
                            if (Cols[k]!=-1) {
                                Kff(Rows[j], Cols[k]) = Kff(Rows[j], Cols[k]) + T(j, k);
                            }
                        }
                    }
                }
            }
        }

        //Kff.print("Kff = ");

        err = arma::norm(Rf);
        //Rf.print("Rf = ");

        arma::vec delta;// = arma::spsolve(Kff, -Rf);
        delta.zeros(ndof);
        //delta.print("delta:");


        // Solve system using CG algorihtm
        int maxiter = 10000000;

        arma::vec r = -Rf - Kff*delta;
        arma::vec rn;
        arma::vec p = r;

        for (int i = 0; i < maxiter; i++) {
            arma::vec Ap = Kff*p;

            double alpha = arma::as_scalar(r.t()*r)/arma::as_scalar(p.t()*Ap);
            delta = delta + alpha*p;
            rn = r - alpha*Ap;
            double residual = arma::norm(r);
            if (residual<1e-8) break;
            double beta = arma::as_scalar(rn.t()*rn)/arma::as_scalar(r.t()*r);
            p=rn+beta*p;
            r = rn;
        }

        if (err>1e-10) {
            //delta.print("delta = ");
            double maxdelta = arma::max(delta);
            if (fabs(maxdelta) > 0.001) {
                delta = delta * 0.001/maxdelta;
            }
        }

        // Update current positions
        for (size_t i=0; i<Vertices.size(); i++) {
            for (int j=0; j<3; j++) {
                if (dofids[i][j]!=-1) {
                    CurrentPositions.at(i)[j] = CurrentPositions.at(i)[j] + delta[dofids[i][j]];
                }
            }
        }

        STATUS("Iteration %i end with err=%f\n", itercount, err);

        itercount ++;

#if EXPORT_SMOOTHING_ANIMATION == 1
        // ************************** DEBUG STUFF
        // Update vertices
        for (unsigned int i=0; i<Vertices.size(); i++) {
            VertexType *v=Vertices.at(i);
            for (int j=0; j<3; j++) {
                if (!v->Fixed[j]) {
                    v->set_c(CurrentPositions.at(i)[j], j);
                }
            }
        }

        if (Mesh!=NULL) {
            FileName.str(""); FileName.clear();
            FileName << "/tmp/Smoothing" << itercount++ << ".vtp";
            Mesh->ExportSurface(FileName.str(), FT_VTK);
        }
        // ************************** /DEBUG STUFF
#endif

#if TEST_MESH_FOR_EACH_SMOOTHING
        TetGenCaller Tetgen;
        Tetgen.Mesh = Mesh;
        Tetgen.TestMesh();
#endif
    }

    if (itercount > 999) {
        STATUS("WARNING: Smoothing did not converge\n", 0);
    }

    // Update vertices
    for (unsigned int i=0; i<Vertices.size(); i++) {
        for (int j=0; j<3; j++) {
            Vertices.at(i)->set_c(CurrentPositions.at(i)[j], j); // TODO: Use array to improve performance
        }
    }

}

std::vector<std::pair<TriangleType *, TriangleType *>> CheckPenetration(std::vector<VertexType *> *Vertices, MeshData *Mesh)
{

    std::vector<std::pair<TriangleType *, TriangleType *>> IntersectingTriangles;
    std::vector<TriangleType *> Triangles;

    for (VertexType *v: *Vertices) {
        for(TriangleType *t: v->Triangles) {
            Triangles.push_back(t);
        }
    }
    std::sort(Triangles.begin(), Triangles.end());
    Triangles.erase(std::unique(Triangles.begin(), Triangles.end()), Triangles.end());
    // return IntersectingTriangles;

    for (TriangleType *t1: Triangles) {

        std::array<double, 3> c = t1->GiveCenterOfMass();
        double d = t1->GiveLongestEdgeLength();

        std::vector<TriangleType *> NearTriangles = Mesh->GetTrianglesAround(c, d*2);

        for (TriangleType *t2: NearTriangles) {
            if (t1!=t2){
                if (Mesh->CheckTrianglePenetration(t1, t2)) {
                    Mesh->CheckTrianglePenetration(t1, t2);
                    IntersectingTriangles.push_back({t1, t2});
                    LOG("Triangles %u and %u intersect!\n", t1->ID, t2->ID);
                    LOG("t1(%u): (%f, %f, %f), (%f, %f, %f), (%f, %f, %f)\n", t1->ID, t1->Vertices.at(0)->get_c(0), t1->Vertices.at(0)->get_c(1), t1->Vertices.at(0)->get_c(2),
                        t1->Vertices.at(1)->get_c(0), t1->Vertices.at(1)->get_c(1), t1->Vertices.at(1)->get_c(2),
                        t1->Vertices.at(2)->get_c(0), t1->Vertices.at(2)->get_c(1), t1->Vertices.at(2)->get_c(2));
                    LOG("t2(%u): (%f, %f, %f), (%f, %f, %f), (%f, %f, %f)\n", t2->ID, t2->Vertices.at(0)->get_c(0), t2->Vertices.at(0)->get_c(1), t2->Vertices.at(0)->get_c(2),
                        t2->Vertices.at(1)->get_c(0), t2->Vertices.at(1)->get_c(1), t2->Vertices.at(1)->get_c(2),
                        t2->Vertices.at(2)->get_c(0), t2->Vertices.at(2)->get_c(1), t2->Vertices.at(2)->get_c(2));
/*
                    std::vector<VertexType*>  StiffenVertices;

                    for (VertexType *v: t1->Vertices) StiffenVertices.push_back(v);
                    for (VertexType *v: t2->Vertices) StiffenVertices.push_back(v);

                    std::sort(StiffenVertices.begin(), StiffenVertices.end());
                    StiffenVertices.erase(std::unique(StiffenVertices.begin(), StiffenVertices.end()), StiffenVertices.end());

                    for (VertexType *v: StiffenVertices) {
                        //v->c_constant = v->c_constant * 0.75;
                        printf("(%u: c)=%f, ", v->ID, v->c_constant);
                    }
                    LOG("\n", 0);*/
                }
            }
        }
    }
    return IntersectingTriangles;
}

void SpringSmooth (std::vector<VertexType*> Vertices, std::vector<bool> Fixed, std::vector<std::vector<VertexType*>> Connections,
                   double c, double alpha, double charlength, bool Automatic_c, voxel2tet::MeshData *Mesh)
{
    int MAX_ITER_COUNT=100000;

    // Create vectors for current and previous positions
    std::vector<std::array<double, 3>> OriginalPositions;
    std::vector<std::array<double, 3>> CurrentPositions;
    std::vector<std::array<double, 3>> PreviousPositions;

    for (unsigned int i=0; i<Vertices.size(); i++) {
        std::array<double, 3> cp;

        for (int j=0; j<3; j++) {
            cp.at(j) = Vertices.at(i)->get_c(j);
        }

        OriginalPositions.push_back(cp);
        CurrentPositions.push_back(cp);
        PreviousPositions.push_back(cp);
    }

    // Create vector-vector for accessing neightbours
    std::vector <std::vector <int>> ConnectionVertexIndex;
    for (std::vector<VertexType*> n: Connections) {
        ConnectionVertexIndex.push_back({});
        for (VertexType *v: n) {
            int VertexIndex = std::distance(Vertices.begin(), std::find(Vertices.begin(), Vertices.end(), v));
            ConnectionVertexIndex.at(ConnectionVertexIndex.size()-1).push_back(VertexIndex);
        }
    }

    // Reset 'c' constant for all vertices
    for (VertexType *v: Vertices) {
        v->c_constant = c;
    }

#if EXPORT_SMOOTHING_ANIMATION == 1
    std::ostringstream FileName;
    if (Mesh!=NULL) {
        FileName << "/tmp/Smoothing" << 0 << ".vtp";
        Mesh->ExportSurface( FileName.str(), FT_VTK );
    }
#endif

#ifdef OPENMP
    threadcount = omp_get_max_threads();
#endif

    bool intersecting = true;
    int intersecting_count=0;

    CheckPenetration(&Vertices, Mesh);

    while (intersecting) {
        intersecting = false;
        int itercount = 0;
        int threadcount = 1;
        double deltamax = 1e8;
        int deltamaxnode;

        while ((itercount < MAX_ITER_COUNT) & (deltamax>(charlength*1e-3))) {

            deltamax=0.0;
            int deltamaxnodes[threadcount];
            double deltamaxvalues[threadcount];

            for (int i=0; i<threadcount; i++) {
                deltamaxnodes[i]=0;
                deltamaxvalues[i]=0.0;
            }

            //#pragma omp parallel default(shared)
            {
                int threadid=0;
#ifdef OPENMP
                threadid=omp_get_thread_num();
#endif

                //#pragma omp for schedule(static, 100)
                for (size_t i=0; i<Vertices.size(); i++) {
                    if (!Fixed[i]) {

                        std::vector<std::array<double, 3>> ConnectionCoords;
                        std::vector<VertexType*> MyConnections = Connections.at(i);

                        for (unsigned k=0; k<MyConnections.size(); k++) {
                            ConnectionCoords.push_back({{PreviousPositions.at(ConnectionVertexIndex.at(i).at(k))[0],
                                                         PreviousPositions.at(ConnectionVertexIndex.at(i).at(k))[1],
                                                         PreviousPositions.at(ConnectionVertexIndex.at(i).at(k))[2]}});
                        }

                        arma::vec xc = {PreviousPositions.at(i)[0], PreviousPositions.at(i)[1], PreviousPositions.at(i)[2]};
                        arma::vec x0 = {OriginalPositions.at(i)[0], OriginalPositions.at(i)[1], OriginalPositions.at(i)[2]};

                        // Compute out-of-balance vector
                        arma::vec R = ComputeOutOfBalance(ConnectionCoords, xc, x0, alpha, Vertices.at(i)->c_constant);
                        double err = arma::norm(R);
                        int iter = 0;

                        // Find equilibrium
                        while ((err>1e-5) & (iter<100000)) {
                            arma::mat T = ComputeAnalyticalTangent(ConnectionCoords, xc, x0, alpha, Vertices.at(i)->c_constant);
                            arma::vec delta = -arma::solve(T, R);
                            xc = xc + 1.*delta;
                            R = ComputeOutOfBalance(ConnectionCoords, xc, x0, alpha, Vertices.at(i)->c_constant);
                            err = arma::norm(R);
                            iter ++;
                        }

                        // If too many iterations, throw an exception and investigate...
                        if (iter>99999) {
                            throw(0);
                        }

                        // Update current position
                        for (int j=0; j<3; j++) {
                            if (!Vertices.at(i)->Fixed[j]) {
                                CurrentPositions.at(i)[j] = xc[j];
                            }
                        }

                        // Update maximum delta
                        arma::vec d={PreviousPositions.at(i)[0]-CurrentPositions.at(i)[0], PreviousPositions.at(i)[1]-CurrentPositions.at(i)[1], PreviousPositions.at(i)[2]-CurrentPositions.at(i)[2]};
                        if (arma::norm(d)>deltamaxvalues[threadid]) {
                            deltamaxvalues[threadid] = arma::norm(d);
                            deltamaxnodes[threadid] = i;
                        }
                    }
                }

                for (unsigned int i=0; i<Vertices.size(); i++) {
                    for (int j=0; j<3; j++) {
                        Vertices.at(i)->set_c(CurrentPositions.at(i)[j], j); // TODO: Use array to improve performance
                    }
                }

                // Update previous positions
                for (size_t i=0; i<Vertices.size(); i++) {
                    for (int j=0; j<3; j++) {
                        if (!Vertices.at(i)->Fixed[j]) {
                            PreviousPositions.at(i)[j] = CurrentPositions.at(i)[j];
                        }
                    }
                }
            }

            deltamaxnode=0;
            for (int i=0; i<threadcount; i++) {
                if (deltamaxvalues[i]>deltamax) {
                    deltamax = deltamaxvalues[i];
                    deltamaxnode = deltamaxnodes[i];
                }
            }

            STATUS("%c[2K\rIteration %u end with deltamax=%f at node %i\r", 27, itercount, deltamax, deltamaxnode);
            fflush(stdout);

            itercount ++;

#if EXPORT_SMOOTHING_ANIMATION == 1
            // ************************** DEBUG STUFF
            // Update vertices
            for (unsigned int i=0; i<Vertices.size(); i++) {
                VertexType *v=Vertices.at(i);
                for (int j=0; j<3; j++) {
                    if (!v->Fixed[j]) {
                        v->set_c(CurrentPositions.at(i)[j], j);
                    }
                }
            }

            if (Mesh!=NULL) {
                FileName.str(""); FileName.clear();
                FileName << "/tmp/Smoothing" << itercount++ << ".vtp";
                Mesh->ExportSurface(FileName.str(), FT_VTK);
            }
            // ************************** /DEBUG STUFF
#endif

#if TEST_MESH_FOR_EACH_SMOOTHING
            TetGenCaller Tetgen;
            Tetgen.Mesh = Mesh;
            Tetgen.TestMesh();
#endif
        }

        // Check for intersecting triangles. If some triangles intersect, stiffen the structure in that area and re-smooth
        std::vector<std::pair<TriangleType *, TriangleType *>> IntersectingTriangles = CheckPenetration(&Vertices, Mesh);
        if (IntersectingTriangles.size()>0) {

            // Pull back nodes until no intersecting triangles remain
            intersecting_count=0;
            while (IntersectingTriangles.size()>0) {

                STATUS("Stiffen and re-smooth surface/edge due to %u intersecting triangles, iteration %u\n", IntersectingTriangles.size(), intersecting_count);

                // Stiffen spring
                std::vector<VertexType *> TriangleVertices;
                for (std::pair<TriangleType *, TriangleType *> p: IntersectingTriangles) {
                    for (VertexType *v: p.first->Vertices) TriangleVertices.push_back(v);
                    for (VertexType *v: p.second->Vertices) TriangleVertices.push_back(v);
                }

                std::sort(TriangleVertices.begin(), TriangleVertices.end());
                TriangleVertices.erase(std::unique(TriangleVertices.begin(), TriangleVertices.end()), TriangleVertices.end());

                for (VertexType *v: TriangleVertices) {
                    arma::vec displacement = {v->get_c(0)-v->originalcoordinates[0], v->get_c(1)-v->originalcoordinates[1], v->get_c(2)-v->originalcoordinates[2]};
                    arma::vec newdisplacement = displacement*.9;

                    for (int i=0; i<3; i++) v->set_c(v->originalcoordinates[i]+newdisplacement[i], i);
                }

                IntersectingTriangles = CheckPenetration(&Vertices, Mesh);
                intersecting_count++;
            }
/*
            intersecting = true;
            STATUS("Stiffen and re-smooth surface/edge due to %u intersecting triangles, iteration %u\n", IntersectingTriangles.size(), intersecting_count);


            if (intersecting_count>50) {
                STATUS("Too many intersection iterations. Just let it be and it might turn out ok anyway due to the mesh coarsening.\n", 0);
                return;
            }

            for (VertexType *v: StiffenVertices) v->c_constant = v->c_constant * .95;

            // Revert to original coordinates
            for (size_t i=0; i<Vertices.size(); i++) {
                for (int j=0; j<3; j++) {
                    Vertices.at(i)->set_c(OriginalPositions.at(i)[j], j);
                }
                PreviousPositions.at(i) = OriginalPositions.at(i);
                CurrentPositions.at(i) = OriginalPositions.at(i);
            }

            intersecting_count++;*/

#if EXPORT_SMOOTHING_ANIMATION == 1
    std::ostringstream FileName;
    if (Mesh!=NULL) {
        FileName << "/tmp/Smoothing" << 0 << ".vtp";
        Mesh->ExportSurface( FileName.str(), FT_VTK );
    }
#endif
        } else {
            intersecting = false;
            // Update vertices
            for (unsigned int i=0; i<Vertices.size(); i++) {
                for (int j=0; j<3; j++) {
                    Vertices.at(i)->set_c(CurrentPositions.at(i)[j], j); // TODO: Use array to improve performance
                }
            }
        }
        STATUS("\n", 0);

    }

}

}
