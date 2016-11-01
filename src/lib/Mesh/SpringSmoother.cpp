#ifdef OPENMP
#include <omp.h>
#endif
#include "SpringSmoother.h"

namespace voxel2tet
{

SpringSmoother::SpringSmoother(double VoxelCharLength, double c, double alpha, double c_factor, bool compute_c) : Smoother ()
{

    this->alpha = alpha;
    this->charlength = VoxelCharLength;
    this->c_factor = c_factor;

    if (compute_c) {
        this->c = this->Compute_c(VoxelCharLength * c_factor, this->alpha);
    } else {
        this->c = c;
    }

}

double SpringSmoother::Compute_c(double l, double alpha)
{
    double c = l; // Initial guess
    double R = exp( pow(l / c, alpha) ) - 1 - l;
    double err = fabs(R);
    int iter = 0;

    while ( err > 1e-8 ) {
        double tangent = -exp( pow(l / c, alpha) ) * alpha * pow(l / c, alpha - 1) * l * pow(c, -2);
        double deltac = -1 / tangent * R;
        c = c + deltac;
        R = exp( pow(l / c, alpha) ) - 1 - l;
        err = fabs(R);
        iter++;
        if (iter > 1000) {
            STATUS("Unable to find a suitable c\n", 0);
            throw(0);
        }
    }

    STATUS("\tUsing alpha=%f, c=%f\n", alpha, c);

    return c;
}

arma :: vec SpringSmoother::ComputeOutOfBalance(std :: vector< arma::vec3 >ConnectionCoords, arma :: vec3 xc, arma :: vec3 x0, double alpha, double c)
{
    arma :: vec F = {
        0., 0., 0.
    };

    // Compute nonlinear part of force
    arma :: vec n0;

    double d0 = arma :: norm(x0 - xc);
    if ( d0 < 1e-8 ) {
        n0 = {
            0., 0., 0.
        };
    } else {
        n0 = ( x0 - xc ) / d0;
    }

    F = ( exp( pow(d0 / c, alpha) ) - 1 ) * n0;

    // Compute linear part of force
    for ( unsigned int i = 0; i < ConnectionCoords.size(); i++ ) {
        arma :: vec xi = {
            ConnectionCoords.at(i) [ 0 ], ConnectionCoords.at(i) [ 1 ], ConnectionCoords.at(i) [ 2 ]
        };
        arma :: vec nj;
        double dj = arma :: norm(xi - xc);
        if ( dj < 1e-8 ) {
            nj = {
                0., 0., 0.,
            };
        } else {
            nj = ( xi - xc ) / dj;
        }
        F = F + dj * nj / ConnectionCoords.size();
    }

    return F;
}

arma :: mat  SpringSmoother::ComputeNumericalTangent(std :: vector< arma::vec3 >ConnectionCoords, arma :: vec xc, arma :: vec x0, double alpha, double c)
{
    double eps = 1e-10;
    arma :: mat Tangent = arma :: zeros< arma :: mat >(3, 3);

    arma :: vec Fval = ComputeOutOfBalance(ConnectionCoords, xc, x0, alpha, c);

    for ( int i = 0; i < 3; i++ ) {
        arma :: vec xi = xc;
        xi [ i ] = xi [ i ] + eps;
        arma :: vec Fvali = ComputeOutOfBalance(ConnectionCoords, xi, x0, alpha, c);
        arma :: vec dF = ( Fvali - Fval ) / eps;
        Tangent.col(i) = dF;
    }

    return Tangent;
}

arma :: mat SpringSmoother::ComputeAnalyticalTangent(std :: vector< arma::vec3 >ConnectionCoords, arma :: vec xc, arma :: vec x0, double alpha, double c)
{
    arma :: mat Tangent = arma :: zeros< arma :: mat >(3, 3);

    arma :: vec a0 = x0 - xc;
    double d0 = arma :: norm(a0);
    arma :: vec n0;

    // Non-linear part
    // We run into numerical trouble if d0=0. However, in the case of d0=0, everyting nonlinear is zero...
    arma :: vec Dexp;
    arma :: mat Dn0;
    if ( d0 < 1e-8 ) {
        n0 = {
            0., 0., 0.
        };
        Dexp = {
            0., 0., 0.
        };
        Dn0 = -arma :: zeros< arma :: mat >(3, 3);
    } else {
        n0 = a0 / d0;
        Dexp = -std :: exp( std :: pow(d0 / c, alpha) ) * alpha / c *pow(d0 / c, alpha - 1) * n0;
        Dn0 = -arma :: eye< arma :: mat >(3, 3) / d0 + a0 *a0.t() / pow(d0, 3);
    }
    Tangent = Tangent + n0 *Dexp.t() + Dn0 * ( std :: exp( std :: pow(d0 / c, alpha) ) - 1 );


    // Linear part
    Tangent = Tangent - arma :: eye< arma :: mat >(3, 3);

    return Tangent;
}

arma :: mat SpringSmoother::ComputeAnalyticalTangentGlobal(std :: vector< arma::vec3 >ConnectionCoords, arma :: vec xc, arma :: vec x0, double alpha, double c)
{
    arma :: mat Tangent = arma :: zeros< arma :: mat >(3, ConnectionCoords.size() * 3 + 3);
    arma :: mat TangentSelf = arma :: zeros< arma :: mat >(3, 3);

    arma :: vec a0 = x0 - xc;
    double d0 = arma :: norm(a0);
    arma :: vec n0;

    // delta x_i part

    // Non-linear part
    // We run into numerical trouble if d0=0. However, in the case of d0=0, everyting nonlinear is zero...
    arma :: vec Dexp;
    arma :: mat Dn0;
    if ( d0 < 1e-8 ) {
        n0 = {
            0., 0., 0.
        };
        Dexp = {
            0., 0., 0.
        };
        Dn0 = -arma :: zeros< arma :: mat >(3, 3);
    } else {
        n0 = a0 / d0;
        Dexp = -std :: exp( std :: pow(d0 / c, alpha) ) * alpha / c *pow(d0 / c, alpha - 1) * n0;
        Dn0 = -arma :: eye< arma :: mat >(3, 3) / d0 + a0 *a0.t() / pow(d0, 3);
    }
    TangentSelf = n0 * Dexp.t() + Dn0 * ( std :: exp( std :: pow(d0 / c, alpha) ) - 1 );

    // Linear part
    TangentSelf = TangentSelf - arma :: eye< arma :: mat >(3, 3) * ConnectionCoords.size();

    // Assemble self part to tangent
    Tangent( arma :: span(0, 2), arma :: span(0, 2) ) = TangentSelf;

    // delta x_j part
    for ( unsigned int i = 0; i < ConnectionCoords.size(); i++ ) {
        for ( int j = 0; j < 3; j++ ) {
            Tangent(j, 3 * ( i + 1 ) + j) = 1.0; // /ConnectionCoords.size();
        }
    }
    return Tangent;
}

void SpringSmoother :: Smooth(std :: vector< VertexType * >Vertices, MeshData *Mesh)
{
    double MAXCHANGE = 1e-4 * charlength;

    std::vector<std::vector<VertexType *>> Connections = this->GetConnectivityVector(Vertices);

    // Create vectors for current and previous positions for all involved vertices (even those vertices connected to a vertex in Vertices vector)
    std :: map<VertexType *, arma::vec3 >OriginalPositions;
    std :: map<VertexType *, arma::vec3 >CurrentPositions;
    std :: map<VertexType *, arma::vec3 >PreviousPositions;

    for (std::vector<VertexType *> VertexList: Connections) {
        for (VertexType *v: VertexList) {
            OriginalPositions[v] = v->get_c_vec();
            CurrentPositions[v] = v->get_c_vec();
            PreviousPositions[v] = v->get_c_vec();
        }
    }

    for (VertexType *v: Vertices) {
        OriginalPositions[v] = v->get_c_vec();
        CurrentPositions[v] = v->get_c_vec();
        PreviousPositions[v] = v->get_c_vec();
    }

    this->CheckPenetration(&Vertices, (MeshManipulations*) Mesh);

    double deltamax = 1e8;

    while (deltamax > MAXCHANGE) {

        deltamax = 0.0;
        size_t iter=0;

        for (VertexType *v: Vertices) {

            std::vector<arma::vec3> ConnectionCoords;
            for (VertexType *cv: Connections[iter]) {
                ConnectionCoords.push_back(cv->get_c_vec());
            }

            arma::vec3 xc = CurrentPositions[v];
            arma::vec3 x0 = OriginalPositions[v];
            arma::vec3 R = ComputeOutOfBalance(ConnectionCoords, xc, x0, alpha, c);

            double err = arma::norm(R);

            while (err>1e-5) {
                arma::mat K = ComputeAnalyticalTangent(ConnectionCoords, xc, x0, alpha, c);
                arma::vec d = -arma::solve(K,R);
                xc = xc + d;
                R = ComputeOutOfBalance(ConnectionCoords, xc, x0, alpha, c);
                err = arma::norm(R);
            }


            for (int i=0; i<3; i++) {
                if (!v->Fixed[i]) {
                    CurrentPositions[v][i] = xc[i];
                    v->set_c(CurrentPositions[v][i], i);
                }
            }

            double delta = arma::norm(CurrentPositions[v]-PreviousPositions[v]);

            for (int i=0; i<3; i++) {
                if (!v->Fixed[i]) {
                    PreviousPositions[v][i] = CurrentPositions[v][i];
                }
            }

            deltamax = std::max(delta, deltamax);
            STATUS("%c[2K\r\tIteration %u end with deltamax=%f\r", 27, iter, deltamax);
            fflush(stdout);
            iter++;
        }
    }
    STATUS("\n", 0);

}


std::ostream &operator<<(std::ostream &stream, const SpringSmoother &Smoother)
{
    stream << "\talpha = " << Smoother.alpha << ", ";
    stream << "c = " << Smoother.c << ", c_factor = " << Smoother.c_factor << "\n";
    return stream;
}


}
