#include <armadillo>
#include <stdio.h>

typedef double vertex [ 3 ];

int testfunction(vertex v)
{
    for ( int i = 0; i < 3; i++ ) {
        printf("%f\n", v [ i ]);
    }
    return 0;
}

int main()
{
    vertex v = {
        1.0, 2.0, 1.0
    };
    testfunction(v);
    arma :: mat A(5, 5);
    arma :: mat B = arma :: randu< arma :: mat >(5, 5);
    A.eye();
    arma :: mat C = A * B.t();
    C.print();
}
