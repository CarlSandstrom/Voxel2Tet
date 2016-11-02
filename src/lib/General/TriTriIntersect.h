#ifndef TRITRIINTERSECT_H
#define TRITRIINTERSECT_H

int coplanar_tri_tri(double N[3], double V0[3], double V1[3], double V2[3], double U0[3], double U1[3], double U2[3]);

int tri_tri_intersect(double V0[3], double V1[3], double V2[3], double U0[3], double U1[3], double U2[3]);

bool point_in_tri(double V0[3], double V1[3], double V2[3], double P[3]);

bool tri_tri_intersect_shared_edge(double s0[3], double s1[3], double u0[3], double u1[3]);

// https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
int triangle_ray_intersection(double V1[3],  // Triangle vertices
                              double V2[3],
                              double V3[3],
                              double O[3], //Ray origin
                              double D[3], //Ray direction
                              float *out);

#endif // TRITRIINTERSECT_H
