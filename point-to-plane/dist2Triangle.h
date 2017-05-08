#include <mex.h>
#include <math.h>
#include "Matrix2D.h"
#include <iostream>
/* --------------------------------------------------------------------- */
inline void vec3sub(const double* v1, const double* v2, double* sub) {
    for (int i = 0; i < 3; i++) {
        sub[i] = v1[i] - v2[i];
    }
}

/* --------------------------------------------------------------------- */
inline void vec3add(const double* v1, const double* v2, double* sum) {
    for (int i = 0; i < 3; i++) {
        sum[i] = v1[i] + v2[i];
    }
}

/* --------------------------------------------------------------------- */
inline double vec3dot(const double* v1, const double* v2) {
    double dot = 0;
    for (int i = 0; i < 3; i++) {
        dot += v1[i]*v2[i];
    }
    return dot;
}

/* --------------------------------------------------------------------- */
inline double* vec3cross(const double* v1, const double* v2) {
	double cross[3];
	cross[0] = v1[1] * v2[2] - v1[2] * v2[1];
	cross[1] = v1[2] * v2[0] - v1[0] * v2[2];
	cross[2] = v1[0] * v2[1] - v1[1] * v2[0];
	return cross;
}

/* --------------------------------------------------------------------- */
inline void vec3copy(const double* v1, double* v2) {
    for (int i = 0; i < 3; i++) {
        v2[i] = v1[i];
    }
}

/* --------------------------------------------------------------------- */
inline double vec3norm(const double* v) {
    return sqrt(vec3dot(v,v));
}

/* --------------------------------------------------------------------- */
inline void vec3normalize(double* v) {
    double n = vec3norm(v);
    for (int i = 0; i < 3; i++)
        v[i] = v[i]/n;
}

/* --------------------------------------------------------------------- */

struct retValues {
	double distance;
	double s, t;
};

retValues dist2Triangle(const double* p0,
	const double* p1,
	const double* p2,
	const double* p);


void mexFunction_simple(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);