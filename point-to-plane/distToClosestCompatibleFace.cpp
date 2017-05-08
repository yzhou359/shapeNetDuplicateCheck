#include <mex.h>
#include <math.h>
#include "Matrix2D.h"
#include "dist2Triangle.h"
#include <iostream>

#define INF (mxGetInf())

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // Make sure we have the correct # of input args
    if (nrhs != 7)
        mexErrMsgTxt("Incorrect # of input arguments.\n");
    
    // Get the size of the input
    int numpoints0 = mxGetN(prhs[0]);
    int numfaces  = mxGetN(prhs[1]);
	int numpoints = mxGetN(prhs[4]);
    
    // inputs: 
	// Ps (vertex positions 3 x numvertices)
	// F (face indices 3 x numfaces)
	// Nf (face normals 3 x numfaces)
	// Nv (vertex normals 3 x points)
	// Qs (query point positions, 3 x numquerypoints)
	// Nq (normals of query points, 3 x numquerypoints)
    Matrix2D<double> Ps((double*) mxGetPr(prhs[0]), 3, numpoints0);
    Matrix2D<double> F ((double*) mxGetPr(prhs[1]), 3, numfaces);
    Matrix2D<double> Nf((double*) mxGetPr(prhs[2]), 3, numfaces);
    Matrix2D<double> Nv((double*) mxGetPr(prhs[3]), 3, numpoints0);
    Matrix2D<double> Qs((double*) mxGetPr(prhs[4]), 3, numpoints);
    Matrix2D<double> Nq((double*) mxGetPr(prhs[5]), 3, numpoints);
    Matrix2D<double> NormalCheck((double*) mxGetPr(prhs[6]), 1, 1);
    
    // outputs
    plhs[0] = mxCreateDoubleMatrix(1, numpoints, mxREAL);
    double* Ds =  mxGetPr(plhs[0]);
    

    // initialize output dists to INF
    int q;
    for (q = 0; q < numpoints; q++) {
        Ds[q] = INF;
    }
    
    int f, i, j, k, c;
    double d, d_try, p[3], p0[3], p1[3], p2[3], Np[3], r;
    retValues ret;
    
    for (q = 0; q < numpoints; q++) 
	{
        d = INF;
        
        p[0] = Qs(0,q);
        p[1] = Qs(1,q);
        p[2] = Qs(2,q);

		//int pidx[3] = { 0,0,0 };
        
        for (f = 0; f < numfaces; f++) 
		{
           // if (Nf(0,f)*Nq(0,q) + Nf(1,f)*Nq(1,q) + Nf(2,f)*Nq(2,q) >= 0) 
			if ( true )
			{
                
                i = ((int) F(0,f)) ;
                j = ((int) F(1,f)) ;
                k = ((int) F(2,f)) ;
                
                //mexPrintf("%d, %f\n", i, F(1,f));
                
                p0[0] = Ps(0,i);
                p0[1] = Ps(1,i);
                p0[2] = Ps(2,i);
                
                p1[0] = Ps(0,j);
                p1[1] = Ps(1,j);
                p1[2] = Ps(2,j);
                
                p2[0] = Ps(0,k);
                p2[1] = Ps(1,k);
                p2[2] = Ps(2,k);

				
                
                ret = dist2Triangle(p0, p1, p2, p);
                d_try = ret.distance;
                
                //mexPrintf("%f %f\n", ret.s, ret.t);
                
                if (d_try < 0)
                    d_try = 0;
                
                if (d_try < d) {
                    
                    r = 1.0 - ret.s - ret.t;
                    
                    for (c = 0; c < 3; c++) {
                        Np[c] = r*Nv(c, i) + ret.s*Nv(c, j) + ret.t*Nv(c, k);
                    }
                    vec3normalize(Np);
                                        
					if (Np[0] * Nq(0, q) + Np[1] * Nq(1, q) + Np[2] * Nq(2, q) >= cos(NormalCheck(0, 0)))
					{
						d = d_try;
						//pidx[0] = i;
						//pidx[1] = j;
						//pidx[2] = k;
					}
                        
                }
            }
        }

        Ds[q] = d;

		//if (q == 963)
		//	std::cout << pidx[0] << ' ' << pidx[1] << ' ' << pidx[2] << std::endl;
    }

}

void mexFunction_simple(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	// Make sure we have the correct # of input args
	if (nrhs != 7)
		mexErrMsgTxt("Incorrect # of input arguments.\n");

	// Get the size of the input
	int numpoints0 = mxGetN(prhs[0]);
	int numfaces = mxGetN(prhs[1]);
	int numpoints = mxGetN(prhs[2]);

	// inputs: 
	// Ps (vertex positions 3 x numvertices)
	// F (face indices 3 x numfaces)
	// Qs (query point positions, 3 x numquerypoints)
	Matrix2D<double> Ps((double*)mxGetPr(prhs[0]), 3, numpoints0);
	Matrix2D<double> F((double*)mxGetPr(prhs[1]), 3, numfaces);
	Matrix2D<double> Qs((double*)mxGetPr(prhs[2]), 3, numpoints);

	// outputs
	plhs[0] = mxCreateDoubleMatrix(1, numpoints, mxREAL);
	double* Ds = mxGetPr(plhs[0]);


	// initialize output dists to INF
	int q;
	for (q = 0; q < numpoints; q++) {
		Ds[q] = INF;
	}

	int f, i, j, k, c;
	double d, d_try, p[3], p0[3], p1[3], p2[3], Np[3], r;
	retValues ret;

	for (q = 0; q < numpoints; q++)
	{
		d = INF;

		p[0] = Qs(0, q);
		p[1] = Qs(1, q);
		p[2] = Qs(2, q);

		for (f = 0; f < numfaces; f++)
		{

			i = ((int)F(0, f));
			j = ((int)F(1, f));
			k = ((int)F(2, f));

			//mexPrintf("%d, %f\n", i, F(1,f));

			p0[0] = Ps(0, i);
			p0[1] = Ps(1, i);
			p0[2] = Ps(2, i);

			p1[0] = Ps(0, j);
			p1[1] = Ps(1, j);
			p1[2] = Ps(2, j);

			p2[0] = Ps(0, k);
			p2[1] = Ps(1, k);
			p2[2] = Ps(2, k);

			ret = dist2Triangle(p0, p1, p2, p);
			d_try = ret.distance;

			//mexPrintf("%f %f\n", ret.s, ret.t);

			if (d_try < 0)
				d_try = 0;

			if (d_try < d)
				d = d_try;

		}

		Ds[q] = d;
	}

}