/*
 * constitutive_problem_3D_S_mex.c
 *
 * Mex wrapper for stress-only evaluation of the 3D Mohr-Coulomb
 * constitutive model.  OpenMP-parallel over integration points.
 *
 * Usage from Octave / MATLAB:
 *   S = constitutive_problem_3D_S_mex(E, c_bar, sin_phi, shear, bulk, lame)
 *
 *   E        – 6 × n_int   (engineering strain)
 *   c_bar    – 1 × n_int   (reduced cohesion)
 *   sin_phi  – 1 × n_int   (sine of reduced friction angle)
 *   shear    – 1 × n_int   (shear modulus)
 *   bulk     – 1 × n_int   (bulk modulus)
 *   lame     – 1 × n_int   (first Lamé parameter)
 *
 *   S        – 6 × n_int   (stress)
 */

#include "mex.h"
#include "constitutive_3D_kernel.h"

#ifdef _OPENMP
#include <omp.h>
#endif

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 6)
        mexErrMsgTxt("constitutive_problem_3D_S_mex: "
                     "6 inputs required (E, c_bar, sin_phi, shear, bulk, lame)");
    if (nlhs > 1)
        mexErrMsgTxt("constitutive_problem_3D_S_mex: 1 output (S)");

    /* ---- Inputs ---- */
    const double *E = mxGetPr(prhs[0]);
    const double *c_bar = mxGetPr(prhs[1]);
    const double *sin_phi = mxGetPr(prhs[2]);
    const double *shear = mxGetPr(prhs[3]);
    const double *bulk = mxGetPr(prhs[4]);
    const double *lam = mxGetPr(prhs[5]);

    mwSize n_int = mxGetN(prhs[0]); /* number of integration points */

    /* ---- Output ---- */
    plhs[0] = mxCreateDoubleMatrix(6, n_int, mxREAL);
    double *S = mxGetPr(plhs[0]);

/* ---- Parallel loop ---- */
#pragma omp parallel for schedule(static)
    for (mwSize p = 0; p < n_int; p++)
    {
        constitutive_3D_point(E + 6 * p,
                              c_bar[p], sin_phi[p],
                              shear[p], bulk[p], lam[p],
                              S + 6 * p, NULL);
    }
}
