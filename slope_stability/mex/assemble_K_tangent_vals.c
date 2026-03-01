/*
 * assemble_K_tangent_vals.c
 *
 * Element-level tangent stiffness assembly for 3D P2 FEM.
 *
 * Computes V_tang = values of K_tangent(Q,Q) at the precomputed sparsity
 * pattern positions, using element-local quadrature:
 *
 *   K_tangent = sum_e sum_q  B_eq' * (w_q * D_eq) * B_eq
 *
 * where B_eq is reconstructed from DPhi1/2/3 at each integration point,
 * D_eq = reshape(DS(:,g), 6, 6), and w_q = WEIGHT(g).
 *
 * The scatter_map (n_local_dof^2 x n_e, int64) maps each local matrix
 * entry to a 1-based index in V_tang.  Entries == 0 are Dirichlet and skipped.
 *
 * Usage (from Octave/MATLAB):
 *   V_tang = assemble_K_tangent_vals(DPhi1, DPhi2, DPhi3, DS, WEIGHT, ...
 *                                     scatter_map, n_q, nnz_out);
 *
 * Compile:
 *   mkoctfile --mex -O2 -fopenmp assemble_K_tangent_vals.c -lgomp
 *   (or)  mex -O CFLAGS='$CFLAGS -fopenmp' -lgomp assemble_K_tangent_vals.c
 */

#include "mex.h"
#include <string.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Fixed dimensions for 3D: n_strain=6, dim=3 */
#define NS 6
#define DIM 3

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* ---- Input validation ---- */
    if (nrhs != 8)
        mexErrMsgIdAndTxt("assemble_K_tangent_vals:nrhs",
                          "Eight inputs required: DPhi1, DPhi2, DPhi3, DS, WEIGHT, scatter_map, n_q, nnz_out");
    if (nlhs > 1)
        mexErrMsgIdAndTxt("assemble_K_tangent_vals:nlhs",
                          "At most one output.");

    /* ---- Unpack inputs ---- */
    const double *DPhi1 = mxGetPr(prhs[0]);  /* n_p x n_int */
    const double *DPhi2 = mxGetPr(prhs[1]);  /* n_p x n_int */
    const double *DPhi3 = mxGetPr(prhs[2]);  /* n_p x n_int */
    const double *DS = mxGetPr(prhs[3]);     /* 36 x n_int  */
    const double *WEIGHT = mxGetPr(prhs[4]); /* 1 x n_int   */

    /* scatter_map is int64 (n_local_dof^2 x n_e) */
    if (!mxIsInt64(prhs[5]))
        mexErrMsgIdAndTxt("assemble_K_tangent_vals:type",
                          "scatter_map must be int64.");
    const int64_t *scatter_map = (const int64_t *)mxGetData(prhs[5]);

    const int n_q = (int)mxGetScalar(prhs[6]);
    const int nnz_out = (int)mxGetScalar(prhs[7]);

    const mwSize n_p = mxGetM(prhs[0]);   /* nodes per element   */
    const mwSize n_int = mxGetN(prhs[0]); /* total int. points   */
    const mwSize n_e = mxGetN(prhs[5]);   /* number of elements  */

    const int n_local_dof = (int)(DIM * n_p);    /* 30 for P2 tet       */
    const int n_ld2 = n_local_dof * n_local_dof; /* 900 */

    /* ---- Allocate output (zeros) ---- */
    plhs[0] = mxCreateDoubleMatrix((mwSize)nnz_out, 1, mxREAL);
    double *V_tang = mxGetPr(plhs[0]);

/* ---- Element loop (OpenMP parallel) ---- */
#pragma omp parallel
    {
        /* Thread-private local arrays on the stack */
        double Ke[900];      /* 30 x 30, column-major */
        double B_eq[6 * 30]; /* NS x n_local_dof, column-major */
        double D_eq[36];     /* NS x NS, column-major */
        double tmp[6 * 30];  /* NS x n_local_dof temp */

#pragma omp for schedule(dynamic, 256)
        for (mwSize e = 0; e < n_e; e++)
        {
            /* Zero local stiffness matrix */
            memset(Ke, 0, sizeof(double) * (size_t)n_ld2);

            const mwSize g_base = e * (mwSize)n_q;

            /* Loop over quadrature points in this element */
            for (int q = 0; q < n_q; q++)
            {
                const mwSize g = g_base + (mwSize)q; /* global int. pt. */
                const double w = WEIGHT[g];

                /* Build D_eq = w * reshape(DS(:,g), 6, 6) */
                const double *ds_col = DS + g * NS * NS;
                for (int k = 0; k < NS * NS; k++)
                    D_eq[k] = w * ds_col[k];

                /* Build B_eq (NS x n_local_dof) from DPhi */
                memset(B_eq, 0, sizeof(double) * NS * (size_t)n_local_dof);
                const double *dp1 = DPhi1 + g * n_p;
                const double *dp2 = DPhi2 + g * n_p;
                const double *dp3 = DPhi3 + g * n_p;

                for (mwSize i = 0; i < n_p; i++)
                {
                    const double dN1 = dp1[i];
                    const double dN2 = dp2[i];
                    const double dN3 = dp3[i];
                    const int c = (int)i * DIM; /* column base */

                    /* eps_11: row 0, col c+0 */
                    B_eq[0 + (c + 0) * NS] = dN1;
                    /* eps_22: row 1, col c+1 */
                    B_eq[1 + (c + 1) * NS] = dN2;
                    /* eps_33: row 2, col c+2 */
                    B_eq[2 + (c + 2) * NS] = dN3;
                    /* gamma_12: row 3, col c+0 and c+1 */
                    B_eq[3 + (c + 0) * NS] = dN2;
                    B_eq[3 + (c + 1) * NS] = dN1;
                    /* gamma_23: row 4, col c+1 and c+2 */
                    B_eq[4 + (c + 1) * NS] = dN3;
                    B_eq[4 + (c + 2) * NS] = dN2;
                    /* gamma_13: row 5, col c+0 and c+2 */
                    B_eq[5 + (c + 0) * NS] = dN3;
                    B_eq[5 + (c + 2) * NS] = dN1;
                }

                /* tmp = D_eq * B_eq   (NS x n_local_dof) */
                for (int j = 0; j < n_local_dof; j++)
                {
                    for (int ii = 0; ii < NS; ii++)
                    {
                        double s = 0.0;
                        for (int kk = 0; kk < NS; kk++)
                            s += D_eq[ii + kk * NS] * B_eq[kk + j * NS];
                        tmp[ii + j * NS] = s;
                    }
                }

                /* Ke += B_eq' * tmp   (n_local_dof x n_local_dof) */
                for (int j = 0; j < n_local_dof; j++)
                {
                    for (int ii = 0; ii < n_local_dof; ii++)
                    {
                        double s = 0.0;
                        for (int kk = 0; kk < NS; kk++)
                            s += B_eq[kk + ii * NS] * tmp[kk + j * NS];
                        Ke[ii + j * n_local_dof] += s;
                    }
                }
            } /* end quadrature point loop */

            /* Scatter Ke into V_tang using atomic adds */
            const int64_t *smap = scatter_map + e * (mwSize)n_ld2;
            for (int k = 0; k < n_ld2; k++)
            {
                const int64_t idx = smap[k];
                if (idx > 0)
                {
#pragma omp atomic update
                    V_tang[idx - 1] += Ke[k];
                }
            }
        } /* end element loop */
    } /* end parallel region */
}
