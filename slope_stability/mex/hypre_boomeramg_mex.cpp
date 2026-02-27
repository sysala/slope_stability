#include "mex.h"
#include "matrix.h"

#include "HYPRE.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_config.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_parcsr_mv.h"
#include "HYPRE_utilities.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

struct Instance
{
    mwSize n = 0;
    mwSize nnz = 0;
    int num_functions = 1;
    int num_interp_vectors = 0;

    HYPRE_IJMatrix A_ij = nullptr;
    HYPRE_ParCSRMatrix A = nullptr;
    HYPRE_IJVector b_ij = nullptr;
    HYPRE_IJVector x_ij = nullptr;
    HYPRE_ParVector b = nullptr;
    HYPRE_ParVector x = nullptr;
    HYPRE_Solver solver = nullptr;

    std::vector<HYPRE_BigInt> rows;
    HYPRE_IJVector *interp_ij = nullptr;
    HYPRE_ParVector *interp_par = nullptr;
};

struct BoomerOptions
{
    int threads = 16;
    int print_level = 0;
    int max_levels = 25;

    int coarsen_type = 8;
    int agg_num_levels = 0;
    int interp_type = 17;
    int p_max_elmts = 4;
    double strong_threshold = 0.5;
    bool has_strong_threshold_r = false;
    double strong_threshold_r = 0.0;

    int relax_type = 8;
    int relax_sweeps = 1;
    int relax_coarse_type = 8;
    int cycle_type = 1;

    int nodal = 4;
    int nodal_diag = 1;
    int nodal_levels = 0;
    int keep_same_sign = 0;

    int interp_vec_variant = 2;
    int interp_vec_qmax = 4;
    bool has_interp_vec_abs_qtrunc = false;
    double interp_vec_abs_qtrunc = 0.0;
    int smooth_interp_vectors = 1;
    int interp_refine = 0;

    bool has_trunc_factor = false;
    double trunc_factor = 0.0;
    bool has_agg_interp_type = false;
    int agg_interp_type = 4;
    bool has_agg_trunc_factor = false;
    double agg_trunc_factor = 0.0;
    bool has_agg_p_max_elmts = false;
    int agg_p_max_elmts = 0;
    bool has_max_coarse_size = false;
    int max_coarse_size = 9;
    bool has_min_coarse_size = false;
    int min_coarse_size = 1;

    bool use_as_preconditioner = true; /* one application */
    int max_iter = 1;
    double tol = 0.0;

    bool has_num_functions = false;
    int num_functions = 1;
};

std::unordered_map<std::string, std::unique_ptr<Instance>> g_instances;
bool g_hypre_initialized = false;
bool g_at_exit_registered = false;

void destroy_instance(Instance &inst)
{
    if (inst.solver)
    {
        if (inst.num_interp_vectors > 0)
        {
            /* Avoid HYPRE teardown paths that may touch user-owned interp vectors. */
            HYPRE_BoomerAMGSetInterpVectors(inst.solver, 0, nullptr);
        }
        HYPRE_BoomerAMGDestroy(inst.solver);
        inst.solver = nullptr;
    }

    if (inst.interp_ij)
    {
        for (int k = 0; k < inst.num_interp_vectors; k++)
        {
            if (inst.interp_ij[k]) { HYPRE_IJVectorDestroy(inst.interp_ij[k]); }
        }
        std::free(inst.interp_ij);
        inst.interp_ij = nullptr;
    }
    if (inst.interp_par)
    {
        std::free(inst.interp_par);
        inst.interp_par = nullptr;
    }

    if (inst.x_ij)
    {
        HYPRE_IJVectorDestroy(inst.x_ij);
        inst.x_ij = nullptr;
        inst.x = nullptr;
    }
    if (inst.b_ij)
    {
        HYPRE_IJVectorDestroy(inst.b_ij);
        inst.b_ij = nullptr;
        inst.b = nullptr;
    }
    if (inst.A_ij)
    {
        HYPRE_IJMatrixDestroy(inst.A_ij);
        inst.A_ij = nullptr;
        inst.A = nullptr;
    }
}

void cleanup_all()
{
    for (auto &kv : g_instances)
    {
        destroy_instance(*kv.second);
    }
    g_instances.clear();

    if (g_hypre_initialized)
    {
        HYPRE_Finalize();
        g_hypre_initialized = false;
    }
}

void at_exit()
{
    cleanup_all();
}

[[noreturn]] void mex_error(const char *id, const std::string &msg)
{
    mexErrMsgIdAndTxt(id, "%s", msg.c_str());
    std::abort();
}

std::string to_id_string(const mxArray *id_arr)
{
    if (mxIsChar(id_arr))
    {
        char *s = mxArrayToString(id_arr);
        if (!s) { mex_error("hypre:badId", "Failed to convert char instance_id."); }
        std::string out(s);
        mxFree(s);
        return out;
    }
    if (mxIsNumeric(id_arr) && mxGetNumberOfElements(id_arr) == 1)
    {
        double v = mxGetScalar(id_arr);
        long long iv = static_cast<long long>(std::llround(v));
        return std::to_string(iv);
    }
    mex_error("hypre:badId", "instance_id must be char or numeric scalar.");
    return std::string();
}

const mxArray *get_field_or_null(const mxArray *opts, const char *name)
{
    if (!opts || mxIsEmpty(opts)) { return nullptr; }
    if (!mxIsStruct(opts)) { mex_error("hypre:badOpts", "opts must be a struct or empty."); }
    return mxGetField(opts, 0, name);
}

bool get_int_field(const mxArray *opts, const char *name, int &dst)
{
    const mxArray *f = get_field_or_null(opts, name);
    if (!f) { return false; }
    if (!(mxIsNumeric(f) || mxIsLogical(f)) || mxGetNumberOfElements(f) != 1)
    {
        mex_error("hypre:badOpts", std::string("opts.") + name + " must be a numeric scalar.");
    }
    dst = static_cast<int>(std::llround(mxGetScalar(f)));
    return true;
}

bool get_double_field(const mxArray *opts, const char *name, double &dst)
{
    const mxArray *f = get_field_or_null(opts, name);
    if (!f) { return false; }
    if (!(mxIsNumeric(f) || mxIsLogical(f)) || mxGetNumberOfElements(f) != 1)
    {
        mex_error("hypre:badOpts", std::string("opts.") + name + " must be a numeric scalar.");
    }
    dst = mxGetScalar(f);
    return true;
}

bool get_bool_field(const mxArray *opts, const char *name, bool &dst)
{
    const mxArray *f = get_field_or_null(opts, name);
    if (!f) { return false; }
    if (mxIsLogicalScalar(f))
    {
        dst = mxIsLogicalScalarTrue(f);
        return true;
    }
    if (mxIsNumeric(f) && mxGetNumberOfElements(f) == 1)
    {
        dst = (mxGetScalar(f) != 0.0);
        return true;
    }
    mex_error("hypre:badOpts", std::string("opts.") + name + " must be logical/numeric scalar.");
    return false;
}

BoomerOptions parse_options(const mxArray *opts, mwSize n, mwSize null_k)
{
    BoomerOptions o;

    get_int_field(opts, "threads", o.threads);
    get_int_field(opts, "print_level", o.print_level);
    get_int_field(opts, "max_levels", o.max_levels);

    get_int_field(opts, "coarsen_type", o.coarsen_type);
    get_int_field(opts, "agg_num_levels", o.agg_num_levels);
    get_int_field(opts, "interp_type", o.interp_type);
    get_int_field(opts, "p_max_elmts", o.p_max_elmts);
    get_double_field(opts, "strong_threshold", o.strong_threshold);
    if (get_double_field(opts, "strong_threshold_r", o.strong_threshold_r)) { o.has_strong_threshold_r = true; }

    get_int_field(opts, "relax_type", o.relax_type);
    get_int_field(opts, "relax_sweeps", o.relax_sweeps);
    get_int_field(opts, "relax_coarse_type", o.relax_coarse_type);
    get_int_field(opts, "cycle_type", o.cycle_type);

    get_int_field(opts, "nodal", o.nodal);
    get_int_field(opts, "nodal_diag", o.nodal_diag);
    get_int_field(opts, "nodal_levels", o.nodal_levels);
    get_int_field(opts, "keep_same_sign", o.keep_same_sign);

    get_int_field(opts, "interp_vec_variant", o.interp_vec_variant);
    get_int_field(opts, "interp_vec_qmax", o.interp_vec_qmax);
    if (get_double_field(opts, "interp_vec_abs_qtrunc", o.interp_vec_abs_qtrunc)) { o.has_interp_vec_abs_qtrunc = true; }
    get_int_field(opts, "smooth_interp_vectors", o.smooth_interp_vectors);
    get_int_field(opts, "interp_refine", o.interp_refine);

    if (get_double_field(opts, "trunc_factor", o.trunc_factor)) { o.has_trunc_factor = true; }
    if (get_int_field(opts, "agg_interp_type", o.agg_interp_type)) { o.has_agg_interp_type = true; }
    if (get_double_field(opts, "agg_trunc_factor", o.agg_trunc_factor)) { o.has_agg_trunc_factor = true; }
    if (get_int_field(opts, "agg_p_max_elmts", o.agg_p_max_elmts)) { o.has_agg_p_max_elmts = true; }
    if (get_int_field(opts, "max_coarse_size", o.max_coarse_size)) { o.has_max_coarse_size = true; }
    if (get_int_field(opts, "min_coarse_size", o.min_coarse_size)) { o.has_min_coarse_size = true; }

    get_bool_field(opts, "use_as_preconditioner", o.use_as_preconditioner);
    get_int_field(opts, "max_iter", o.max_iter);
    get_double_field(opts, "tol", o.tol);

    if (get_int_field(opts, "num_functions", o.num_functions)) { o.has_num_functions = true; }

    if (!o.has_num_functions)
    {
        if (null_k > 0 && (n % 3 == 0))
        {
            o.num_functions = 3;
        }
        else
        {
            o.num_functions = 1;
        }
    }

    if (o.use_as_preconditioner)
    {
        o.max_iter = 1;
        o.tol = 0.0;
    }

    return o;
}

std::vector<HYPRE_Int> build_dof_func(const mxArray *opts, mwSize n, int num_functions)
{
    std::vector<HYPRE_Int> mapping(n);

    const mxArray *dof_field = get_field_or_null(opts, "dof_func");
    if (dof_field && !mxIsEmpty(dof_field))
    {
        if (!mxIsNumeric(dof_field))
        {
            mex_error("hypre:badOpts", "opts.dof_func must be numeric.");
        }
        if (mxGetNumberOfElements(dof_field) != n)
        {
            mex_error("hypre:badOpts", "opts.dof_func must have length n.");
        }
        const double *p = mxGetPr(dof_field);
        for (mwSize i = 0; i < n; i++)
        {
            int v = static_cast<int>(std::llround(p[i]));
            if (v < 0 || v >= num_functions)
            {
                mex_error("hypre:badOpts", "opts.dof_func values must be in [0, num_functions-1].");
            }
            mapping[i] = static_cast<HYPRE_Int>(v);
        }
        return mapping;
    }

    if (num_functions <= 1)
    {
        std::fill(mapping.begin(), mapping.end(), static_cast<HYPRE_Int>(0));
        return mapping;
    }

    for (mwSize i = 0; i < n; i++)
    {
        mapping[i] = static_cast<HYPRE_Int>(i % static_cast<mwSize>(num_functions));
    }
    return mapping;
}

void build_rowwise_from_sparse(const mxArray *A,
                               std::vector<HYPRE_Int> &row_nnz,
                               std::vector<HYPRE_Int> &row_idx,
                               std::vector<HYPRE_BigInt> &rows,
                               std::vector<HYPRE_BigInt> &cols,
                               std::vector<HYPRE_Complex> &vals)
{
    if (!mxIsSparse(A) || mxIsComplex(A) || !mxIsDouble(A))
    {
        mex_error("hypre:badA", "A must be a real double sparse matrix.");
    }

    mwSize m = mxGetM(A);
    mwSize n = mxGetN(A);
    if (m != n)
    {
        mex_error("hypre:badA", "A must be square.");
    }

    const mwIndex *jc = mxGetJc(A);
    const mwIndex *ir = mxGetIr(A);
    const double *pr = mxGetPr(A);
    mwSize nnz = static_cast<mwSize>(jc[n]);

    std::vector<HYPRE_Int> row_ptr(n + 1, 0);
    for (mwSize j = 0; j < n; j++)
    {
        for (mwIndex p = jc[j]; p < jc[j + 1]; p++)
        {
            row_ptr[ir[p] + 1]++;
        }
    }
    for (mwSize i = 0; i < n; i++)
    {
        row_ptr[i + 1] += row_ptr[i];
    }

    row_nnz.resize(n);
    row_idx.resize(n);
    rows.resize(n);
    for (mwSize i = 0; i < n; i++)
    {
        row_nnz[i] = row_ptr[i + 1] - row_ptr[i];
        row_idx[i] = row_ptr[i];
        rows[i] = static_cast<HYPRE_BigInt>(i);
    }

    cols.resize(nnz);
    vals.resize(nnz);
    std::vector<HYPRE_Int> next(row_ptr.begin(), row_ptr.end());

    for (mwSize j = 0; j < n; j++)
    {
        for (mwIndex p = jc[j]; p < jc[j + 1]; p++)
        {
            mwSize i = ir[p];
            HYPRE_Int dst = next[i]++;
            cols[dst] = static_cast<HYPRE_BigInt>(j);
            vals[dst] = static_cast<HYPRE_Complex>(pr[p]);
        }
    }
}

double setup_instance(Instance &inst,
                      const mxArray *A,
                      const mxArray *Z,
                      const mxArray *opts)
{
    const mwSize n = mxGetN(A);
    inst.n = n;
    inst.num_interp_vectors = (Z && !mxIsEmpty(Z)) ? static_cast<int>(mxGetN(Z)) : 0;
    if (Z && !mxIsEmpty(Z))
    {
        if (!mxIsDouble(Z) || mxIsComplex(Z))
        {
            mex_error("hypre:badZ", "null space matrix must be real double.");
        }
        if (mxGetM(Z) != n)
        {
            mex_error("hypre:badZ", "null space must have size n x k.");
        }
    }

    if (!g_hypre_initialized)
    {
        HYPRE_Init();
        g_hypre_initialized = true;
    }

    BoomerOptions o = parse_options(opts, n, static_cast<mwSize>(inst.num_interp_vectors));
    inst.num_functions = o.num_functions;

#ifdef _OPENMP
    if (o.threads > 0)
    {
        omp_set_num_threads(o.threads);
    }
#endif

    std::vector<HYPRE_Int> row_nnz;
    std::vector<HYPRE_Int> row_idx;
    std::vector<HYPRE_BigInt> rows;
    std::vector<HYPRE_BigInt> cols;
    std::vector<HYPRE_Complex> vals;
    build_rowwise_from_sparse(A, row_nnz, row_idx, rows, cols, vals);
    inst.nnz = vals.size();
    inst.rows = rows;

    MPI_Comm comm = 0; /* single process, MPI disabled build */
    HYPRE_IJMatrixCreate(comm, 0, static_cast<HYPRE_BigInt>(n - 1), 0, static_cast<HYPRE_BigInt>(n - 1), &inst.A_ij);
    HYPRE_IJMatrixSetObjectType(inst.A_ij, HYPRE_PARCSR);
    HYPRE_IJMatrixSetRowSizes(inst.A_ij, row_nnz.data());
    HYPRE_IJMatrixInitialize(inst.A_ij);
    HYPRE_IJMatrixSetValues2(inst.A_ij, static_cast<HYPRE_Int>(n), row_nnz.data(), rows.data(),
                             row_idx.data(), cols.data(), vals.data());
    HYPRE_IJMatrixAssemble(inst.A_ij);
    HYPRE_IJMatrixGetObject(inst.A_ij, reinterpret_cast<void **>(&inst.A));

    std::vector<HYPRE_Complex> zeros(n, 0.0);

    HYPRE_IJVectorCreate(comm, 0, static_cast<HYPRE_BigInt>(n - 1), &inst.b_ij);
    HYPRE_IJVectorSetObjectType(inst.b_ij, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(inst.b_ij);
    HYPRE_IJVectorSetValues(inst.b_ij, static_cast<HYPRE_Int>(n), rows.data(), zeros.data());
    HYPRE_IJVectorAssemble(inst.b_ij);
    HYPRE_IJVectorGetObject(inst.b_ij, reinterpret_cast<void **>(&inst.b));

    HYPRE_IJVectorCreate(comm, 0, static_cast<HYPRE_BigInt>(n - 1), &inst.x_ij);
    HYPRE_IJVectorSetObjectType(inst.x_ij, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(inst.x_ij);
    HYPRE_IJVectorSetValues(inst.x_ij, static_cast<HYPRE_Int>(n), rows.data(), zeros.data());
    HYPRE_IJVectorAssemble(inst.x_ij);
    HYPRE_IJVectorGetObject(inst.x_ij, reinterpret_cast<void **>(&inst.x));

    if (inst.num_interp_vectors > 0)
    {
        const double *z_ptr = mxGetPr(Z);
        inst.interp_ij = static_cast<HYPRE_IJVector *>(std::calloc(inst.num_interp_vectors, sizeof(HYPRE_IJVector)));
        inst.interp_par = static_cast<HYPRE_ParVector *>(std::calloc(inst.num_interp_vectors, sizeof(HYPRE_ParVector)));
        if (!inst.interp_ij || !inst.interp_par)
        {
            mex_error("hypre:alloc", "Failed to allocate interpolation vector arrays.");
        }
        for (int k = 0; k < inst.num_interp_vectors; k++)
        {
            HYPRE_IJVectorCreate(comm, 0, static_cast<HYPRE_BigInt>(n - 1), &inst.interp_ij[k]);
            HYPRE_IJVectorSetObjectType(inst.interp_ij[k], HYPRE_PARCSR);
            HYPRE_IJVectorInitialize(inst.interp_ij[k]);
            const HYPRE_Complex *z_col = reinterpret_cast<const HYPRE_Complex *>(z_ptr + static_cast<mwSize>(k) * n);
            HYPRE_IJVectorSetValues(inst.interp_ij[k], static_cast<HYPRE_Int>(n), rows.data(), z_col);
            HYPRE_IJVectorAssemble(inst.interp_ij[k]);
            HYPRE_IJVectorGetObject(inst.interp_ij[k], reinterpret_cast<void **>(&inst.interp_par[k]));
        }
    }

    HYPRE_BoomerAMGCreate(&inst.solver);

    HYPRE_BoomerAMGSetPrintLevel(inst.solver, o.print_level);
    HYPRE_BoomerAMGSetMaxLevels(inst.solver, o.max_levels);
    HYPRE_BoomerAMGSetCoarsenType(inst.solver, o.coarsen_type);
    HYPRE_BoomerAMGSetAggNumLevels(inst.solver, o.agg_num_levels);
    HYPRE_BoomerAMGSetInterpType(inst.solver, o.interp_type);
    HYPRE_BoomerAMGSetPMaxElmts(inst.solver, o.p_max_elmts);
    HYPRE_BoomerAMGSetStrongThreshold(inst.solver, o.strong_threshold);
    if (o.has_strong_threshold_r)
    {
        HYPRE_BoomerAMGSetStrongThresholdR(inst.solver, o.strong_threshold_r);
    }

    HYPRE_BoomerAMGSetRelaxType(inst.solver, o.relax_type);
    HYPRE_BoomerAMGSetNumSweeps(inst.solver, o.relax_sweeps);
    HYPRE_BoomerAMGSetCycleType(inst.solver, o.cycle_type);
    HYPRE_BoomerAMGSetCycleRelaxType(inst.solver, o.relax_coarse_type, 3);

    HYPRE_BoomerAMGSetNumFunctions(inst.solver, o.num_functions);
    if (o.num_functions > 1)
    {
        std::vector<HYPRE_Int> dof_mapping = build_dof_func(opts, n, o.num_functions);
        HYPRE_Int *mapping_heap = static_cast<HYPRE_Int *>(std::malloc(n * sizeof(HYPRE_Int)));
        if (!mapping_heap)
        {
            mex_error("hypre:alloc", "Failed to allocate dof_func mapping.");
        }
        std::memcpy(mapping_heap, dof_mapping.data(), n * sizeof(HYPRE_Int));
        /* Hypre takes ownership and frees this pointer on solver destroy. */
        HYPRE_BoomerAMGSetDofFunc(inst.solver, mapping_heap);
    }

    HYPRE_BoomerAMGSetNodal(inst.solver, o.nodal);
    HYPRE_BoomerAMGSetNodalDiag(inst.solver, o.nodal_diag);
    HYPRE_BoomerAMGSetNodalLevels(inst.solver, o.nodal_levels);
    HYPRE_BoomerAMGSetKeepSameSign(inst.solver, o.keep_same_sign);

    if (o.has_trunc_factor) { HYPRE_BoomerAMGSetTruncFactor(inst.solver, o.trunc_factor); }
    if (o.has_agg_interp_type) { HYPRE_BoomerAMGSetAggInterpType(inst.solver, o.agg_interp_type); }
    if (o.has_agg_trunc_factor) { HYPRE_BoomerAMGSetAggTruncFactor(inst.solver, o.agg_trunc_factor); }
    if (o.has_agg_p_max_elmts) { HYPRE_BoomerAMGSetAggPMaxElmts(inst.solver, o.agg_p_max_elmts); }
    if (o.has_max_coarse_size) { HYPRE_BoomerAMGSetMaxCoarseSize(inst.solver, o.max_coarse_size); }
    if (o.has_min_coarse_size) { HYPRE_BoomerAMGSetMinCoarseSize(inst.solver, o.min_coarse_size); }

    if (inst.num_interp_vectors > 0)
    {
        HYPRE_BoomerAMGSetInterpVectors(inst.solver, inst.num_interp_vectors, inst.interp_par);
        HYPRE_BoomerAMGSetInterpVecVariant(inst.solver, o.interp_vec_variant);
        HYPRE_BoomerAMGSetInterpVecQMax(inst.solver, o.interp_vec_qmax);
        if (o.has_interp_vec_abs_qtrunc)
        {
            HYPRE_BoomerAMGSetInterpVecAbsQTrunc(inst.solver, o.interp_vec_abs_qtrunc);
        }
        HYPRE_BoomerAMGSetSmoothInterpVectors(inst.solver, o.smooth_interp_vectors);
        HYPRE_BoomerAMGSetInterpRefine(inst.solver, o.interp_refine);
    }

    HYPRE_BoomerAMGSetMaxIter(inst.solver, o.max_iter);
    HYPRE_BoomerAMGSetTol(inst.solver, o.tol);

#ifdef _OPENMP
    double t0 = omp_get_wtime();
#else
    double t0 = 0.0;
#endif
    HYPRE_BoomerAMGSetup(inst.solver, inst.A, inst.b, inst.x);
#ifdef _OPENMP
    double t_setup = omp_get_wtime() - t0;
#else
    double t_setup = mxGetNaN();
#endif
    return t_setup;
}

void command_setup(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 5)
    {
        mex_error("hypre:usage", "setup usage: hypre_boomeramg_mex('setup', A, Z, opts, instance_id)");
    }
    const mxArray *A = prhs[1];
    const mxArray *Z = prhs[2];
    const mxArray *opts = prhs[3];
    std::string id = to_id_string(prhs[4]);

    auto it = g_instances.find(id);
    if (it != g_instances.end())
    {
        destroy_instance(*it->second);
        g_instances.erase(it);
    }

    auto inst = std::make_unique<Instance>();
    double t_setup = setup_instance(*inst, A, Z, opts);
    g_instances[id] = std::move(inst);

    if (!g_at_exit_registered)
    {
        mexAtExit(at_exit);
        g_at_exit_registered = true;
    }
    if (g_instances.size() == 1)
    {
        mexLock();
    }

    if (nlhs > 0)
    {
        const char *fields[] = {"instance_id", "n", "nnz", "num_functions", "num_interp_vectors", "setup_time_seconds"};
        mxArray *out = mxCreateStructMatrix(1, 1, 6, fields);
        mxSetField(out, 0, "instance_id", mxCreateString(id.c_str()));
        mxSetField(out, 0, "n", mxCreateDoubleScalar(static_cast<double>(g_instances[id]->n)));
        mxSetField(out, 0, "nnz", mxCreateDoubleScalar(static_cast<double>(g_instances[id]->nnz)));
        mxSetField(out, 0, "num_functions", mxCreateDoubleScalar(static_cast<double>(g_instances[id]->num_functions)));
        mxSetField(out, 0, "num_interp_vectors", mxCreateDoubleScalar(static_cast<double>(g_instances[id]->num_interp_vectors)));
        mxSetField(out, 0, "setup_time_seconds", mxCreateDoubleScalar(t_setup));
        plhs[0] = out;
    }
}

void command_apply(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 3)
    {
        mex_error("hypre:usage", "apply usage: hypre_boomeramg_mex('apply', rhs, instance_id)");
    }
    if (nlhs > 1)
    {
        mex_error("hypre:usage", "apply returns one output vector.");
    }

    const mxArray *rhs = prhs[1];
    std::string id = to_id_string(prhs[2]);

    auto it = g_instances.find(id);
    if (it == g_instances.end())
    {
        mex_error("hypre:noInstance", "Unknown instance_id. Call setup first.");
    }
    Instance &inst = *it->second;

    if (!mxIsDouble(rhs) || mxIsComplex(rhs))
    {
        mex_error("hypre:badRhs", "rhs must be a real double vector.");
    }
    if (mxGetNumberOfElements(rhs) != inst.n)
    {
        mex_error("hypre:badRhs", "rhs length does not match preconditioner size.");
    }

    const HYPRE_Complex *rhs_ptr = reinterpret_cast<const HYPRE_Complex *>(mxGetPr(rhs));
    HYPRE_IJVectorSetValues(inst.b_ij, static_cast<HYPRE_Int>(inst.n), inst.rows.data(), rhs_ptr);
    HYPRE_IJVectorAssemble(inst.b_ij);

    HYPRE_IJVectorSetConstantValues(inst.x_ij, 0.0);
    HYPRE_IJVectorAssemble(inst.x_ij);

    HYPRE_BoomerAMGSolve(inst.solver, inst.A, inst.b, inst.x);

    if (nlhs > 0)
    {
        plhs[0] = mxCreateDoubleMatrix(inst.n, 1, mxREAL);
        HYPRE_Complex *out = reinterpret_cast<HYPRE_Complex *>(mxGetPr(plhs[0]));
        HYPRE_IJVectorGetValues(inst.x_ij, static_cast<HYPRE_Int>(inst.n), inst.rows.data(), out);
    }
}

void command_clear(int nrhs, const mxArray *prhs[])
{
    if (nrhs == 1)
    {
        cleanup_all();
        while (mexIsLocked()) { mexUnlock(); }
        return;
    }
    if (nrhs != 2)
    {
        mex_error("hypre:usage", "clear usage: hypre_boomeramg_mex('clear') or hypre_boomeramg_mex('clear', instance_id)");
    }
    std::string id = to_id_string(prhs[1]);
    auto it = g_instances.find(id);
    if (it != g_instances.end())
    {
        destroy_instance(*it->second);
        g_instances.erase(it);
    }
    if (g_instances.empty())
    {
        cleanup_all();
        while (mexIsLocked()) { mexUnlock(); }
    }
}

void command_info(int nlhs, mxArray *plhs[])
{
    if (nlhs < 1)
    {
        return;
    }
    const char *fields[] = {"instance_id", "n", "nnz", "num_functions", "num_interp_vectors"};
    mxArray *arr = mxCreateStructMatrix(g_instances.size(), 1, 5, fields);
    mwIndex idx = 0;
    for (const auto &kv : g_instances)
    {
        const Instance &inst = *kv.second;
        mxSetField(arr, idx, "instance_id", mxCreateString(kv.first.c_str()));
        mxSetField(arr, idx, "n", mxCreateDoubleScalar(static_cast<double>(inst.n)));
        mxSetField(arr, idx, "nnz", mxCreateDoubleScalar(static_cast<double>(inst.nnz)));
        mxSetField(arr, idx, "num_functions", mxCreateDoubleScalar(static_cast<double>(inst.num_functions)));
        mxSetField(arr, idx, "num_interp_vectors", mxCreateDoubleScalar(static_cast<double>(inst.num_interp_vectors)));
        idx++;
    }
    plhs[0] = arr;
}

} // namespace

extern "C" void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs < 1 || !mxIsChar(prhs[0]))
    {
        mexErrMsgIdAndTxt("hypre:usage",
                          "Usage: hypre_boomeramg_mex('setup'| 'apply' | 'clear' | 'info', ...)");
    }

    char *cmd_c = mxArrayToString(prhs[0]);
    if (!cmd_c)
    {
        mexErrMsgIdAndTxt("hypre:usage", "Failed to parse command.");
    }
    std::string cmd(cmd_c);
    mxFree(cmd_c);

    try
    {
        if (cmd == "setup")
        {
            command_setup(nlhs, plhs, nrhs, prhs);
            return;
        }
        if (cmd == "apply")
        {
            command_apply(nlhs, plhs, nrhs, prhs);
            return;
        }
        if (cmd == "clear")
        {
            command_clear(nrhs, prhs);
            return;
        }
        if (cmd == "info")
        {
            command_info(nlhs, plhs);
            return;
        }
    }
    catch (const std::exception &e)
    {
        mexErrMsgIdAndTxt("hypre:exception", "%s", e.what());
    }

    mexErrMsgIdAndTxt("hypre:usage", "Unknown command '%s'.", cmd.c_str());
}
