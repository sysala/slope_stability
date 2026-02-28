# `sparsersb` Usage and Prebuilt `B' * diag(d) * B` Assembly

This note focuses on two things:

1. How to use `sparsersb` effectively in Octave.
2. How to prebuild indexing so repeated `A(d) = B' * diag(d) * B` assembly is fast when `B` is fixed and only `d` changes.

## 1. `sparsersb` basics

`sparsersb` is an Octave package backed by `librsb` (recursive sparse block format with multithreaded sparse kernels).

Load package:

```octave
pkg load sparsersb
```

Create RSB matrices:

```octave
% from native sparse
A_rsb = sparsersb(A_sparse);

% from triplets (i,j,v)
A_rsb = sparsersb(iA, jA, vA, m, n);
```

Use like regular sparse:

```octave
y = A_rsb * x;       % SpMV
C = A_rsb * B_rsb;   % SpMM
```

Runtime threading:

```bash
OMP_NUM_THREADS=16 benchmark_octave/local/bin/octave-rsb --quiet
```

## 2. Problem setup for repeated `B' * diag(d) * B`

Given fixed sparse `B` (`m x n`) and changing vectors `d` (`m x 1`), you need many assemblies:

```text
A(d) = B' * diag(d) * B
```

Do not rebuild the sparsity pattern every iteration. Precompute row-wise pair contributions once, then reuse.

## 3. Prebuild indexing (one-time)

```octave
function P = prepare_BtDB_pattern_upper(B)
  [iB, jB, vB] = find(B);
  m = rows(B);
  n = columns(B);

  [iB, p] = sort(iB);
  jB = jB(p);
  vB = vB(p);

  row_counts = accumarray(iB, 1, [m, 1]);
  row_ptr = [1; cumsum(row_counts) + 1];

  pair_counts = row_counts .* (row_counts + 1) / 2;
  total_pairs = sum(pair_counts);

  Iraw = zeros(total_pairs, 1);
  Jraw = zeros(total_pairs, 1);
  Rraw = zeros(total_pairs, 1);
  Wraw = zeros(total_pairs, 1);

  pos = 1;
  for r = 1:m
    k = row_counts(r);
    if k == 0, continue; end

    idx = row_ptr(r):(row_ptr(r + 1) - 1);
    cols = jB(idx);
    vals = vB(idx);

    [CJ, CI] = meshgrid(cols, cols);
    [VJ, VI] = meshgrid(vals, vals);
    mask = triu(true(k));  % upper triangle only

    t = nnz(mask);
    rr = pos:(pos + t - 1);
    Iraw(rr) = CI(mask);
    Jraw(rr) = CJ(mask);
    Rraw(rr) = r;
    Wraw(rr) = VI(mask) .* VJ(mask);
    pos = pos + t;
  end

  Iraw = Iraw(1:pos-1);
  Jraw = Jraw(1:pos-1);
  Rraw = Rraw(1:pos-1);
  Wraw = Wraw(1:pos-1);

  [IJu, ~, g] = unique([Iraw Jraw], "rows");
  Iu = IJu(:, 1);
  Ju = IJu(:, 2);

  off = find(Iu ~= Ju);

  P = struct();
  P.n = n;
  P.Iu = Iu;
  P.Ju = Ju;
  P.g = g;
  P.Rraw = Rraw;
  P.Wraw = Wraw;
  P.off = off;
  P.Ifull = [Iu; Ju(off)];
  P.Jfull = [Ju; Iu(off)];
end
```

## 4. Repeated assembly (timed hot path)

The measured hot lines are:

```octave
vals_pair = P.Wraw .* d(P.Rraw);
vals_u = accumarray(P.g, vals_pair, [numel(P.Iu), 1]);
vals_full = [vals_u; vals_u(P.off)];
A = sparsersb(P.Ifull, P.Jfull, vals_full, P.n, P.n);
```

Recommended faster constructor for symmetric unique upper-triangle input:

```octave
vals_pair = P.Wraw .* d(P.Rraw);
vals_u = accumarray(P.g, vals_pair, [numel(P.Iu), 1]);
A = sparsersb(P.Iu, P.Ju, vals_u, P.n, P.n, "unique", "sym");
```

`"unique"` tells constructor that `(i,j)` entries are already unique (no duplicate summation needed).  
`"sym"` tells constructor matrix is symmetric (upper triangle provided).

## 5. When this strategy helps

Use prebuilt indexing when:

- `B` is fixed for many Newton/continuation steps.
- `d` changes often.
- each row of `B` has small support (`k` small), so precompute cost `sum(k_r^2)` stays manageable.

If rows are very dense, precompute memory/work can blow up.

## 6. Bench scripts in this repo

Line-by-line assembly timing (`full` vs `sym_unique`):

```bash
OMP_NUM_THREADS=16 benchmark_octave/local/bin/octave-rsb --quiet --eval \
"addpath('benchmark_octave'); \
 run_btDb_prebuilt_line_breakdown('profile','medium','constructor','full','repeats',5,'warmup',2); \
 run_btDb_prebuilt_line_breakdown('profile','medium','constructor','sym_unique','repeats',5,'warmup',2);"
```

Full kernel benchmark with `sparsersb` backend:

```bash
OMP_NUM_THREADS=16 benchmark_octave/local/bin/octave-rsb --quiet --eval \
"addpath('benchmark_octave'); \
 run_benchmarks('profile','medium','sparse_backend','sparsersb','repeats',3,'warmup',1);"
```

## 7. Practical recommendation

For your repeated tangent assembly, use prebuilt indexing with:

```octave
A = sparsersb(P.Iu, P.Ju, vals_u, P.n, P.n, "unique", "sym");
```

This avoids duplicate reduction in constructor and avoids explicit lower-triangle materialization.
