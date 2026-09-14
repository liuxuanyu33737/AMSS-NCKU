#ifndef ABE_WTIME_PROF_H
#define ABE_WTIME_PROF_H

// Coarse-grained MPI_Wtime instrumentation for the ABE evolution path
// (Evolve -> RecursiveStep -> Step).
//
// This header keeps per-MPI-rank accumulators (call count, total time,
// min/max single-call time) for a fixed set of coarse phases and provides a
// single report function that MPI-reduces the statistics once at the end of
// the evolution.
//
// It is intentionally self-contained (no extra object files) and must be
// included from a single translation unit so the file-scope accumulators are
// unique to that TU.  The accumulators are plain static data: on each MPI rank
// this is a separate process, so no inter-rank locking is required.

#include <mpi.h>
#include <cstdio>

#ifndef ENABLE_ABE_PROFILING
#define ENABLE_ABE_PROFILING 0
#endif

#if ENABLE_ABE_PROFILING
// forward declaration of the level-4 surf_MassPAng report (defined in
// surface_integral.C); reduced/printed right after the level-3 report below.
void abe_l4_surf_mass_report(int myrank, int nprocs);

#define ABE_PROF_NUM_PHASES 7

enum AbeProfPhase
{
  ABE_PROF_COMPUTE_RHS = 0, // f_compute_rhs_bssn (and _ss when shell is enabled)
  ABE_PROF_RK4         = 1, // f_rungekutta4_rout / f_rungekutta4_scalar
  ABE_PROF_BOUNDARY    = 2, // sommerfeld / lowerboundset physical-boundary work
  ABE_PROF_ENFORCE_GA  = 3, // f_enforce_ga
  ABE_PROF_SYNC        = 4, // Parallel::Sync (MPI communication)
  ABE_PROF_AMR         = 5, // RestrictProlong / Regrid_Onelevel
  ABE_PROF_OTHER       = 6  // explicitly timed remaining Step work
};

struct AbeProfEntry
{
  const char *name;
  long long count;
  double total;
  double tmin;
  double tmax;
};

static AbeProfEntry g_abe_prof[ABE_PROF_NUM_PHASES] =
{
  {"compute_rhs_bssn", 0, 0.0, 1.0e300, 0.0},
  {"rk4",              0, 0.0, 1.0e300, 0.0},
  {"boundary",         0, 0.0, 1.0e300, 0.0},
  {"enforce_ga",       0, 0.0, 1.0e300, 0.0},
  {"sync",             0, 0.0, 1.0e300, 0.0},
  {"amr",              0, 0.0, 1.0e300, 0.0},
  {"other",            0, 0.0, 1.0e300, 0.0}
};

inline double abe_prof_begin()
{
  return MPI_Wtime();
}

inline void abe_prof_end(int phase, double t0)
{
  const double dt = MPI_Wtime() - t0;
  AbeProfEntry &e = g_abe_prof[phase];
  ++e.count;
  e.total += dt;
  if (dt < e.tmin) e.tmin = dt;
  if (dt > e.tmax) e.tmax = dt;
}

// ---------------------------------------------------------------------------
// Level-2 (fine-grained) instrumentation
//
// The coarse phase ABE_PROF_OTHER still lumps together several unrelated Step
// code blocks.  The fine phases below split it into concrete code segments and
// also time the MPI_Allreduce error checks that were never timed before.  Every
// fine event records its refinement level, so the report can show the level
// distribution as well as per-rank / per-call outliers.
//
// The header is included from a single translation unit (bssn_class.C), so all
// accumulators below are local to that TU, exactly like the coarse ones above.
// ---------------------------------------------------------------------------

#define ABE_FINE_NUM_PHASES 7
#define ABE_FINE_MAXLEV     32

enum AbeFinePhase
{
  ABE_FINE_BH_PORG_RHS     = 0, // compute_Porg_rhs: puncture-position interpolation
  ABE_FINE_BH_EULER_PRED   = 1, // BH Euler predictor step + symmetry/NaN handling
  ABE_FINE_ANALYSIS        = 2, // AnalysisStuff (only when lev == a_lev)
  ABE_FINE_SWAP_PRE_COR    = 3, // swapList(SynchList_pre,SynchList_cor): RK time-level swap
  ABE_FINE_SWAP_STATE_OLD  = 4, // swapList(StateList,...,SynchList_cor): end-of-Step swap
  ABE_FINE_BH_FINAL_COPY   = 5, // BH Porg0 <- Porg1 at the end of each Step
  ABE_FINE_ERR_ALLREDUCE   = 6  // MPI_Allreduce error checks (Step-level global sync)
};

struct AbeFineEntry
{
  const char *name;
  long long count;    // number of events on this MPI rank
  double total;       // accumulated time on this MPI rank
  double tmin;        // shortest single event on this MPI rank
  double tmax;        // longest  single event on this MPI rank
  int    maxlev;      // refinement level of the longest single event
  long long cnt_lev[ABE_FINE_MAXLEV];
  double      tot_lev[ABE_FINE_MAXLEV];
  double      max_lev[ABE_FINE_MAXLEV];
};

static AbeFineEntry g_abe_fine[ABE_FINE_NUM_PHASES] =
{
  {"bh_porg_rhs",     0, 0.0, 1.0e300, 0.0, -1},
  {"bh_euler_pred",   0, 0.0, 1.0e300, 0.0, -1},
  {"analysis",        0, 0.0, 1.0e300, 0.0, -1},
  {"swap_pre_cor",    0, 0.0, 1.0e300, 0.0, -1},
  {"swap_state_old",  0, 0.0, 1.0e300, 0.0, -1},
  {"bh_final_copy",   0, 0.0, 1.0e300, 0.0, -1},
  {"err_allreduce",   0, 0.0, 1.0e300, 0.0, -1}
};

inline double abe_fine_begin()
{
  return MPI_Wtime();
}

inline void abe_fine_update(int phase, int lev, double dt)
{
  if (phase < 0 || phase >= ABE_FINE_NUM_PHASES)
    return;
  AbeFineEntry &e = g_abe_fine[phase];
  ++e.count;
  e.total += dt;
  if (dt < e.tmin) e.tmin = dt;
  if (dt > e.tmax)
  {
    e.tmax = dt;
    e.maxlev = lev;
  }
  const int lv = (lev < 0) ? 0 : ((lev >= ABE_FINE_MAXLEV) ? ABE_FINE_MAXLEV - 1 : lev);
  ++e.cnt_lev[lv];
  e.tot_lev[lv] += dt;
  if (dt > e.max_lev[lv]) e.max_lev[lv] = dt;
}

inline void abe_fine_end(int phase, int lev, double t0)
{
  abe_fine_update(phase, lev, MPI_Wtime() - t0);
}

// Combined end call for code regions that previously produced a single coarse
// 'other' sample: the coarse ABE_PROF_OTHER statistics are kept unchanged and,
// with the very same duration, the matching level-2 event is recorded.
inline void abe_prof_other_fine_end(int fine_phase, int lev, double t0)
{
  const double dt = MPI_Wtime() - t0;
  AbeProfEntry &c = g_abe_prof[ABE_PROF_OTHER];
  ++c.count;
  c.total += dt;
  if (dt < c.tmin) c.tmin = dt;
  if (dt > c.tmax) c.tmax = dt;
  abe_fine_update(fine_phase, lev, dt);
}

inline void abe_fine_report(int myrank, int nprocs)
{
  const int nf   = ABE_FINE_NUM_PHASES;
  const int nlv  = ABE_FINE_MAXLEV;
  const int nrec = nf * nlv;

  long long loc_cnt[nf];
  double loc_total[nf], loc_min[nf], loc_max[nf];
  int loc_maxlev[nf];
  long long loc_cntlev[nrec];
  double loc_totlev[nrec], loc_maxlev_time[nrec];

  for (int i = 0; i < nf; ++i)
  {
    const AbeFineEntry &e = g_abe_fine[i];
    loc_cnt[i] = e.count;
    loc_total[i] = e.total;
    loc_min[i] = (e.count > 0) ? e.tmin : 0.0;
    loc_max[i] = e.tmax;
    loc_maxlev[i] = e.maxlev;
    for (int l = 0; l < nlv; ++l)
    {
      loc_cntlev[i * nlv + l]      = e.cnt_lev[l];
      loc_totlev[i * nlv + l]      = e.tot_lev[l];
      loc_maxlev_time[i * nlv + l] = e.max_lev[l];
    }
  }

  long long gl_cnt_sum[nf];
  double gl_total_sum[nf], gl_total_max[nf], gl_min[nf], gl_max[nf];
  long long gl_cntlev[nrec];
  double gl_totlev[nrec], gl_maxlev_time[nrec];

  MPI_Reduce(loc_cnt,   gl_cnt_sum,   nf,   MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_total, gl_total_sum, nf,   MPI_DOUBLE,        MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_min,   gl_min,       nf,   MPI_DOUBLE,        MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_max,   gl_max,       nf,   MPI_DOUBLE,        MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_cntlev,      gl_cntlev,      nrec, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_totlev,      gl_totlev,      nrec, MPI_DOUBLE,        MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_maxlev_time, gl_maxlev_time, nrec, MPI_DOUBLE,        MPI_MAX, 0, MPI_COMM_WORLD);

  double *all_total = new double[nprocs * nf];
  int *all_maxlev = new int[nprocs * nf];
  MPI_Gather(loc_total,  nf, MPI_DOUBLE, all_total,  nf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(loc_maxlev, nf, MPI_INT,    all_maxlev, nf, MPI_INT,    0, MPI_COMM_WORLD);

  if (myrank == 0)
  {
    printf("\n===== ABE MPI_Wtime level-2 profile: breakdown of the coarse 'other' phase (nprocs = %d) =====\n", nprocs);
    printf("%-16s %12s %14s %14s %14s %14s %14s %7s\n",
           "name", "count", "total(s)", "avg(s)", "min(s)", "max(s)", "maxrank(s)", "maxlev");
    for (int i = 0; i < nf; ++i)
    {
      if (gl_cnt_sum[i] <= 0)
        continue;
      const double avg = gl_total_sum[i] / (double)gl_cnt_sum[i];
      gl_total_max[i] = 0.0;
      int wrank = 0;
      for (int r = 0; r < nprocs; ++r)
      {
        if (all_total[r * nf + i] > gl_total_max[i])
        {
          gl_total_max[i] = all_total[r * nf + i];
          wrank = r;
        }
      }
      const int wlev = (wrank >= 0 && wrank < nprocs) ? all_maxlev[wrank * nf + i] : -1;
      printf("%-16s %12lld %14.6f %14.6e %14.6e %14.6e %14.6f %7d\n",
             g_abe_fine[i].name, gl_cnt_sum[i], gl_total_sum[i], avg,
             gl_min[i], gl_max[i], gl_total_max[i], wlev);
    }
    printf("\nPer-level detail (count / total(s) / avg(s) / max(s) of each fine phase):\n");
    for (int i = 0; i < nf; ++i)
    {
      if (gl_cnt_sum[i] <= 0)
        continue;
      printf("  [%s]\n", g_abe_fine[i].name);
      for (int l = 0; l < nlv; ++l)
      {
        if (gl_cntlev[i * nlv + l] <= 0)
          continue;
        const double lavg = gl_totlev[i * nlv + l] / (double)gl_cntlev[i * nlv + l];
        printf("    lev %2d : count %12lld  total %14.6f  avg %14.6e  max %14.6e\n",
               l, gl_cntlev[i * nlv + l], gl_totlev[i * nlv + l], lavg,
               gl_maxlev_time[i * nlv + l]);
      }
    }
    printf("===== end ABE MPI_Wtime level-2 profile =====\n\n");
  }

  delete[] all_total;
  delete[] all_maxlev;
}

// ---------------------------------------------------------------------------
// Level-3 instrumentation: breakdown of one 'analysis' event (AnalysisStuff)
//
// Level-2 (ABE_FINE_ANALYSIS) already times the whole AnalysisStuff call from
// the Step side.  The phases below split the inside of AnalysisStuff so the
// ~10 s long tail reported for 'analysis' can be attributed:
//
//   compute_psi4   - whole-grid EB-method Weyl-scalar computation
//                    (bssn_class::Compute_Psi4; contains a Parallel::Sync)
//   surf_wave      - per-extraction-radius psi4 spherical-harmonic expansion
//                    (Waveshell->surf_Wave; has an internal MPI_Allreduce)
//   surf_masspang  - per-extraction-radius ADM mass/linear/angular momentum
//                    (Waveshell->surf_MassPAng; internal MPI_Allreduce)
//   psi4_write     - Psi4Monitor->writefile at each radius (rank-local I/O)
//   map_write      - MAPMonitor->writefile at each radius  (rank-local I/O)
//   bh_write       - BHMonitor->writefile of the puncture positions (I/O)
//
// Recording layout is identical to level-2 (per-rank count / total / min /
// max / maxlev / per-level arrays); abe_l3_report() reduces it once at the end
// of the evolution.  It only adds two MPI_Wtime() calls around existing
// function calls, so it neither changes the calculation nor the MPI order.
// ---------------------------------------------------------------------------

#define ABE_L3_NUM_PHASES 6
#define ABE_L3_MAXLEV     ABE_FINE_MAXLEV

enum AbeL3Phase
{
  ABE_L3_COMPUTE_PSI4 = 0, // bssn_class::Compute_Psi4 (whole-grid Weyl scalar)
  ABE_L3_SURF_WAVE    = 1, // per-radius surf_Wave      (internal MPI_Allreduce)
  ABE_L3_SURF_MASS    = 2, // per-radius surf_MassPAng  (internal MPI_Allreduce)
  ABE_L3_PSI4_WRITE   = 3, // Psi4Monitor->writefile per radius (rank-local I/O)
  ABE_L3_MAP_WRITE    = 4, // MAPMonitor->writefile per radius  (rank-local I/O)
  ABE_L3_BH_WRITE     = 5  // BHMonitor->writefile per analysis (rank-local I/O)
};

static AbeFineEntry g_abe_l3[ABE_L3_NUM_PHASES] =
{
  {"compute_psi4",   0, 0.0, 1.0e300, 0.0, -1},
  {"surf_wave",      0, 0.0, 1.0e300, 0.0, -1},
  {"surf_masspang",  0, 0.0, 1.0e300, 0.0, -1},
  {"psi4_write",     0, 0.0, 1.0e300, 0.0, -1},
  {"map_write",      0, 0.0, 1.0e300, 0.0, -1},
  {"bh_write",       0, 0.0, 1.0e300, 0.0, -1}
};

inline double abe_l3_begin()
{
  return MPI_Wtime();
}

inline void abe_l3_update(int phase, int lev, double dt)
{
  if (phase < 0 || phase >= ABE_L3_NUM_PHASES)
    return;
  AbeFineEntry &e = g_abe_l3[phase];
  ++e.count;
  e.total += dt;
  if (dt < e.tmin) e.tmin = dt;
  if (dt > e.tmax)
  {
    e.tmax = dt;
    e.maxlev = lev;
  }
  const int lv = (lev < 0) ? 0 : ((lev >= ABE_L3_MAXLEV) ? ABE_L3_MAXLEV - 1 : lev);
  ++e.cnt_lev[lv];
  e.tot_lev[lv] += dt;
  if (dt > e.max_lev[lv]) e.max_lev[lv] = dt;
}

inline void abe_l3_end(int phase, int lev, double t0)
{
  abe_l3_update(phase, lev, MPI_Wtime() - t0);
}

inline void abe_l3_report(int myrank, int nprocs)
{
  const int nf   = ABE_L3_NUM_PHASES;
  const int nlv  = ABE_L3_MAXLEV;
  const int nrec = nf * nlv;

  long long loc_cnt[nf];
  double loc_total[nf], loc_min[nf], loc_max[nf];
  int loc_maxlev[nf];
  long long loc_cntlev[nrec];
  double loc_totlev[nrec], loc_maxlev_time[nrec];

  for (int i = 0; i < nf; ++i)
  {
    const AbeFineEntry &e = g_abe_l3[i];
    loc_cnt[i] = e.count;
    loc_total[i] = e.total;
    loc_min[i] = (e.count > 0) ? e.tmin : 0.0;
    loc_max[i] = e.tmax;
    loc_maxlev[i] = e.maxlev;
    for (int l = 0; l < nlv; ++l)
    {
      loc_cntlev[i * nlv + l]      = e.cnt_lev[l];
      loc_totlev[i * nlv + l]      = e.tot_lev[l];
      loc_maxlev_time[i * nlv + l] = e.max_lev[l];
    }
  }

  long long gl_cnt_sum[nf];
  double gl_total_sum[nf], gl_total_max[nf], gl_min[nf], gl_max[nf];
  long long gl_cntlev[nrec];
  double gl_totlev[nrec], gl_maxlev_time[nrec];

  MPI_Reduce(loc_cnt,   gl_cnt_sum,   nf,   MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_total, gl_total_sum, nf,   MPI_DOUBLE,        MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_min,   gl_min,       nf,   MPI_DOUBLE,        MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_max,   gl_max,       nf,   MPI_DOUBLE,        MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_cntlev,      gl_cntlev,      nrec, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_totlev,      gl_totlev,      nrec, MPI_DOUBLE,        MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_maxlev_time, gl_maxlev_time, nrec, MPI_DOUBLE,        MPI_MAX, 0, MPI_COMM_WORLD);

  double *all_total = new double[nprocs * nf];
  int *all_maxlev = new int[nprocs * nf];
  MPI_Gather(loc_total,  nf, MPI_DOUBLE, all_total,  nf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(loc_maxlev, nf, MPI_INT,    all_maxlev, nf, MPI_INT,    0, MPI_COMM_WORLD);

  if (myrank == 0)
  {
    printf("\n===== ABE MPI_Wtime level-3 profile: breakdown of analysis (AnalysisStuff) (nprocs = %d) =====\n", nprocs);
    printf("%-16s %12s %14s %14s %14s %14s %14s %7s\n",
           "name", "count", "total(s)", "avg(s)", "min(s)", "max(s)", "maxrank(s)", "maxlev");
    for (int i = 0; i < nf; ++i)
    {
      if (gl_cnt_sum[i] <= 0)
        continue;
      const double avg = gl_total_sum[i] / (double)gl_cnt_sum[i];
      gl_total_max[i] = 0.0;
      int wrank = 0;
      for (int r = 0; r < nprocs; ++r)
      {
        if (all_total[r * nf + i] > gl_total_max[i])
        {
          gl_total_max[i] = all_total[r * nf + i];
          wrank = r;
        }
      }
      const int wlev = (wrank >= 0 && wrank < nprocs) ? all_maxlev[wrank * nf + i] : -1;
      printf("%-16s %12lld %14.6f %14.6e %14.6e %14.6e %14.6f %7d\n",
             g_abe_l3[i].name, gl_cnt_sum[i], gl_total_sum[i], avg,
             gl_min[i], gl_max[i], gl_total_max[i], wlev);
    }
    printf("\nPer-level detail (count / total(s) / avg(s) / max(s) of each level-3 phase):\n");
    for (int i = 0; i < nf; ++i)
    {
      if (gl_cnt_sum[i] <= 0)
        continue;
      printf("  [%s]\n", g_abe_l3[i].name);
      for (int l = 0; l < nlv; ++l)
      {
        if (gl_cntlev[i * nlv + l] <= 0)
          continue;
        const double lavg = gl_totlev[i * nlv + l] / (double)gl_cntlev[i * nlv + l];
        printf("    lev %2d : count %12lld  total %14.6f  avg %14.6e  max %14.6e\n",
               l, gl_cntlev[i * nlv + l], gl_totlev[i * nlv + l], lavg,
               gl_maxlev_time[i * nlv + l]);
      }
    }
    printf("===== end ABE MPI_Wtime level-3 profile =====\n\n");
  }

  delete[] all_total;
  delete[] all_maxlev;
}

inline void abe_prof_report(int myrank, int nprocs)
{
  const int n = ABE_PROF_NUM_PHASES;

  long long loc_count[n];
  double loc_total[n];
  double loc_min[n];
  double loc_max[n];

  for (int i = 0; i < n; ++i)
  {
    loc_count[i] = g_abe_prof[i].count;
    loc_total[i] = g_abe_prof[i].total;
    loc_min[i]   = (g_abe_prof[i].count > 0) ? g_abe_prof[i].tmin : 0.0;
    loc_max[i]   = g_abe_prof[i].tmax;
  }

  long long gl_count_sum[n];
  long long gl_count_max[n];
  double gl_total_sum[n];
  double gl_total_max[n];
  double gl_min[n];
  double gl_max[n];

  MPI_Reduce(loc_count, gl_count_sum, n, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_count, gl_count_max, n, MPI_LONG_LONG_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_total, gl_total_sum, n, MPI_DOUBLE,        MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_total, gl_total_max, n, MPI_DOUBLE,        MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_min,   gl_min,       n, MPI_DOUBLE,        MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_max,   gl_max,       n, MPI_DOUBLE,        MPI_MAX, 0, MPI_COMM_WORLD);

  if (myrank == 0)
  {
    printf("\n===== ABE MPI_Wtime coarse-grained profile (nprocs = %d) =====\n", nprocs);
    printf("%-18s %12s %14s %14s %14s %14s %14s\n",
           "phase", "count(sum)", "total(s)", "avg(s)", "min(s)", "max(s)", "maxrank(s)");
    for (int i = 0; i < n; ++i)
    {
      const double avg = (gl_count_sum[i] > 0) ? (gl_total_sum[i] / (double)gl_count_sum[i]) : 0.0;
      printf("%-18s %12lld %14.6f %14.6e %14.6e %14.6e %14.6f\n",
             g_abe_prof[i].name,
             gl_count_sum[i],
             gl_total_sum[i],
             avg,
             gl_min[i],
             gl_max[i],
             gl_total_max[i]);
    }
    printf("===== end ABE MPI_Wtime coarse-grained profile =====\n\n");
  }
  abe_fine_report(myrank, nprocs);
  abe_l3_report(myrank, nprocs);
  abe_l4_surf_mass_report(myrank, nprocs);
}

#else

enum AbeProfPhase
{
  ABE_PROF_COMPUTE_RHS = 0, ABE_PROF_RK4, ABE_PROF_BOUNDARY,
  ABE_PROF_ENFORCE_GA, ABE_PROF_SYNC, ABE_PROF_AMR, ABE_PROF_OTHER
};

enum AbeFinePhase
{
  ABE_FINE_BH_PORG_RHS = 0, ABE_FINE_BH_EULER_PRED, ABE_FINE_ANALYSIS,
  ABE_FINE_SWAP_PRE_COR, ABE_FINE_SWAP_STATE_OLD, ABE_FINE_BH_FINAL_COPY,
  ABE_FINE_ERR_ALLREDUCE
};

enum AbeL3Phase
{
  ABE_L3_COMPUTE_PSI4 = 0, ABE_L3_SURF_WAVE, ABE_L3_SURF_MASS,
  ABE_L3_PSI4_WRITE, ABE_L3_MAP_WRITE, ABE_L3_BH_WRITE
};

#define abe_prof_begin()                         (0.0)
#define abe_prof_end(phase, t0)                  ((void)0)
#define abe_fine_begin()                         (0.0)
#define abe_fine_end(phase, lev, t0)             ((void)0)
#define abe_prof_other_fine_end(phase, lev, t0)  ((void)0)
#define abe_l3_begin()                           (0.0)
#define abe_l3_end(phase, lev, t0)               ((void)0)
#define abe_prof_report(myrank, nprocs)          ((void)0)

#endif

#endif // ABE_WTIME_PROF_H
