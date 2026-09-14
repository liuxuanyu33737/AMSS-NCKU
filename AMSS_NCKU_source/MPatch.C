
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <new>
#include <vector>
#include <climits>
using namespace std;

#include "misc.h"
#include "MPatch.h"
#include "Parallel.h"
#include "fmisc.h"
#include "v7_wave_profile.h"

#ifndef VERIFY_BATCH17
#define VERIFY_BATCH17 0
#endif
#ifndef ENABLE_ABE_PROFILING
#define ENABLE_ABE_PROFILING 0
#endif
#if ENABLE_ABE_PROFILING
#define ABE_PROFILE_WTIME() (MPI_Wtime)()
//|===========================================================================
//| Level-5 profiling of Patch::Interp_Points(): MPI_Wtime phase timers inside
//| the two Interp_Points() overloads below.  Profiling only - no algorithm or
//| numerical change.  A per-rank mark is armed by the profiled caller
//| (surf_MassPAng cgh overload in surface_integral.C); only marked calls take
//| the samples.  All other callers of Interp_Points() (surf_Wave, Parallel
//| prolongation, ...) skip the timing completely.  Statistics are merged and
//| printed at the end of the run from abe_l5_interp_report(), invoked by
//| abe_l4_surf_mass_report().
//|===========================================================================

enum
{
  ABE_L5_SETUP = 0,   // entry -> start of the point loop (allocs, DH/llb/uub)
  ABE_L5_POINT_LOOP,  // whole per-point loop over the NN requested points
  ABE_L5_SEARCH,      // derived residual POINT_LOOP - INTERP (block scan, per-
                      // point prep, points owned by remote patches)
  ABE_L5_INTERP,      // f_global_interp() kernels on locally owned points
  ABE_L5_AR_SHELLF,   // MPI_Allreduce of shellf (surface values)
  ABE_L5_AR_WEIGHT,   // MPI_Allreduce of weight (per-point ownership)
  ABE_L5_POST,        // weight normalization / final Shellf fill
  ABE_L5_CLEANUP,     // delete[] of local temporaries
  ABE_L5_TOTAL,       // whole Interp_Points call (entry->exit MPI_Wtime)
  ABE_L5_NPHASE
};

static int abe_l5_mark_pending = 0;     // marked Interp_Points calls queued here
static long long abe_l5_calls = 0;      // profiled Interp_Points executions
static long long abe_l5_pts_total = 0;  // NN of profiled calls (pre-Allreduce)
static long long abe_l5_pts_local = 0;  // points interpolated locally
static double abe_l5s_cnt[ABE_L5_NPHASE] = {0.0};
static double abe_l5s_tot[ABE_L5_NPHASE] = {0.0};
static double abe_l5s_mn[ABE_L5_NPHASE];
static double abe_l5s_mx[ABE_L5_NPHASE] = {0.0};

static void abe_l5s_rec(const int ph, const double dt)
{
  if (ph < 0 || ph >= ABE_L5_NPHASE || dt < 0.0)
    return;
  abe_l5s_cnt[ph] += 1.0;
  abe_l5s_tot[ph] += dt;
  if (abe_l5s_cnt[ph] == 1.0)
    abe_l5s_mn[ph] = dt;
  else if (dt < abe_l5s_mn[ph])
    abe_l5s_mn[ph] = dt;
  if (dt > abe_l5s_mx[ph])
    abe_l5s_mx[ph] = dt;
}
// Caller-side: arm profiling for the next Interp_Points() call on this rank.
void abe_l5_interp_mark(void)
{
  abe_l5_mark_pending++;
}

static int abe_l5_interp_consume(void)
{
  if (abe_l5_mark_pending <= 0)
    return 0;
  abe_l5_mark_pending--;
  abe_l5_calls++;
  return 1;
}

//|===========================================================================
//| Level-6 profiling of f_global_interp() (Fortran global_interp_ in
//| fmisc.f90, Cell-center build).  The Fortran module abe_l6_gi_prof times
//| the internal phases of every global_interp() call while abe_l6_gi_arm(1)
//| is in force; the two Interp_Points() overloads below arm/disarm it only
//| for the same Level-5-marked calls coming from surf_MassPAng, so all other
//| callers of f_global_interp() keep near-zero overhead.  The accumulators
//| are fetched here, reduced over MPI_COMM_WORLD and printed on rank 0 right
//| after the Level-5 table.  Profiling only - no algorithm/numerical change.
//|===========================================================================

extern "C" void abe_l6_gi_arm(int on);
extern "C" void abe_l6_gi_get(int ph, double *cnt, double *tot, double *mn, double *mx);
extern "C" void abe_l7_p3_arm(int on);
extern "C" void abe_l7_p3_get(int ph, double *cnt, double *tot, double *mn, double *mx);
extern "C" void abe_l7_p3_get_polint(double *cnt);

enum
{
  ABE_L7_P3_PREP = 0, // polin3() entry/setup before the first polint block
  ABE_L7_P3_ZDIR,     // polint along x3a (3rd axis / z): m*n calls per polin3
  ABE_L7_P3_YDIR,     // polint along x2a (2nd axis / y): m calls per polin3
  ABE_L7_P3_XDIR,     // polint along x1a (1st axis / x): 1 call per polin3
  ABE_L7_P3_TOTAL,    // whole polin3() call
  ABE_L7_P3_NPHASE
};

// Fetch/reduce/print the Level-7 polin3() phase accumulators (prototype;
// definition below, called at the end of abe_l6_gi_report()).
static void abe_l7_p3_report(int myrank, int nprocs);

enum
{
  ABE_L6_GI_INDEX = 0, // index/bounds prep + box clamp + normalized coords
  ABE_L6_GI_DECIDE3D,  // decide3d(): stencil fetch (symmetry / SoA reflect)
  ABE_L6_GI_POLIN3,    // polin3(): 3-D polynomial interpolation (polint nest)
  ABE_L6_GI_TOTAL,     // whole f_global_interp() call
  ABE_L6_GI_NPHASE
};

// Fetch the Fortran accumulators, reduce over ranks, print on rank 0.
static void abe_l6_gi_report(int myrank, int nprocs)
{
  const int n = ABE_L6_GI_NPHASE;

  double loc_cnt[n], gl_cnt[n], loc_tot[n], gl_tot[n];
  double loc_mn[n], gl_mn[n], loc_mx[n], gl_mx[n], gl_totmax[n];

  for (int i = 0; i < n; ++i)
  {
    double cnt = 0.0, tot = 0.0, mn = 1.0e300, mx = 0.0;
    abe_l6_gi_get(i + 1, &cnt, &tot, &mn, &mx); // Fortran phases are 1-based
    loc_cnt[i] = cnt;
    loc_tot[i] = tot;
    loc_mn[i]  = (cnt > 0.0) ? mn : 1.0e300;
    loc_mx[i]  = (cnt > 0.0) ? mx : 0.0;
  }

  MPI_Reduce(loc_cnt, gl_cnt, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_tot, gl_tot, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_tot, gl_totmax, n, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_mn, gl_mn, n, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_mx, gl_mx, n, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (myrank == 0)
  {
    static const char *l6_name[n] =
    {
      "gi_index",     // index/bounds prep, box clamp, coordinates
      "gi_decide3d",  // decide3d() stencil fetch (symmetry / SoA)
      "gi_polin3",    // polin3() interpolation core (nested polint)
      "gi_total"      // whole f_global_interp() call
    };
    printf("\n===== Level-6 profile: f_global_interp() phases (surf_MassPAng-marked Interp_Points) =====\n");
    printf("nprocs=%d  counts are per f_global_interp() call\n", nprocs);
    printf("%-22s %12s %14s %14s %14s %14s %14s\n",
           "phase", "count(sum)", "total(s)", "avg(s)", "min(s)", "max(s)", "maxrank(s)");
    for (int i = 0; i < n; ++i)
    {
      if (gl_cnt[i] <= 0.0)
        continue;
      const double avg = gl_tot[i] / gl_cnt[i];
      printf("%-22s %12.0f %14.6f %14.6e %14.6e %14.6e %14.6f\n",
             l6_name[i], gl_cnt[i], gl_tot[i], avg, gl_mn[i], gl_mx[i], gl_totmax[i]);
    }
    printf("===== end Level-6 profile =====\n\n");
  }
  // Level-7: same surf_MassPAng-marked window, one level deeper - the
  // polin3() internal phases inside f_global_interp().  All ranks enter here
  // together, right after the Level-6 table.
  abe_l7_p3_report(myrank, nprocs);
}
//|===========================================================================
//| Level-7 profiling of polin3() (Fortran module abe_l7_p3_prof in fmisc.f90,
//| Cell-center build).  polin3() follows zyx order: it interpolates along x3a
//| (third axis / z) on m*n grid lines, then along x2a (second axis / y) on m
//| lines, then a single final polint along x1a (first axis / x).  The module
//| records each block's wall time per polin3() call while abe_l7_p3_arm(1) is
//| in force - the same surf_MassPAng-marked window as Level-6 - plus the total
//| polint() invocation count, so per-call overhead can be separated from the
//| polint() body cost.  polin3()/polint() are algorithmically untouched.
//|===========================================================================
static void abe_l7_p3_report(int myrank, int nprocs)
{
  const int n = ABE_L7_P3_NPHASE;

  double loc_cnt[n], gl_cnt[n], loc_tot[n], gl_tot[n];
  double loc_mn[n], gl_mn[n], loc_mx[n], gl_mx[n], gl_totmax[n];
  double loc_plnt = 0.0, gl_plnt = 0.0;

  for (int i = 0; i < n; ++i)
  {
    double cnt = 0.0, tot = 0.0, mn = 1.0e300, mx = 0.0;
    abe_l7_p3_get(i + 1, &cnt, &tot, &mn, &mx); // Fortran phases are 1-based
    loc_cnt[i] = cnt;
    loc_tot[i] = tot;
    loc_mn[i]  = (cnt > 0.0) ? mn : 1.0e300;
    loc_mx[i]  = (cnt > 0.0) ? mx : 0.0;
  }
  abe_l7_p3_get_polint(&loc_plnt);

  MPI_Reduce(loc_cnt, gl_cnt, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_tot, gl_tot, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_tot, gl_totmax, n, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_mn, gl_mn, n, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_mx, gl_mx, n, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&loc_plnt, &gl_plnt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (myrank == 0)
  {
    static const char *l7_name[n] =
    {
      "p3_prep",    // polin3() entry/setup before the first polint block
      "p3_z_x3a",   // polint along x3a (3rd axis / z): m*n calls per polin3
      "p3_y_x2a",   // polint along x2a (2nd axis / y): m calls per polin3
      "p3_x_x1a",   // polint along x1a (1st axis / x): 1 call per polin3
      "p3_total"    // whole polin3() call
    };
    const double gl_p3calls = gl_cnt[ABE_L7_P3_TOTAL];
    printf("\n===== Level-7 profile: polin3() internal phases (surf_MassPAng-marked Interp_Points) =====\n");
    printf("nprocs=%d  counts are per polin3() call", nprocs);
    if (gl_p3calls > 0.0)
      printf("  polint_calls(sum)=%.0f  polint_calls/polin3(avg)=%.4f (order 6 -> 36+6+1=43 static)",
             gl_plnt, gl_plnt / gl_p3calls);
    printf("\n");
    printf("%-22s %12s %14s %14s %14s %14s %14s\n",
           "phase", "count(sum)", "total(s)", "avg(s)", "min(s)", "max(s)", "maxrank(s)");
    for (int i = 0; i < n; ++i)
    {
      if (gl_cnt[i] <= 0.0)
        continue;
      const double avg = gl_tot[i] / gl_cnt[i];
      printf("%-22s %12.0f %14.6f %14.6e %14.6e %14.6e %14.6f\n",
             l7_name[i], gl_cnt[i], gl_tot[i], avg, gl_mn[i], gl_mx[i], gl_totmax[i]);
    }
    printf("===== end Level-7 profile =====\n\n");
  }
}

// Merge the per-rank Level-5 accumulators and print the report on rank 0.
// All MPI ranks call this together (invoked from abe_l4_surf_mass_report(),
// i.e. once at the end of the run, after the Level-4 report).
void abe_l5_interp_report(int myrank, int nprocs)
{
  const int n = ABE_L5_NPHASE;

  long long loc_cnt[n], gl_cnt[n];
  double loc_tot[n], gl_tot[n], loc_mn[n], gl_mn[n];
  double loc_mx[n], gl_mx[n], gl_totmax[n];

  for (int i = 0; i < n; ++i)
  {
    loc_cnt[i] = (long long)abe_l5s_cnt[i];
    loc_tot[i] = abe_l5s_tot[i];
    loc_mn[i]  = (abe_l5s_cnt[i] > 0.0) ? abe_l5s_mn[i] : 1.0e300;
    loc_mx[i]  = (abe_l5s_cnt[i] > 0.0) ? abe_l5s_mx[i] : 0.0;
  }

  long long loc_calls = abe_l5_calls, gl_calls = 0;
  long long loc_totpts = abe_l5_pts_total, gl_totpts = 0;
  long long loc_locpts = abe_l5_pts_local, gl_locpts = 0;

  MPI_Reduce(loc_cnt, gl_cnt, n, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_tot, gl_tot, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_tot, gl_totmax, n, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_mn, gl_mn, n, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_mx, gl_mx, n, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&loc_calls, &gl_calls, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&loc_totpts, &gl_totpts, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&loc_locpts, &gl_locpts, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  if (myrank == 0)
  {
    static const char *l5_name[n] =
    {
      "interp_setup",
      "interp_point_loop",
      "interp_block_search",
      "interp_global_interp",
      "interp_ar_shellf",
      "interp_ar_weight",
      "interp_post",
      "interp_cleanup",
      "interp_total"
    };
    printf("\n===== Level-5 profile: Patch::Interp_Points phases (surf_MassPAng-marked calls) =====\n");
    printf("nprocs=%d  marked_interp_calls(sum)=%lld  points_fed(sum)=%lld  points_interp_local(sum)=%lld\n",
           nprocs, gl_calls, gl_totpts, gl_locpts);
    printf("%-22s %12s %14s %14s %14s %14s %14s\n",
           "phase", "count(sum)", "total(s)", "avg(s)", "min(s)", "max(s)", "maxrank(s)");
    for (int i = 0; i < n; ++i)
    {
      if (gl_cnt[i] <= 0)
        continue;
      const double avg = gl_tot[i] / (double)gl_cnt[i];
      printf("%-22s %12lld %14.6f %14.6e %14.6e %14.6e %14.6f\n",
             l5_name[i], gl_cnt[i], gl_tot[i], avg, gl_mn[i], gl_mx[i], gl_totmax[i]);
    }
    printf("===== end Level-5 profile =====\n\n");
  }
  // Level-6: fetch the Fortran-side f_global_interp() phase accumulators,
  // reduce over ranks and print directly after the Level-5 table.
  abe_l6_gi_report(myrank, nprocs);
}

#else

#define ABE_PROFILE_WTIME()       (0.0)
#define abe_l5_interp_consume()   (0)
#define abe_l5s_rec(phase, dt)    ((void)0)
#define abe_l6_gi_arm(on)         ((void)0)
#define abe_l7_p3_arm(on)         ((void)0)
static long long abe_l5_pts_total = 0;
static long long abe_l5_pts_local = 0;
enum
{
  ABE_L5_SETUP = 0, ABE_L5_POINT_LOOP, ABE_L5_SEARCH, ABE_L5_INTERP,
  ABE_L5_AR_SHELLF, ABE_L5_AR_WEIGHT, ABE_L5_POST, ABE_L5_CLEANUP,
  ABE_L5_TOTAL
};

#endif


Patch::Patch(int DIM, int *shapei, double *bboxi, int levi, bool buflog, int Symmetry) : lev(levi)
{

  int hbuffer_width = buffer_width;
  if (lev == 0)
    hbuffer_width = CS_width; // specific for shell-box coulping

  if (DIM != dim)
  {
    cout << "dimension is not consistent in Patch construction" << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  for (int i = 0; i < dim; i++)
  {
    shape[i] = shapei[i];
    bbox[i] = bboxi[i];
    bbox[dim + i] = bboxi[dim + i];
    lli[i] = uui[i] = 0;
    if (buflog)
    {
      double DH;
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
      DH = (bbox[dim + i] - bbox[i]) / (shape[i] - 1);
#else
#ifdef Cell
      DH = (bbox[dim + i] - bbox[i]) / shape[i];
#else
#error Not define Vertex nor Cell
#endif
#endif
      uui[i] = hbuffer_width;
      bbox[dim + i] = bbox[dim + i] + uui[i] * DH;
      shape[i] = shape[i] + uui[i];
    }
  }

  if (buflog)
  {
    if (DIM != 3)
    {
      cout << "Symmetry in Patch construction only support 3 yet but dim = " << DIM << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    double tmpb, DH;
    if (Symmetry > 0)
    {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
      DH = (bbox[5] - bbox[2]) / (shape[2] - 1);
#else
#ifdef Cell
      DH = (bbox[5] - bbox[2]) / shape[2];
#else
#error Not define Vertex nor Cell
#endif
#endif
      tmpb = Mymax(0, bbox[2] - hbuffer_width * DH);
      lli[2] = int((bbox[2] - tmpb) / DH + 0.4);
      bbox[2] = bbox[2] - lli[2] * DH;
      shape[2] = shape[2] + lli[2];
      if (lli[2] < hbuffer_width)
      {
        if (feq(bbox[2], 0, DH / 2))
          lli[2] = 0;
        else
        {
          cout << "Code mistake for lli[2] = " << lli[2] << ", bbox[2] = " << bbox[2] << endl;
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
      }
      if (Symmetry > 1)
      {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
        DH = (bbox[3] - bbox[0]) / (shape[0] - 1);
#else
#ifdef Cell
        DH = (bbox[3] - bbox[0]) / shape[0];
#else
#error Not define Vertex nor Cell
#endif
#endif
        tmpb = Mymax(0, bbox[0] - hbuffer_width * DH);
        lli[0] = int((bbox[0] - tmpb) / DH + 0.4);
        bbox[0] = bbox[0] - lli[0] * DH;
        shape[0] = shape[0] + lli[0];
        if (lli[0] < hbuffer_width)
        {
          if (feq(bbox[0], 0, DH / 2))
            lli[0] = 0;
          else
          {
            cout << "Code mistake for lli[0] = " << lli[0] << ", bbox[0] = " << bbox[0] << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
          }
        }
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
        DH = (bbox[4] - bbox[1]) / (shape[1] - 1);
#else
#ifdef Cell
        DH = (bbox[4] - bbox[1]) / shape[1];
#else
#error Not define Vertex nor Cell
#endif
#endif
        tmpb = Mymax(0, bbox[1] - hbuffer_width * DH);
        lli[1] = int((bbox[1] - tmpb) / DH + 0.4);
        bbox[1] = bbox[1] - lli[1] * DH;
        shape[1] = shape[1] + lli[1];
        if (lli[1] < hbuffer_width)
        {
          if (feq(bbox[1], 0, DH / 2))
            lli[1] = 0;
          else
          {
            cout << "Code mistake for lli[1] = " << lli[1] << ", bbox[1] = " << bbox[1] << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
          }
        }
      }
      else
      {
        for (int i = 0; i < 2; i++)
        {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
          DH = (bbox[dim + i] - bbox[i]) / (shape[i] - 1);
#else
#ifdef Cell
          DH = (bbox[dim + i] - bbox[i]) / shape[i];
#else
#error Not define Vertex nor Cell
#endif
#endif
          lli[i] = hbuffer_width;
          bbox[i] = bbox[i] - lli[i] * DH;
          shape[i] = shape[i] + lli[i];
        }
      }
    }
    else
    {
      for (int i = 0; i < dim; i++)
      {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
        DH = (bbox[dim + i] - bbox[i]) / (shape[i] - 1);
#else
#ifdef Cell
        DH = (bbox[dim + i] - bbox[i]) / shape[i];
#else
#error Not define Vertex nor Cell
#endif
#endif
        lli[i] = hbuffer_width;
        bbox[i] = bbox[i] - lli[i] * DH;
        shape[i] = shape[i] + lli[i];
      }
    }
  }

  blb = ble = 0;
}
Patch::~Patch()
{
}
// buflog 1: with buffer points; 0 without
void Patch::checkPatch(bool buflog)
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == 0)
  {
    if (buflog)
    {
      cout << " belong to level " << lev << endl;
      cout << " shape: [";
      for (int i = 0; i < dim; i++)
      {
        cout << shape[i];
        if (i < dim - 1)
          cout << ",";
        else
          cout << "]";
      }
      cout << " resolution: [";
      for (int i = 0; i < dim; i++)
      {
        cout << getdX(i);
        if (i < dim - 1)
          cout << ",";
        else
          cout << "]" << endl;
      }
      cout << " range:" << "(";
      for (int i = 0; i < dim; i++)
      {
        cout << bbox[i] << ":" << bbox[dim + i];
        if (i < dim - 1)
          cout << ",";
        else
          cout << ")" << endl;
      }
    }
    else
    {
      cout << " belong to level " << lev << endl;
      cout << " shape: [";
      for (int i = 0; i < dim; i++)
      {
        cout << shape[i] - lli[i] - uui[i];
        if (i < dim - 1)
          cout << ",";
        else
          cout << "]";
      }
      cout << " resolution: [";
      for (int i = 0; i < dim; i++)
      {
        cout << getdX(i);
        if (i < dim - 1)
          cout << ",";
        else
          cout << "]" << endl;
      }
      cout << " range:" << "(";
      for (int i = 0; i < dim; i++)
      {
        cout << bbox[i] + lli[i] * getdX(i) << ":" << bbox[dim + i] - uui[i] * getdX(i);
        if (i < dim - 1)
          cout << ",";
        else
          cout << ")" << endl;
      }
    }
  }
}
void Patch::checkPatch(bool buflog, const int out_rank)
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == out_rank)
  {
    cout << " out_rank = " << out_rank << endl;
    if (buflog)
    {
      cout << " belong to level " << lev << endl;
      cout << " shape: [";
      for (int i = 0; i < dim; i++)
      {
        cout << shape[i];
        if (i < dim - 1)
          cout << ",";
        else
          cout << "]";
      }
      cout << " resolution: [";
      for (int i = 0; i < dim; i++)
      {
        cout << getdX(i);
        if (i < dim - 1)
          cout << ",";
        else
          cout << "]" << endl;
      }
      cout << " range:" << "(";
      for (int i = 0; i < dim; i++)
      {
        cout << bbox[i] << ":" << bbox[dim + i];
        if (i < dim - 1)
          cout << ",";
        else
          cout << ")" << endl;
      }
    }
    else
    {
      cout << " belong to level " << lev << endl;
      cout << " shape: [";
      for (int i = 0; i < dim; i++)
      {
        cout << shape[i] - lli[i] - uui[i];
        if (i < dim - 1)
          cout << ",";
        else
          cout << "]";
      }
      cout << " resolution: [";
      for (int i = 0; i < dim; i++)
      {
        cout << getdX(i);
        if (i < dim - 1)
          cout << ",";
        else
          cout << "]" << endl;
      }
      cout << " range:" << "(";
      for (int i = 0; i < dim; i++)
      {
        cout << bbox[i] + lli[i] * getdX(i) << ":" << bbox[dim + i] - uui[i] * getdX(i);
        if (i < dim - 1)
          cout << ",";
        else
          cout << ")" << endl;
      }
    }
  }
}
void Patch::Interp_Points(MyList<var> *VarList,
                          int NN, double **XX,
                          double *Shellf, int Symmetry)
{
  const int wave_prof_on = v7_wave_interp_consume();
#if ENABLE_V7_WAVE_PROFILING
  const double wave_interp_t0 = wave_prof_on ? MPI_Wtime() : 0.0;
  double wave_interp_comm = 0.0;
#endif
  // NOTE: we do not Synchnize variables here, make sure of that before calling this routine
  // Level-5 phase timers for this Interp_Points() call (profiling only).  The
  // per-rank mark was armed by surf_MassPAng; unmarked callers skip all timing.
  const int l5_on = abe_l5_interp_consume();
  // Level-6: arm f_global_interp() phase capture in fmisc.f90 while this
  // marked Interp_Points() call runs (profiling only; unmarked callers
  // never reach this point).
  if (l5_on)
  {
    abe_l6_gi_arm(1);
    abe_l7_p3_arm(1);
  }
  double l5_t0 = 0.0, l5_t_setup = 0.0, l5_t_ploop = 0.0;
  double l5_interp = 0.0, l5_ar1a = 0.0, l5_ar1b = 0.0;
  double l5_ar2a = 0.0, l5_ar2b = 0.0, l5_t_post = 0.0, l5_t_clean = 0.0;
  long long l5_nloc = 0;
  if (l5_on)
    l5_t0 = ABE_PROFILE_WTIME();
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  int ordn = 2 * ghost_width;
  MyList<var> *varl;
  int num_var = 0;
  varl = VarList;
  while (varl)
  {
    num_var++;
    varl = varl->next;
  }

  double *shellf;
  shellf = new double[NN * num_var];
  memset(shellf, 0, sizeof(double) * NN * num_var);

  // we use weight to monitor code, later some day we can move it for optimization
  int *weight;
  weight = new int[NN];
  memset(weight, 0, sizeof(int) * NN);

  double *DH, *llb, *uub;
  DH = new double[dim];

  for (int i = 0; i < dim; i++)
  {
    DH[i] = getdX(i);
  }
  llb = new double[dim];
  uub = new double[dim];

  if (l5_on)
    l5_t_setup = ABE_PROFILE_WTIME();
  for (int j = 0; j < NN; j++) // run along points
  {
    double pox[dim];
    for (int i = 0; i < dim; i++)
    {
      pox[i] = XX[i][j];
      if (myrank == 0 && (XX[i][j] < bbox[i] + lli[i] * DH[i] || XX[i][j] > bbox[dim + i] - uui[i] * DH[i]))
      {
        cout << "Patch::Interp_Points: point (";
        for (int k = 0; k < dim; k++)
        {
          cout << XX[k][j];
          if (k < dim - 1)
            cout << ",";
          else
            cout << ") is out of current Patch." << endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }

    MyList<Block> *Bp = blb;
    bool notfind = true;
    while (notfind && Bp) // run along Blocks
    {
      Block *BP = Bp->data;

      bool flag = true;
      for (int i = 0; i < dim; i++)
      {
// NOTE: our dividing structure is (exclude ghost)
// -1 0
//       1  2
// so (0,1) does not belong to any part for vertex structure
// here we put (0,0.5) to left part and (0.5,1) to right part
// BUT for cell structure the bbox is (-1.5,0.5) and (0.5,2.5), there is no missing region at all
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
        llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + (ghost_width - 0.5) * DH[i];
        uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - (ghost_width - 0.5) * DH[i];
#else
#ifdef Cell
        llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + ghost_width * DH[i];
        uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - ghost_width * DH[i];
#else
#error Not define Vertex nor Cell
#endif
#endif
        if (XX[i][j] - llb[i] < -DH[i] / 2 || XX[i][j] - uub[i] > DH[i] / 2)
        {
          flag = false;
          break;
        }
      }

      if (flag)
      {
        notfind = false;
        if (myrank == BP->rank)
        {
          //---> interpolation
          varl = VarList;
          int k = 0;
          double l5_ti = 0.0;
          if (l5_on)
            l5_ti = ABE_PROFILE_WTIME();
          while (varl) // run along variables
          {
            //              shellf[j*num_var+k] = Parallel::global_interp(dim,BP->shape,BP->X,BP->fgfs[varl->data->sgfn],
            //	  		                                    pox,ordn,varl->data->SoA,Symmetry);
            f_global_interp(BP->shape, BP->X[0], BP->X[1], BP->X[2], BP->fgfs[varl->data->sgfn], shellf[j * num_var + k],
                            pox[0], pox[1], pox[2], ordn, varl->data->SoA, Symmetry);
            varl = varl->next;
            k++;
          }
          weight[j] = 1;
          if (l5_on)
          {
            l5_interp += ABE_PROFILE_WTIME() - l5_ti;
            l5_nloc++;
          }
        }
      }
      if (Bp == ble)
        break;
      Bp = Bp->next;
    }
  }

  if (l5_on)
  {
    l5_t_ploop = ABE_PROFILE_WTIME();
    l5_ar1a = ABE_PROFILE_WTIME();
  }
#if ENABLE_V7_WAVE_PROFILING
  const double wave_comm_t0 = wave_prof_on ? MPI_Wtime() : 0.0;
#endif
  MPI_Allreduce(shellf, Shellf, NN * num_var, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#if ENABLE_V7_WAVE_PROFILING
  if (wave_prof_on) wave_interp_comm += MPI_Wtime() - wave_comm_t0;
#endif
  if (l5_on)
    l5_ar1b = ABE_PROFILE_WTIME();
  int *Weight;
  Weight = new int[NN];
  if (l5_on)
    l5_ar2a = ABE_PROFILE_WTIME();
#if ENABLE_V7_WAVE_PROFILING
  const double wave_comm_t1 = wave_prof_on ? MPI_Wtime() : 0.0;
#endif
  MPI_Allreduce(weight, Weight, NN, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#if ENABLE_V7_WAVE_PROFILING
  if (wave_prof_on) wave_interp_comm += MPI_Wtime() - wave_comm_t1;
#endif
  if (l5_on)
  {
    l5_ar2b = ABE_PROFILE_WTIME();
    l5_t_post = ABE_PROFILE_WTIME();
  }

  //  misc::tillherecheck("print me");

  for (int i = 0; i < NN; i++)
  {
    if (Weight[i] > 1)
    {
      if (myrank == 0)
        cout << "WARNING: Patch::Interp_Points meets multiple weight" << endl;
      for (int j = 0; j < num_var; j++)
        Shellf[j + i * num_var] = Shellf[j + i * num_var] / Weight[i];
    }
    else if (Weight[i] == 0 && myrank == 0)
    {
      cout << "ERROR: Patch::Interp_Points fails to find point (";
      for (int j = 0; j < dim; j++)
      {
        cout << XX[j][i];
        if (j < dim - 1)
          cout << ",";
        else
          cout << ")";
      }
      cout << " on Patch (";
      for (int j = 0; j < dim; j++)
      {
        cout << bbox[j] << "+" << lli[j] * getdX(j);
        if (j < dim - 1)
          cout << ",";
        else
          cout << ")--";
      }
      cout << "(";
      for (int j = 0; j < dim; j++)
      {
        cout << bbox[dim + j] << "-" << uui[j] * getdX(j);
        if (j < dim - 1)
          cout << ",";
        else
          cout << ")" << endl;
      }
#if 0
       checkBlock();
#else
      cout << "splited domains:" << endl;
      {
        MyList<Block> *Bp = blb;
        while (Bp)
        {
          Block *BP = Bp->data;

          for (int i = 0; i < dim; i++)
          {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
            llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + (ghost_width - 0.5) * DH[i];
            uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - (ghost_width - 0.5) * DH[i];
#else
#ifdef Cell
            llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + ghost_width * DH[i];
            uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - ghost_width * DH[i];
#else
#error Not define Vertex nor Cell
#endif
#endif
          }
          cout << "(";
          for (int j = 0; j < dim; j++)
          {
            cout << llb[j] << ":" << uub[j];
            if (j < dim - 1)
              cout << ",";
            else
              cout << ")" << endl;
          }
          if (Bp == ble)
            break;
          Bp = Bp->next;
        }
      }
#endif
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  if (l5_on)
    l5_t_clean = ABE_PROFILE_WTIME();
  delete[] shellf;
  delete[] weight;
  delete[] Weight;
  delete[] DH;
  delete[] llb;
  delete[] uub;
#if ENABLE_V7_WAVE_PROFILING
  if (wave_prof_on) {
    const double elapsed = MPI_Wtime() - wave_interp_t0;
    v7_wave_record(V7_WAVE_INTERP_COMPUTE, elapsed - wave_interp_comm);
    v7_wave_record(V7_WAVE_INTERP_COMM, wave_interp_comm);
  }
#endif
  if (l5_on)
  {
    const double l5_t_end = ABE_PROFILE_WTIME();
    const double l5_dt_total = l5_t_end - l5_t0;
    const double l5_dt_ploop = l5_t_ploop - l5_t_setup;
    const double l5_dt_search = (l5_dt_ploop - l5_interp > 0.0) ? l5_dt_ploop - l5_interp : 0.0;
    abe_l5s_rec(ABE_L5_TOTAL, l5_dt_total);
    abe_l5s_rec(ABE_L5_SETUP, l5_t_setup - l5_t0);
    abe_l5s_rec(ABE_L5_POINT_LOOP, l5_dt_ploop);
    abe_l5s_rec(ABE_L5_SEARCH, l5_dt_search);
    abe_l5s_rec(ABE_L5_INTERP, l5_interp);
    abe_l5s_rec(ABE_L5_AR_SHELLF, l5_ar1b - l5_ar1a);
    abe_l5s_rec(ABE_L5_AR_WEIGHT, l5_ar2b - l5_ar2a);
    abe_l5s_rec(ABE_L5_POST, l5_t_clean - l5_t_post);
    abe_l5s_rec(ABE_L5_CLEANUP, l5_t_end - l5_t_clean);
    abe_l5_pts_total += (long long)NN;
    abe_l5_pts_local += l5_nloc;
    // Level-6/7: disarm f_global_interp()/polin3() phase capture in fmisc.f90
    // (this marked Interp_Points() call ended).
    abe_l6_gi_arm(0);
    abe_l7_p3_arm(0);
  }
}
void Patch::Interp_Points(MyList<var> *VarList,
                          int NN, double **XX,
                          double *Shellf, int Symmetry, MPI_Comm Comm_here)
{
  // NOTE: we do not Synchnize variables here, make sure of that before calling this routine
  // Level-5 phase timers for this Interp_Points() call (profiling only).  The
  // per-rank mark was armed by surf_MassPAng; unmarked callers skip all timing.
  const int l5_on = abe_l5_interp_consume();
  // Level-6: arm f_global_interp() phase capture in fmisc.f90 while this
  // marked Interp_Points() call runs (profiling only; unmarked callers
  // never reach this point).
  if (l5_on)
  {
    abe_l6_gi_arm(1);
    abe_l7_p3_arm(1);
  }
  double l5_t0 = 0.0, l5_t_setup = 0.0, l5_t_ploop = 0.0;
  double l5_interp = 0.0, l5_ar1a = 0.0, l5_ar1b = 0.0;
  double l5_ar2a = 0.0, l5_ar2b = 0.0, l5_t_post = 0.0, l5_t_clean = 0.0;
  long long l5_nloc = 0;
  if (l5_on)
    l5_t0 = ABE_PROFILE_WTIME();
  int myrank, lmyrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_rank(Comm_here, &lmyrank);

  int ordn = 2 * ghost_width;
  MyList<var> *varl;
  int num_var = 0;
  varl = VarList;
  while (varl)
  {
    num_var++;
    varl = varl->next;
  }

  double *shellf;
  shellf = new double[NN * num_var];
  memset(shellf, 0, sizeof(double) * NN * num_var);

  // we use weight to monitor code, later some day we can move it for optimization
  int *weight;
  weight = new int[NN];
  memset(weight, 0, sizeof(int) * NN);

  double *DH, *llb, *uub;
  DH = new double[dim];

  for (int i = 0; i < dim; i++)
  {
    DH[i] = getdX(i);
  }
  llb = new double[dim];
  uub = new double[dim];

  if (l5_on)
    l5_t_setup = ABE_PROFILE_WTIME();
  for (int j = 0; j < NN; j++) // run along points
  {
    double pox[dim];
    for (int i = 0; i < dim; i++)
    {
      pox[i] = XX[i][j];
      if (lmyrank == 0 && (XX[i][j] < bbox[i] + lli[i] * DH[i] || XX[i][j] > bbox[dim + i] - uui[i] * DH[i]))
      {
        cout << "Patch::Interp_Points: point (";
        for (int k = 0; k < dim; k++)
        {
          cout << XX[k][j];
          if (k < dim - 1)
            cout << ",";
          else
            cout << ") is out of current Patch." << endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }

    MyList<Block> *Bp = blb;
    bool notfind = true;
    while (notfind && Bp) // run along Blocks
    {
      Block *BP = Bp->data;

      bool flag = true;
      for (int i = 0; i < dim; i++)
      {
// NOTE: our dividing structure is (exclude ghost)
// -1 0
//       1  2
// so (0,1) does not belong to any part for vertex structure
// here we put (0,0.5) to left part and (0.5,1) to right part
// BUT for cell structure the bbox is (-1.5,0.5) and (0.5,2.5), there is no missing region at all
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
        llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + (ghost_width - 0.5) * DH[i];
        uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - (ghost_width - 0.5) * DH[i];
#else
#ifdef Cell
        llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + ghost_width * DH[i];
        uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - ghost_width * DH[i];
#else
#error Not define Vertex nor Cell
#endif
#endif
        if (XX[i][j] - llb[i] < -DH[i] / 2 || XX[i][j] - uub[i] > DH[i] / 2)
        {
          flag = false;
          break;
        }
      }

      if (flag)
      {
        notfind = false;
        if (myrank == BP->rank)
        {
          //---> interpolation
          varl = VarList;
          int k = 0;
          double l5_ti = 0.0;
          if (l5_on)
            l5_ti = ABE_PROFILE_WTIME();
          while (varl) // run along variables
          {
            //              shellf[j*num_var+k] = Parallel::global_interp(dim,BP->shape,BP->X,BP->fgfs[varl->data->sgfn],
            //	  		                                    pox,ordn,varl->data->SoA,Symmetry);
            f_global_interp(BP->shape, BP->X[0], BP->X[1], BP->X[2], BP->fgfs[varl->data->sgfn], shellf[j * num_var + k],
                            pox[0], pox[1], pox[2], ordn, varl->data->SoA, Symmetry);
            varl = varl->next;
            k++;
          }
          weight[j] = 1;
          if (l5_on)
          {
            l5_interp += ABE_PROFILE_WTIME() - l5_ti;
            l5_nloc++;
          }
        }
      }
      if (Bp == ble)
        break;
      Bp = Bp->next;
    }
  }

  if (l5_on)
  {
    l5_t_ploop = ABE_PROFILE_WTIME();
    l5_ar1a = ABE_PROFILE_WTIME();
  }
  MPI_Allreduce(shellf, Shellf, NN * num_var, MPI_DOUBLE, MPI_SUM, Comm_here);
  if (l5_on)
    l5_ar1b = ABE_PROFILE_WTIME();
  int *Weight;
  Weight = new int[NN];
  if (l5_on)
    l5_ar2a = ABE_PROFILE_WTIME();
  MPI_Allreduce(weight, Weight, NN, MPI_INT, MPI_SUM, Comm_here);
  if (l5_on)
  {
    l5_ar2b = ABE_PROFILE_WTIME();
    l5_t_post = ABE_PROFILE_WTIME();
  }

  //  misc::tillherecheck("print me");
  //  if(lmyrank == 0) cout<<"myrank = "<<myrank<<"print me"<<endl;

  for (int i = 0; i < NN; i++)
  {
    if (Weight[i] > 1)
    {
      if (lmyrank == 0)
        cout << "WARNING: Patch::Interp_Points meets multiple weight" << endl;
      for (int j = 0; j < num_var; j++)
        Shellf[j + i * num_var] = Shellf[j + i * num_var] / Weight[i];
    }
#if 0 // for not involved levels, this may fail     
     else if(Weight[i] == 0 && lmyrank == 0)
     {
       cout<<"ERROR: Patch::Interp_Points fails to find point (";
       for(int j=0;j<dim;j++)
       {
	  cout<<XX[j][i];
	  if(j<dim-1) cout<<",";
	  else        cout<<")";
       }
       cout<<" on Patch (";
       for(int j=0;j<dim;j++)
       {
	  cout<<bbox[j]<<"+"<<lli[j]*getdX(j);
	  if(j<dim-1) cout<<",";
	  else        cout<<")--";
       }
       cout<<"(";
       for(int j=0;j<dim;j++)
       {
	  cout<<bbox[dim+j]<<"-"<<uui[j]*getdX(j);
	  if(j<dim-1) cout<<",";
	  else        cout<<")"<<endl;
       }
#if 0
       checkBlock();
#else
  cout<<"splited domains:"<<endl;
  {
     MyList<Block> *Bp=blb;
     while(Bp)
     {
	Block *BP=Bp->data;

	for(int i=0;i<dim;i++)
	{
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
          llb[i] = (feq(BP->bbox[i]    ,bbox[i]    ,DH[i]/2)) ? BP->bbox[i]+lli[i]*DH[i]     : BP->bbox[i]    +(ghost_width-0.5)*DH[i];
          uub[i] = (feq(BP->bbox[dim+i],bbox[dim+i],DH[i]/2)) ? BP->bbox[dim+i]-uui[i]*DH[i] : BP->bbox[dim+i]-(ghost_width-0.5)*DH[i];
#else
#ifdef Cell
          llb[i] = (feq(BP->bbox[i]    ,bbox[i]    ,DH[i]/2)) ? BP->bbox[i]+lli[i]*DH[i]     : BP->bbox[i]    +ghost_width*DH[i];
          uub[i] = (feq(BP->bbox[dim+i],bbox[dim+i],DH[i]/2)) ? BP->bbox[dim+i]-uui[i]*DH[i] : BP->bbox[dim+i]-ghost_width*DH[i];
#else
#error Not define Vertex nor Cell
#endif
#endif 
	}       
       cout<<"(";
       for(int j=0;j<dim;j++)
       {
	  cout<<llb[j]<<":"<<uub[j];
	  if(j<dim-1) cout<<",";
	  else        cout<<")"<<endl;
       }
	if(Bp == ble) break;
	Bp=Bp->next;
     }
  }
#endif       
       MPI_Abort(MPI_COMM_WORLD,1);
     }
#endif
  }

  if (l5_on)
    l5_t_clean = ABE_PROFILE_WTIME();
  delete[] shellf;
  delete[] weight;
  delete[] Weight;
  delete[] DH;
  delete[] llb;
  delete[] uub;
  if (l5_on)
  {
    const double l5_t_end = ABE_PROFILE_WTIME();
    const double l5_dt_total = l5_t_end - l5_t0;
    const double l5_dt_ploop = l5_t_ploop - l5_t_setup;
    const double l5_dt_search = (l5_dt_ploop - l5_interp > 0.0) ? l5_dt_ploop - l5_interp : 0.0;
    abe_l5s_rec(ABE_L5_TOTAL, l5_dt_total);
    abe_l5s_rec(ABE_L5_SETUP, l5_t_setup - l5_t0);
    abe_l5s_rec(ABE_L5_POINT_LOOP, l5_dt_ploop);
    abe_l5s_rec(ABE_L5_SEARCH, l5_dt_search);
    abe_l5s_rec(ABE_L5_INTERP, l5_interp);
    abe_l5s_rec(ABE_L5_AR_SHELLF, l5_ar1b - l5_ar1a);
    abe_l5s_rec(ABE_L5_AR_WEIGHT, l5_ar2b - l5_ar2a);
    abe_l5s_rec(ABE_L5_POST, l5_t_clean - l5_t_post);
    abe_l5s_rec(ABE_L5_CLEANUP, l5_t_end - l5_t_clean);
    abe_l5_pts_total += (long long)NN;
    abe_l5_pts_local += l5_nloc;
    // Level-6/7: disarm f_global_interp()/polin3() phase capture in fmisc.f90
    // (this marked Interp_Points() call ended).
    abe_l6_gi_arm(0);
    abe_l7_p3_arm(0);
  }
}
void Patch::checkBlock()
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == 0)
  {
    MyList<Block> *BP = blb;
    while (BP)
    {
      BP->data->checkBlock();
      if (BP == ble)
        break;
      BP = BP->next;
    }
  }
}

void Patch::Interp_Points_ReduceScatter(MyList<var> *VarList,
                                        int NN, double **XX,
                                        double *Shellf_local,
                                        int *Weight_local,
                                        int Symmetry, MPI_Comm Comm_here,
                                        int local_start, int local_count)
{
  const int l5_on = abe_l5_interp_consume();
  if (l5_on)
  {
    abe_l6_gi_arm(1);
    abe_l7_p3_arm(1);
  }
  const double t0 = ABE_PROFILE_WTIME();
  int myrank, lmyrank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_rank(Comm_here, &lmyrank);
  MPI_Comm_size(Comm_here, &comm_size);

  if (local_start < 0 || local_count < 0 || local_start + local_count > NN)
    MPI_Abort(MPI_COMM_WORLD, 1);

  int num_var = 0;
  for (MyList<var> *varl = VarList; varl; varl = varl->next)
    ++num_var;

  int ordn = 2 * ghost_width;
  double *local_shellf = new double[NN * num_var];
  int *local_weight = new int[NN];
  memset(local_shellf, 0, sizeof(double) * NN * num_var);
  memset(local_weight, 0, sizeof(int) * NN);

  double *DH = new double[dim];
  double *llb = new double[dim];
  double *uub = new double[dim];
  for (int i = 0; i < dim; ++i)
    DH[i] = getdX(i);

  const double t_setup = ABE_PROFILE_WTIME();
  for (int j = 0; j < NN; ++j)
  {
    double pox[dim];
    for (int i = 0; i < dim; ++i)
    {
      pox[i] = XX[i][j];
      if (lmyrank == 0 &&
          (XX[i][j] < bbox[i] + lli[i] * DH[i] ||
           XX[i][j] > bbox[dim + i] - uui[i] * DH[i]))
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MyList<Block> *Bp = blb;
    bool notfind = true;
    while (notfind && Bp)
    {
      Block *BP = Bp->data;
      bool flag = true;
      for (int i = 0; i < dim; ++i)
      {
#ifdef Vertex
        llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2))
                   ? BP->bbox[i] + lli[i] * DH[i]
                   : BP->bbox[i] + (ghost_width - 0.5) * DH[i];
        uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2))
                   ? BP->bbox[dim + i] - uui[i] * DH[i]
                   : BP->bbox[dim + i] - (ghost_width - 0.5) * DH[i];
#else
        llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2))
                   ? BP->bbox[i] + lli[i] * DH[i]
                   : BP->bbox[i] + ghost_width * DH[i];
        uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2))
                   ? BP->bbox[dim + i] - uui[i] * DH[i]
                   : BP->bbox[dim + i] - ghost_width * DH[i];
#endif
        if (XX[i][j] - llb[i] < -DH[i] / 2 ||
            XX[i][j] - uub[i] > DH[i] / 2)
        {
          flag = false;
          break;
        }
      }

      if (flag)
      {
        notfind = false;
        /* BP->rank is a MPI_COMM_WORLD rank, even for Comm_here. */
        if (myrank == BP->rank)
        {
          MyList<var> *varl = VarList;
          int k = 0;
          while (varl)
          {
            f_global_interp(BP->shape, BP->X[0], BP->X[1], BP->X[2],
                            BP->fgfs[varl->data->sgfn],
                            local_shellf[j * num_var + k],
                            pox[0], pox[1], pox[2], ordn,
                            varl->data->SoA, Symmetry);
            varl = varl->next;
            ++k;
          }
          local_weight[j] = 1;
        }
      }
      if (Bp == ble)
        break;
      Bp = Bp->next;
    }
  }

  int *recvcounts = new int[comm_size];
  int *weight_recvcounts = new int[comm_size];
  const int mp = NN / comm_size;
  const int Lp = NN - comm_size * mp;
  const int expected_start = (lmyrank < Lp)
                               ? lmyrank * (mp + 1)
                               : Lp * (mp + 1) + (lmyrank - Lp) * mp;
  const int expected_count = mp + (lmyrank < Lp ? 1 : 0);
  if (local_start != expected_start || local_count != expected_count)
    MPI_Abort(MPI_COMM_WORLD, 1);
  long long sum_shellf = 0;
  long long sum_weight = 0;
  for (int r = 0; r < comm_size; ++r)
  {
    const int point_count = mp + (r < Lp ? 1 : 0);
    recvcounts[r] = point_count * num_var;
    weight_recvcounts[r] = point_count;
    sum_shellf += recvcounts[r];
    sum_weight += weight_recvcounts[r];
  }
  if (sum_shellf != (long long)NN * num_var || sum_weight != NN)
    MPI_Abort(MPI_COMM_WORLD, 1);

  const double t_ar_shellf = ABE_PROFILE_WTIME();
  MPI_Reduce_scatter(local_shellf, Shellf_local, recvcounts,
                     MPI_DOUBLE, MPI_SUM, Comm_here);
  const double t_ar_weight = ABE_PROFILE_WTIME();
  MPI_Reduce_scatter(local_weight, Weight_local, weight_recvcounts,
                     MPI_INT, MPI_SUM, Comm_here);
  const double t_post = ABE_PROFILE_WTIME();

  for (int i = 0; i < local_count; ++i)
  {
    if (Weight_local[i] > 1)
    {
      if (lmyrank == 0)
        cout << "WARNING: Patch::Interp_Points_ReduceScatter meets multiple weight" << endl;
      for (int j = 0; j < num_var; ++j)
        Shellf_local[i * num_var + j] /= Weight_local[i];
    }
    else if (Weight_local[i] == 0 && lmyrank == 0)
    {
      cout << "ERROR: Patch::Interp_Points_ReduceScatter fails to find point "
           << (local_start + i) << endl;
    }
  }

  delete[] recvcounts;
  delete[] weight_recvcounts;
  delete[] local_shellf;
  delete[] local_weight;
  delete[] DH;
  delete[] llb;
  delete[] uub;

  if (l5_on)
  {
    abe_l5s_rec(ABE_L5_TOTAL, ABE_PROFILE_WTIME() - t0);
    abe_l5s_rec(ABE_L5_SETUP, t_setup - t0);
    abe_l5s_rec(ABE_L5_POINT_LOOP, t_ar_shellf - t_setup);
    abe_l5s_rec(ABE_L5_INTERP, t_ar_shellf - t_setup);
    abe_l5s_rec(ABE_L5_AR_SHELLF, t_ar_weight - t_ar_shellf);
    abe_l5s_rec(ABE_L5_AR_WEIGHT, t_post - t_ar_weight);
    abe_l5s_rec(ABE_L5_POST, ABE_PROFILE_WTIME() - t_post);
    abe_l5_pts_total += (long long)NN;
    abe_l5_pts_local += (long long)local_count;
    abe_l6_gi_arm(0);
    abe_l7_p3_arm(0);
  }
}

/* Sparse owner-computes uses the same first-match Block predicate as the
 * Interp_Points overloads above.  Block metadata is replicated on all ranks;
 * only the selected Block owner has X/fgfs data. */
static Block *patch_sparse_first_matching_block(Patch *patch,
                                                const double *point,
                                                const double *DH,
                                                double *llb, double *uub)
{
  MyList<Block> *Bp = patch->blb;
  while (Bp)
  {
    Block *BP = Bp->data;
    bool flag = true;
    for (int i = 0; i < dim; ++i)
    {
#ifdef Vertex
      llb[i] = (feq(BP->bbox[i], patch->bbox[i], DH[i] / 2))
                   ? BP->bbox[i] + patch->lli[i] * DH[i]
                   : BP->bbox[i] + (ghost_width - 0.5) * DH[i];
      uub[i] = (feq(BP->bbox[dim + i], patch->bbox[dim + i], DH[i] / 2))
                   ? BP->bbox[dim + i] - patch->uui[i] * DH[i]
                   : BP->bbox[dim + i] - (ghost_width - 0.5) * DH[i];
#else
      llb[i] = (feq(BP->bbox[i], patch->bbox[i], DH[i] / 2))
                   ? BP->bbox[i] + patch->lli[i] * DH[i]
                   : BP->bbox[i] + ghost_width * DH[i];
      uub[i] = (feq(BP->bbox[dim + i], patch->bbox[dim + i], DH[i] / 2))
                   ? BP->bbox[dim + i] - patch->uui[i] * DH[i]
                   : BP->bbox[dim + i] - ghost_width * DH[i];
#endif
      if (point[i] - llb[i] < -DH[i] / 2 || point[i] - uub[i] > DH[i] / 2)
      {
        flag = false;
        break;
      }
    }
    if (flag)
      return BP; // first matching Block wins, exactly as Interp_Points
    if (Bp == patch->ble)
      break;
    Bp = Bp->next;
  }
  return 0;
}

void Patch::Interp_Points_SparseOwner(MyList<var> *VarList,
                                      int NN, double **XX,
                                      double *Shellf_local,
                                      int *Weight_local,
                                      int Symmetry, MPI_Comm Comm_here,
                                      int local_start, int local_count)
{
  const int l5_on = abe_l5_interp_consume();
  if (l5_on)
  {
    abe_l6_gi_arm(1);
    abe_l7_p3_arm(1);
  }
  const double t0 = ABE_PROFILE_WTIME();

  int world_rank, local_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_rank(Comm_here, &local_rank);
  MPI_Comm_size(Comm_here, &comm_size);
  if (NN < 0 || local_start < 0 || local_count < 0 ||
      local_start > NN || local_count > NN - local_start || comm_size <= 0)
    MPI_Abort(MPI_COMM_WORLD, 1);
  const int mp = NN / comm_size;
  const int Lp = NN - comm_size * mp;
  const int expected_start = (local_rank < Lp)
                               ? local_rank * (mp + 1)
                               : Lp * (mp + 1) + (local_rank - Lp) * mp;
  const int expected_count = mp + (local_rank < Lp ? 1 : 0);
  if (local_start != expected_start || local_count != expected_count)
    MPI_Abort(MPI_COMM_WORLD, 1);

  int num_var = 0;
  for (MyList<var> *varl = VarList; varl; varl = varl->next)
    ++num_var;
  if (num_var <= 0 || (NN > 0 && num_var > INT_MAX / NN))
    MPI_Abort(MPI_COMM_WORLD, 1);

  for (int i = 0; i < local_count; ++i)
  {
    Weight_local[i] = 0;
    for (int k = 0; k < num_var; ++k)
      Shellf_local[i * num_var + k] = 0.0;
  }

  std::vector<int> comm_world_ranks(comm_size);
  MPI_Allgather(&world_rank, 1, MPI_INT, comm_world_ranks.data(), 1,
                MPI_INT, Comm_here);
  if (comm_world_ranks[local_rank] != world_rank)
    MPI_Abort(MPI_COMM_WORLD, 1);

  std::vector<double> DH(dim), llb(dim), uub(dim);
  for (int i = 0; i < dim; ++i)
    DH[i] = getdX(i);

  std::vector<int> owner_local(local_count, -1);
  std::vector<int> send_counts(comm_size, 0);
  int unresolved_consumer = 0;
  for (int local_n = 0; local_n < local_count; ++local_n)
  {
    const int n = local_start + local_n;
    double point[dim];
    bool patch_ok = true;
    for (int i = 0; i < dim; ++i)
    {
      point[i] = XX[i][n];
      if (point[i] < bbox[i] + lli[i] * DH[i] ||
          point[i] > bbox[dim + i] - uui[i] * DH[i])
        patch_ok = false;
    }
    Block *BP = patch_ok
                    ? patch_sparse_first_matching_block(this, point, DH.data(),
                                                        llb.data(), uub.data())
                    : 0;
    if (!BP)
    {
      ++unresolved_consumer;
      continue;
    }

    int destination = -1;
    for (int r = 0; r < comm_size; ++r)
      if (comm_world_ranks[r] == BP->rank)
      {
        destination = r;
        break;
      }
    if (destination < 0)
    {
      ++unresolved_consumer;
      continue;
    }
    owner_local[local_n] = destination;
    ++send_counts[destination];
  }

  std::vector<int> send_displs(comm_size, 0), recv_counts(comm_size, 0),
      recv_displs(comm_size, 0);
  for (int r = 1; r < comm_size; ++r)
    send_displs[r] = send_displs[r - 1] + send_counts[r - 1];
  const int total_sent = send_displs[comm_size - 1] + send_counts[comm_size - 1];
  if (total_sent > local_count || total_sent + unresolved_consumer != local_count)
    MPI_Abort(MPI_COMM_WORLD, 1);

  std::vector<int> send_indices(total_sent), send_local_slots(total_sent);
  std::vector<int> send_cursor = send_displs;
  for (int local_n = 0; local_n < local_count; ++local_n)
  {
    const int destination = owner_local[local_n];
    if (destination < 0)
      continue;
    const int pos = send_cursor[destination]++;
    send_indices[pos] = local_start + local_n;
    send_local_slots[pos] = local_n;
  }

  const double t_setup = ABE_PROFILE_WTIME();
  MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1,
               MPI_INT, Comm_here);
  for (int r = 1; r < comm_size; ++r)
    recv_displs[r] = recv_displs[r - 1] + recv_counts[r - 1];
  const int total_received = recv_displs[comm_size - 1] + recv_counts[comm_size - 1];
  if (total_received < 0 || total_received > INT_MAX / num_var ||
      total_sent > INT_MAX / num_var)
    MPI_Abort(MPI_COMM_WORLD, 1);

  std::vector<int> recv_indices(total_received);
  MPI_Alltoallv(send_indices.data(), send_counts.data(), send_displs.data(), MPI_INT,
                recv_indices.data(), recv_counts.data(), recv_displs.data(), MPI_INT,
                Comm_here);
  const double t_request_done = ABE_PROFILE_WTIME();

  std::vector<double> owner_results((size_t)total_received * num_var, 0.0);
  std::vector<int> owner_weights(total_received, 0);
  int owner_mismatch = 0;
  double interpolation_time = 0.0;
  int ordn = 2 * ghost_width;
  for (int j = 0; j < total_received; ++j)
  {
    const int n = recv_indices[j];
    if (n < 0 || n >= NN)
    {
      ++owner_mismatch;
      continue;
    }
    double point[dim];
    for (int i = 0; i < dim; ++i)
      point[i] = XX[i][n];
    Block *BP = patch_sparse_first_matching_block(this, point, DH.data(),
                                                  llb.data(), uub.data());
    if (!BP || BP->rank != world_rank)
    {
      ++owner_mismatch;
      continue;
    }
    const double ti = l5_on ? ABE_PROFILE_WTIME() : 0.0;
    bool batch_done = false;
#ifdef Cell
    if (num_var == 17 && ordn == 6)
    {
      double *fields[17];
      double soa[17][3];
      double values[17];
      MyList<var> *varl = VarList;
      for (int k = 0; k < 17; ++k)
      {
        if (!varl)
          MPI_Abort(MPI_COMM_WORLD, 1);
        fields[k] = BP->fgfs[varl->data->sgfn];
        for (int d = 0; d < 3; ++d)
          soa[k][d] = varl->data->SoA[d];
        varl = varl->next;
      }
      int batch_status = 1;
      f_global_interp_batch17(BP->shape, BP->X[0], BP->X[1], BP->X[2],
                              fields[0], fields[1], fields[2], fields[3],
                              fields[4], fields[5], fields[6], fields[7],
                              fields[8], fields[9], fields[10], fields[11],
                              fields[12], fields[13], fields[14], fields[15],
                              fields[16], values,
                              point[0], point[1], point[2], ordn,
                              &soa[0][0], Symmetry, batch_status);
      if (batch_status == 0)
      {
#if VERIFY_BATCH17
        // Correctness-only dual evaluation.  At most twelve representative
        // owner points are recorded per rank over the whole run: the first,
        // middle, and last received point plus the first point seen for each
        // reflection-mask class.  This block is absent from normal builds.
        static bool verify_seen_mask[8] = {false, false, false, false,
                                           false, false, false, false};
        static bool verify_seen_position[3] = {false, false, false};
        static int verify_sample_count = 0;
        int verify_mask = 0;
        for (int d = 0; d < 3; ++d)
        {
          const double delta = BP->X[d][1] - BP->X[d][0];
          int center = (int)((point[d] - BP->X[d][0]) / delta + 0.4) + 1;
          int lower = center - ordn / 2 + 1;
          int minimum = 1;
          if ((Symmetry == 2 && d < 2 && fabs(BP->X[d][0]) < delta) ||
              (Symmetry != 0 && d == 2 && fabs(BP->X[d][0]) < delta))
            minimum = -ordn / 2 + 1;
          if (lower < minimum)
            lower = minimum;
          int upper = lower + ordn - 1;
          if (upper > BP->shape[d])
            lower = BP->shape[d] + 1 - ordn;
          if (lower <= 0)
            verify_mask |= 1 << d;
        }
        int position_class = -1;
        if (j == 0)
          position_class = 0;
        else if (j == total_received / 2)
          position_class = 1;
        else if (j == total_received - 1)
          position_class = 2;
        const bool positional = position_class >= 0 &&
                                !verify_seen_position[position_class];
        const bool select_sample = verify_sample_count < 12 &&
                                   (positional || !verify_seen_mask[verify_mask]);
        if (select_sample)
        {
          verify_seen_mask[verify_mask] = true;
          if (position_class >= 0)
            verify_seen_position[position_class] = true;
          double reference[17];
          for (int k = 0; k < 17; ++k)
            f_global_interp(BP->shape, BP->X[0], BP->X[1], BP->X[2],
                            fields[k], reference[k],
                            point[0], point[1], point[2], ordn,
                            &soa[k][0], Symmetry);

          int block_ordinal = 0;
          for (MyList<Block> *scan = blb; scan; scan = scan->next)
          {
            if (scan->data == BP)
              break;
            ++block_ordinal;
            if (scan == ble)
              break;
          }
          char verify_name[80];
          sprintf(verify_name, "batch17_verify_rank%06d.dat", world_rank);
          ios_base::openmode mode = ios::out;
          if (verify_sample_count > 0)
            mode |= ios::app;
          ofstream verify_out(verify_name, mode);
          if (verify_sample_count == 0)
            verify_out << "# global_index x y z block_ordinal owner_rank "
                          "reflection_mask var_index var_name soa_x soa_y soa_z "
                          "reference batch abs_error rel_error\n";
          static const char *verify_var_names[17] = {
              "Sfx_rhs", "Sfy_rhs", "Sfz_rhs", "chi", "trK",
              "gxx", "gxy", "gxz", "gyy", "gyz", "gzz",
              "Axx", "Axy", "Axz", "Ayy", "Ayz", "Azz"};
          verify_out << setprecision(17);
          for (int k = 0; k < 17; ++k)
          {
            const double abs_error = fabs(values[k] - reference[k]);
            const double rel_error = abs_error / max(fabs(reference[k]), 1.0e-300);
            verify_out << n << ' ' << point[0] << ' ' << point[1] << ' '
                       << point[2] << ' ' << block_ordinal << ' ' << world_rank
                       << ' ' << verify_mask << ' ' << k << ' '
                       << verify_var_names[k] << ' ' << soa[k][0] << ' '
                       << soa[k][1] << ' ' << soa[k][2] << ' '
                       << reference[k] << ' ' << values[k] << ' '
                       << abs_error << ' ' << rel_error << '\n';
          }
          ++verify_sample_count;
        }
#endif
        for (int k = 0; k < 17; ++k)
          owner_results[j * num_var + k] = values[k];
        batch_done = true;
      }
    }
#endif
    if (!batch_done)
    {
      MyList<var> *varl = VarList;
      int k = 0;
      while (varl)
      {
        f_global_interp(BP->shape, BP->X[0], BP->X[1], BP->X[2],
                        BP->fgfs[varl->data->sgfn],
                        owner_results[j * num_var + k],
                        point[0], point[1], point[2], ordn,
                        varl->data->SoA, Symmetry);
        varl = varl->next;
        ++k;
      }
    }
    if (l5_on)
      interpolation_time += ABE_PROFILE_WTIME() - ti;
    owner_weights[j] = 1;
  }
  if (owner_mismatch)
    cerr << "ERROR: Patch::Interp_Points_SparseOwner owner validation failed on rank "
         << world_rank << " for " << owner_mismatch << " point(s)" << endl;

  std::vector<int> result_send_counts(comm_size), result_send_displs(comm_size),
      result_recv_counts(comm_size), result_recv_displs(comm_size);
  for (int r = 0; r < comm_size; ++r)
  {
    result_send_counts[r] = recv_counts[r] * num_var;
    result_send_displs[r] = recv_displs[r] * num_var;
    result_recv_counts[r] = send_counts[r] * num_var;
    result_recv_displs[r] = send_displs[r] * num_var;
  }
  std::vector<double> consumer_results((size_t)total_sent * num_var, 0.0);
  MPI_Alltoallv(owner_results.data(), result_send_counts.data(),
                result_send_displs.data(), MPI_DOUBLE,
                consumer_results.data(), result_recv_counts.data(),
                result_recv_displs.data(), MPI_DOUBLE, Comm_here);
  const double t_result_done = ABE_PROFILE_WTIME();

  std::vector<int> consumer_weights(total_sent, 0);
  MPI_Alltoallv(owner_weights.data(), recv_counts.data(), recv_displs.data(), MPI_INT,
                consumer_weights.data(), send_counts.data(), send_displs.data(), MPI_INT,
                Comm_here);
  const double t_weight_done = ABE_PROFILE_WTIME();

  for (int pos = 0; pos < total_sent; ++pos)
  {
    const int local_n = send_local_slots[pos];
    Weight_local[local_n] = consumer_weights[pos];
    for (int k = 0; k < num_var; ++k)
      Shellf_local[local_n * num_var + k] = consumer_results[pos * num_var + k];
  }

  for (int local_n = 0; local_n < local_count; ++local_n)
  {
    if (Weight_local[local_n] > 1)
    {
      cerr << "WARNING: Patch::Interp_Points_SparseOwner meets multiple weight" << endl;
      for (int k = 0; k < num_var; ++k)
        Shellf_local[local_n * num_var + k] /= Weight_local[local_n];
    }
    else if (Weight_local[local_n] == 0)
    {
      cerr << "ERROR: Patch::Interp_Points_SparseOwner fails to find point "
           << (local_start + local_n) << endl;
    }
  }
  if (l5_on)
  {
    const double tend = ABE_PROFILE_WTIME();
    const double point_loop = t_request_done - t_setup;
    abe_l5s_rec(ABE_L5_TOTAL, tend - t0);
    abe_l5s_rec(ABE_L5_SETUP, t_setup - t0);
    abe_l5s_rec(ABE_L5_POINT_LOOP, point_loop);
    abe_l5s_rec(ABE_L5_SEARCH, point_loop);
    abe_l5s_rec(ABE_L5_INTERP, interpolation_time);
    abe_l5s_rec(ABE_L5_AR_SHELLF, t_result_done - t_request_done);
    abe_l5s_rec(ABE_L5_AR_WEIGHT, t_weight_done - t_result_done);
    abe_l5s_rec(ABE_L5_POST, tend - t_weight_done);
    abe_l5_pts_total += (long long)NN;
    abe_l5_pts_local += (long long)total_received;
    abe_l6_gi_arm(0);
    abe_l7_p3_arm(0);
  }
}

 double Patch::getdX(int dir)
{
  if (dir < 0 || dir >= dim)
  {
    cout << "Patch::getdX: error input dir = " << dir << ", this Patch has direction (0," << dim - 1 << ")" << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  double h;
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
  if (shape[dir] == 1)
  {
    cout << "Patch::getdX: for direction " << dir << ", this Patch has only one point. Can not determine dX for vertex center grid." << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  h = (bbox[dim + dir] - bbox[dir]) / (shape[dir] - 1);
#else
#ifdef Cell
  h = (bbox[dim + dir] - bbox[dir]) / shape[dir];
#else
#error Not define Vertex nor Cell
#endif
#endif
  return h;
}
bool Patch::Interp_ONE_Point(MyList<var> *VarList, double *XX,
                             double *Shellf, int Symmetry)
{
  // NOTE: we do not Synchnize variables here, make sure of that before calling this routine
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  int ordn = 2 * ghost_width;
  MyList<var> *varl;
  int num_var = 0;
  varl = VarList;
  while (varl)
  {
    num_var++;
    varl = varl->next;
  }

  double *shellf;
  shellf = new double[num_var];
  memset(shellf, 0, sizeof(double) * num_var);

  double *DH, *llb, *uub;
  DH = new double[dim];

  for (int i = 0; i < dim; i++)
  {
    DH[i] = getdX(i);
  }
  llb = new double[dim];
  uub = new double[dim];

  double pox[dim];
  for (int i = 0; i < dim; i++)
  {
    pox[i] = XX[i];
    // has excluded the buffer points
    if (XX[i] < bbox[i] + lli[i] * DH[i] - DH[i] / 100 || XX[i] > bbox[dim + i] - uui[i] * DH[i] + DH[i] / 100)
    {
      delete[] shellf;
      delete[] DH;
      delete[] llb;
      delete[] uub;
      return false; // out of current patch,
                    // remember to delete the allocated arrays before return!!!
    }
  }

  MyList<Block> *Bp = blb;
  bool notfind = true;
  while (notfind && Bp) // run along Blocks
  {
    Block *BP = Bp->data;

    bool flag = true;
    for (int i = 0; i < dim; i++)
    {
// NOTE: our dividing structure is (exclude ghost)
// -1 0
//       1  2
// so (0,1) does not belong to any part for vertex structure
// here we put (0,0.5) to left part and (0.5,1) to right part
// BUT for cell structure the bbox is (-1.5,0.5) and (0.5,2.5), there is no missing region at all
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
      llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + (ghost_width - 0.5) * DH[i];
      uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - (ghost_width - 0.5) * DH[i];
#else
#ifdef Cell
      llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + ghost_width * DH[i];
      uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - ghost_width * DH[i];
#else
#error Not define Vertex nor Cell
#endif
#endif
      if (XX[i] - llb[i] < -DH[i] / 2 || XX[i] - uub[i] > DH[i] / 2)
      {
        flag = false;
        break;
      }
    }

    if (flag)
    {
      notfind = false;
      if (myrank == BP->rank)
      {
// test old code
#if 0
#define floorint(a) ((a) < 0 ? int(a) - 1 : int(a))
//---> interpolation
                int ixl,iyl,izl,ixu,iyu,izu;
	    	double Delx,Dely,Delz;

		ixl = 1+floorint((pox[0]-BP->X[0][0])/DH[0]);
	   	iyl = 1+floorint((pox[1]-BP->X[1][0])/DH[1]);
	   	izl = 1+floorint((pox[2]-BP->X[2][0])/DH[2]);

		int nn=ordn/2;

		ixl = ixl-nn;
		iyl = iyl-nn;
		izl = izl-nn;
	   
		int tmi;
		tmi = (Symmetry==2)?-1:0;
		if(ixl<tmi) ixl=tmi;
	   	if(iyl<tmi) iyl=tmi;
		tmi = (Symmetry>0)?-1:0;
	   	if(izl<tmi) izl=tmi;
      
	   	if(ixl+ordn>BP->shape[0]) ixl=BP->shape[0]-ordn;
	   	if(iyl+ordn>BP->shape[1]) iyl=BP->shape[1]-ordn;
	   	if(izl+ordn>BP->shape[2]) izl=BP->shape[2]-ordn;
// support cell center
		if(ixl>=0) Delx = ( pox[0] - BP->X[0][ixl] )/ DH[0];
		else       Delx = ( pox[0] + BP->X[0][0] )/ DH[0];
                if(iyl>=0) Dely = ( pox[1] - BP->X[1][iyl] )/ DH[1];
		else       Dely = ( pox[1] + BP->X[1][0] )/ DH[1];
                if(izl>=0) Delz = ( pox[2] - BP->X[2][izl] )/ DH[2];
		else       Delz = ( pox[2] + BP->X[2][0] )/ DH[2];
//change to fortran index
                ixl++;
	   	iyl++;
	   	izl++;
	   	ixu = ixl + ordn - 1;
	   	iyu = iyl + ordn - 1;
	   	izu = izl + ordn - 1;
	    	varl=VarList;
		int j=0;
	    	while(varl)
		{
                 f_interp_2(BP->shape,BP->fgfs[varl->data->sgfn],shellf[j],ixl,ixu,iyl,iyu,izl,izu,Delx,Dely,Delz,
                                     ordn,varl->data->SoA,Symmetry);
		 varl=varl->next;
		 j++;
		} //varl
#else
        //---> interpolation
        varl = VarList;
        int k = 0;
        while (varl) // run along variables
        {
          //              shellf[j*num_var+k] = Parallel::global_interp(dim,BP->shape,BP->X,BP->fgfs[varl->data->sgfn],
          //	  		                                    pox,ordn,varl->data->SoA,Symmetry);
          f_global_interp(BP->shape, BP->X[0], BP->X[1], BP->X[2], BP->fgfs[varl->data->sgfn], shellf[k],
                          pox[0], pox[1], pox[2], ordn, varl->data->SoA, Symmetry);
          varl = varl->next;
          k++;
        }
#endif
      }
    }
    if (Bp == ble)
      break;
    Bp = Bp->next;
  }

  if (notfind && myrank == 0)
  {
    cout << "ERROR: Patch::Interp_Points fails to find point (";
    for (int j = 0; j < dim; j++)
    {
      cout << XX[j];
      if (j < dim - 1)
        cout << ",";
      else
        cout << ")";
    }
    cout << " on Patch (";
    for (int j = 0; j < dim; j++)
    {
      cout << bbox[j] << "+" << lli[j] * getdX(j);
      if (j < dim - 1)
        cout << ",";
      else
        cout << ")--";
    }
    cout << "(";
    for (int j = 0; j < dim; j++)
    {
      cout << bbox[dim + j] << "-" << uui[j] * getdX(j);
      if (j < dim - 1)
        cout << ",";
      else
        cout << ")" << endl;
    }
#if 0
       checkBlock();
#else
    cout << "splited domains:" << endl;
    {
      MyList<Block> *Bp = blb;
      while (Bp)
      {
        Block *BP = Bp->data;

        for (int i = 0; i < dim; i++)
        {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
          llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + (ghost_width - 0.5) * DH[i];
          uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - (ghost_width - 0.5) * DH[i];
#else
#ifdef Cell
          llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + ghost_width * DH[i];
          uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - ghost_width * DH[i];
#else
#error Not define Vertex nor Cell
#endif
#endif
        }
        cout << "(";
        for (int j = 0; j < dim; j++)
        {
          cout << llb[j] << ":" << uub[j];
          if (j < dim - 1)
            cout << ",";
          else
            cout << ")" << endl;
        }
        if (Bp == ble)
          break;
        Bp = Bp->next;
      }
    }
#endif
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  MPI_Allreduce(shellf, Shellf, num_var, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  delete[] shellf;
  delete[] DH;
  delete[] llb;
  delete[] uub;

  return true;
}
bool Patch::Interp_ONE_Point(MyList<var> *VarList, double *XX,
                             double *Shellf, int Symmetry, MPI_Comm Comm_here)
{
  // NOTE: we do not Synchnize variables here, make sure of that before calling this routine
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  int ordn = 2 * ghost_width;
  MyList<var> *varl;
  int num_var = 0;
  varl = VarList;
  while (varl)
  {
    num_var++;
    varl = varl->next;
  }

  double *shellf;
  shellf = new double[num_var];
  memset(shellf, 0, sizeof(double) * num_var);

  double *DH, *llb, *uub;
  DH = new double[dim];

  for (int i = 0; i < dim; i++)
  {
    DH[i] = getdX(i);
  }
  llb = new double[dim];
  uub = new double[dim];

  double pox[dim];
  for (int i = 0; i < dim; i++)
  {
    pox[i] = XX[i];
    // has excluded the buffer points
    if (XX[i] < bbox[i] + lli[i] * DH[i] - DH[i] / 100 || XX[i] > bbox[dim + i] - uui[i] * DH[i] + DH[i] / 100)
    {
      delete[] shellf;
      delete[] DH;
      delete[] llb;
      delete[] uub;
      return false; // out of current patch,
                    // remember to delete the allocated arrays before return!!!
    }
  }

  MyList<Block> *Bp = blb;
  bool notfind = true;
  while (notfind && Bp) // run along Blocks
  {
    Block *BP = Bp->data;

    bool flag = true;
    for (int i = 0; i < dim; i++)
    {
// NOTE: our dividing structure is (exclude ghost)
// -1 0
//       1  2
// so (0,1) does not belong to any part for vertex structure
// here we put (0,0.5) to left part and (0.5,1) to right part
// BUT for cell structure the bbox is (-1.5,0.5) and (0.5,2.5), there is no missing region at all
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
      llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + (ghost_width - 0.5) * DH[i];
      uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - (ghost_width - 0.5) * DH[i];
#else
#ifdef Cell
      llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + ghost_width * DH[i];
      uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - ghost_width * DH[i];
#else
#error Not define Vertex nor Cell
#endif
#endif
      if (XX[i] - llb[i] < -DH[i] / 2 || XX[i] - uub[i] > DH[i] / 2)
      {
        flag = false;
        break;
      }
    }

    if (flag)
    {
      notfind = false;
      if (myrank == BP->rank)
      {
// test old code
#if 0
#define floorint(a) ((a) < 0 ? int(a) - 1 : int(a))
//---> interpolation
                int ixl,iyl,izl,ixu,iyu,izu;
	    	double Delx,Dely,Delz;

		ixl = 1+floorint((pox[0]-BP->X[0][0])/DH[0]);
	   	iyl = 1+floorint((pox[1]-BP->X[1][0])/DH[1]);
	   	izl = 1+floorint((pox[2]-BP->X[2][0])/DH[2]);

		int nn=ordn/2;

		ixl = ixl-nn;
		iyl = iyl-nn;
		izl = izl-nn;
	   
		int tmi;
		tmi = (Symmetry==2)?-1:0;
		if(ixl<tmi) ixl=tmi;
	   	if(iyl<tmi) iyl=tmi;
		tmi = (Symmetry>0)?-1:0;
	   	if(izl<tmi) izl=tmi;
      
	   	if(ixl+ordn>BP->shape[0]) ixl=BP->shape[0]-ordn;
	   	if(iyl+ordn>BP->shape[1]) iyl=BP->shape[1]-ordn;
	   	if(izl+ordn>BP->shape[2]) izl=BP->shape[2]-ordn;
// support cell center
		if(ixl>=0) Delx = ( pox[0] - BP->X[0][ixl] )/ DH[0];
		else       Delx = ( pox[0] + BP->X[0][0] )/ DH[0];
                if(iyl>=0) Dely = ( pox[1] - BP->X[1][iyl] )/ DH[1];
		else       Dely = ( pox[1] + BP->X[1][0] )/ DH[1];
                if(izl>=0) Delz = ( pox[2] - BP->X[2][izl] )/ DH[2];
		else       Delz = ( pox[2] + BP->X[2][0] )/ DH[2];
//change to fortran index
                ixl++;
	   	iyl++;
	   	izl++;
	   	ixu = ixl + ordn - 1;
	   	iyu = iyl + ordn - 1;
	   	izu = izl + ordn - 1;
	    	varl=VarList;
		int j=0;
	    	while(varl)
		{
                 f_interp_2(BP->shape,BP->fgfs[varl->data->sgfn],shellf[j],ixl,ixu,iyl,iyu,izl,izu,Delx,Dely,Delz,
                                     ordn,varl->data->SoA,Symmetry);
		 varl=varl->next;
		 j++;
		} //varl
#else
        //---> interpolation
        varl = VarList;
        int k = 0;
        while (varl) // run along variables
        {
          //              shellf[j*num_var+k] = Parallel::global_interp(dim,BP->shape,BP->X,BP->fgfs[varl->data->sgfn],
          //	  		                                    pox,ordn,varl->data->SoA,Symmetry);
          f_global_interp(BP->shape, BP->X[0], BP->X[1], BP->X[2], BP->fgfs[varl->data->sgfn], shellf[k],
                          pox[0], pox[1], pox[2], ordn, varl->data->SoA, Symmetry);
          varl = varl->next;
          k++;
        }
#endif
      }
    }
    if (Bp == ble)
      break;
    Bp = Bp->next;
  }

  if (notfind && myrank == 0)
  {
    cout << "ERROR: Patch::Interp_Points fails to find point (";
    for (int j = 0; j < dim; j++)
    {
      cout << XX[j];
      if (j < dim - 1)
        cout << ",";
      else
        cout << ")";
    }
    cout << " on Patch (";
    for (int j = 0; j < dim; j++)
    {
      cout << bbox[j] << "+" << lli[j] * getdX(j);
      if (j < dim - 1)
        cout << ",";
      else
        cout << ")--";
    }
    cout << "(";
    for (int j = 0; j < dim; j++)
    {
      cout << bbox[dim + j] << "-" << uui[j] * getdX(j);
      if (j < dim - 1)
        cout << ",";
      else
        cout << ")" << endl;
    }
#if 0
       checkBlock();
#else
    cout << "splited domains:" << endl;
    {
      MyList<Block> *Bp = blb;
      while (Bp)
      {
        Block *BP = Bp->data;

        for (int i = 0; i < dim; i++)
        {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
          llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + (ghost_width - 0.5) * DH[i];
          uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - (ghost_width - 0.5) * DH[i];
#else
#ifdef Cell
          llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + ghost_width * DH[i];
          uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - ghost_width * DH[i];
#else
#error Not define Vertex nor Cell
#endif
#endif
        }
        cout << "(";
        for (int j = 0; j < dim; j++)
        {
          cout << llb[j] << ":" << uub[j];
          if (j < dim - 1)
            cout << ",";
          else
            cout << ")" << endl;
        }
        if (Bp == ble)
          break;
        Bp = Bp->next;
      }
    }
#endif
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  MPI_Allreduce(shellf, Shellf, num_var, MPI_DOUBLE, MPI_SUM, Comm_here);

  delete[] shellf;
  delete[] DH;
  delete[] llb;
  delete[] uub;

  return true;
}
// find maximum of abstract value, XX store position for maximum, Shellf store maximum themselvs
void Patch::Find_Maximum(MyList<var> *VarList, double *XX,
                         double *Shellf)
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  MyList<var> *varl;
  int num_var = 0;
  varl = VarList;
  while (varl)
  {
    num_var++;
    varl = varl->next;
  }

  double *shellf, *xx;
  shellf = new double[num_var];
  xx = new double[dim * num_var];
  memset(shellf, 0, sizeof(double) * num_var);
  memset(xx, 0, sizeof(double) * dim * num_var);

  double *DH;
  int *llb, *uub;
  DH = new double[dim];

  for (int i = 0; i < dim; i++)
  {
    DH[i] = getdX(i);
  }

  llb = new int[dim];
  uub = new int[dim];

  MyList<Block> *Bp = blb;
  while (Bp) // run along Blocks
  {
    Block *BP = Bp->data;

    if (myrank == BP->rank)
    {

      for (int i = 0; i < dim; i++)
      {
        llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? lli[i] : ghost_width;
        uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? uui[i] : ghost_width;
      }

      varl = VarList;
      int k = 0;
      double tmp, tmpx[dim];
      while (varl) // run along variables
      {
        f_find_maximum(BP->shape, BP->X[0], BP->X[1], BP->X[2], BP->fgfs[varl->data->sgfn], tmp, tmpx, llb, uub);
        if (tmp > shellf[k])
        {
          shellf[k] = tmp;
          for (int i = 0; i < dim; i++)
            xx[dim * k + i] = tmpx[i];
        }
        varl = varl->next;
        k++;
      }
    }

    if (Bp == ble)
      break;
    Bp = Bp->next;
  }

  struct mloc
  {
    double val;
    int rank;
  };

  mloc *IN, *OUT;
  IN = new mloc[num_var];
  OUT = new mloc[num_var];
  for (int i = 0; i < num_var; i++)
  {
    IN[i].val = shellf[i];
    IN[i].rank = myrank;
  }

  MPI_Allreduce(IN, OUT, num_var, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

  for (int i = 0; i < num_var; i++)
  {
    Shellf[i] = OUT[i].val;
    if (myrank != OUT[i].rank)
      for (int k = 0; k < 3; k++)
        xx[3 * i + k] = 0;
  }

  MPI_Allreduce(xx, XX, dim * num_var, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  delete[] IN;
  delete[] OUT;
  delete[] shellf;
  delete[] xx;
  delete[] DH;
  delete[] llb;
  delete[] uub;
}
void Patch::Find_Maximum(MyList<var> *VarList, double *XX,
                         double *Shellf, MPI_Comm Comm_here)
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  MyList<var> *varl;
  int num_var = 0;
  varl = VarList;
  while (varl)
  {
    num_var++;
    varl = varl->next;
  }

  double *shellf, *xx;
  shellf = new double[num_var];
  xx = new double[dim * num_var];
  memset(shellf, 0, sizeof(double) * num_var);
  memset(xx, 0, sizeof(double) * dim * num_var);

  double *DH;
  int *llb, *uub;
  DH = new double[dim];

  for (int i = 0; i < dim; i++)
  {
    DH[i] = getdX(i);
  }

  llb = new int[dim];
  uub = new int[dim];

  MyList<Block> *Bp = blb;
  while (Bp) // run along Blocks
  {
    Block *BP = Bp->data;

    if (myrank == BP->rank)
    {

      for (int i = 0; i < dim; i++)
      {
        llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? lli[i] : ghost_width;
        uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? uui[i] : ghost_width;
      }

      varl = VarList;
      int k = 0;
      double tmp, tmpx[dim];
      while (varl) // run along variables
      {
        f_find_maximum(BP->shape, BP->X[0], BP->X[1], BP->X[2], BP->fgfs[varl->data->sgfn], tmp, tmpx, llb, uub);
        if (tmp > shellf[k])
        {
          shellf[k] = tmp;
          for (int i = 0; i < dim; i++)
            xx[dim * k + i] = tmpx[i];
        }
        varl = varl->next;
        k++;
      }
    }

    if (Bp == ble)
      break;
    Bp = Bp->next;
  }

  struct mloc
  {
    double val;
    int rank;
  };

  mloc *IN, *OUT;
  IN = new mloc[num_var];
  OUT = new mloc[num_var];
  for (int i = 0; i < num_var; i++)
  {
    IN[i].val = shellf[i];
    IN[i].rank = myrank;
  }

  MPI_Allreduce(IN, OUT, num_var, MPI_DOUBLE_INT, MPI_MAXLOC, Comm_here);

  for (int i = 0; i < num_var; i++)
  {
    Shellf[i] = OUT[i].val;
    if (myrank != OUT[i].rank)
      for (int k = 0; k < 3; k++)
        xx[3 * i + k] = 0;
  }

  MPI_Allreduce(xx, XX, dim * num_var, MPI_DOUBLE, MPI_SUM, Comm_here);

  delete[] IN;
  delete[] OUT;
  delete[] shellf;
  delete[] xx;
  delete[] DH;
  delete[] llb;
  delete[] uub;
}
// if the given point locates in the present Patch return true
// otherwise return false
bool Patch::Find_Point(double *XX)
{
  double *DH;
  DH = new double[dim];

  for (int i = 0; i < dim; i++)
  {
    DH[i] = getdX(i);
  }

  for (int i = 0; i < dim; i++)
  {
    // has excluded the buffer points
    if (XX[i] < bbox[i] + lli[i] * DH[i] - DH[i] / 100 || XX[i] > bbox[dim + i] - uui[i] * DH[i] + DH[i] / 100)
    {
      delete[] DH;
      return false; // out of current patch,
                    // remember to delete the allocated arrays before return!!!
    }
  }

  delete[] DH;

  return true;
}
