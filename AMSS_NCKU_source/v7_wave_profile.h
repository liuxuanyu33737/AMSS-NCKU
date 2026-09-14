#ifndef V7_WAVE_PROFILE_H
#define V7_WAVE_PROFILE_H

#ifndef ENABLE_V7_WAVE_PROFILING
#define ENABLE_V7_WAVE_PROFILING 0
#endif

#if ENABLE_V7_WAVE_PROFILING
enum V7WavePhase {
  V7_WAVE_TOTAL = 0, V7_WAVE_INTERP_TOTAL, V7_WAVE_INTERP_COMPUTE,
  V7_WAVE_INTERP_COMM, V7_WAVE_ANGULAR_COEFF_INTEGRAL,
  V7_WAVE_OUTPUT_COMM, V7_WAVE_OTHER, V14_ANALYSIS_TOTAL, V14_MASS_TOTAL, V7_WAVE_NPHASE
};
class V14AnalysisScope {
  int phase_; double start_;
public:
  explicit V14AnalysisScope(int phase);
  ~V14AnalysisScope();
};
void v7_wave_record(int phase, double seconds);
void v7_wave_interp_mark(void);
int v7_wave_interp_consume(void);
void v7_wave_report(int myrank, int nprocs);
#else
#define v7_wave_record(phase, seconds) ((void)0)
#define v7_wave_interp_mark()          ((void)0)
#define v7_wave_interp_consume()       (0)
#define v7_wave_report(rank, nprocs)   ((void)0)
#endif

#endif
