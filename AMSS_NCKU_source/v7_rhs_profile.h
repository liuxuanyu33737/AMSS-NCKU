#ifndef V7_RHS_PROFILE_H
#define V7_RHS_PROFILE_H

#ifndef ENABLE_V7_RHS_PROFILING
#define ENABLE_V7_RHS_PROFILING 0
#endif

#if defined(__cplusplus) && ENABLE_V7_RHS_PROFILING
extern "C" void v7_rhs_profile_report_();
#endif

#endif
