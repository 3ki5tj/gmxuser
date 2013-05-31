#ifndef GMX_UNUSED__
#define GMX_UNUSED__

#include <stdio.h>

#define CALL_UNUSED_MDRUN_H() \
  CALL_UNUSED_FORCE_H() \
  CALL_UNUSED_VEC_H()

#define CALL_UNUSED_PME_H() \
  CALL_UNUSED_GMXCOMPLEX_H()

#define CALL_UNUSED_STATUTIL_H() \
  CALL_UNUSED_READINP_H()

/* define macros that calls the unused variables
 * in various headers of GROMACS
 * since they are macros, there's no need to
 * include unnecessary headers in order to make
 * macros in this file sensible */

#define CALL_UNUSED_VEC_H() \
  if (gmx_numzero(1.0) == -1) {  \
    /* maths.h included in vec.h */ \
    rvec r1, r2, r3; \
    m_rvec_to_one(1, &r2); \
    m_rveccopy(1, &r2, &r1); \
    calc_lll(r1, r2); \
    m_rvecsub(1, &r1, &r2, &r3); \
    m_rvecadd(1, &r1, &r3, &r2); \
  }

#define CALL_UNUSED_FORCE_H() \
  if (sepdvdlformat[0] == '\0')  printf("%c", sepdvdlformat[0]);

#define CALL_UNUSED_GROMPP_H() \
  if (ds[0][0] == -1)  printf("%d", ds[0][0]);

#define CALL_UNUSED_GMXCOMPLEX_H() { \
  t_complex a; \
  a = cadd(cnul, rcexp(0.0f)); \
  a = cdiv(csub(cnul, a), a); }

#define CALL_UNUSED_READINP_H()  \
  if (argtp[0][0] == '\b') printf("%s", "");

#endif /* GMX_UNUSED__ */

