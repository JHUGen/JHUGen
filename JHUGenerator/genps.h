#ifndef genps_h
#define genps_h

#include <stdio.h>
#include <math.h>
#include <time.h>

extern void           gensing_(int *, double *, double *, double [][4], int *, int *, double *, int *);
extern void inline    swapmom_(double *, double *);
extern void           genps_( int *, double *, double *, double *, double [][4], double *);

static inline void    boostback(double *, double *, double *, double *);
static inline double  split(double *, double * , double *, double *, double *);
static void           checkOSrel(int *, double *, double *, double [][4], double);
static void           checkEMcon(int *, double *, double *, double [][4], double);
static inline double  kappa(double, double , double );
static inline double  sq(double);

void sincos(double, double*, double*);

#endif
