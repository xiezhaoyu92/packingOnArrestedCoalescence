//
//  misc.h
//  ellipsoidRelaxation
//
//  Created by Chris Burke on 2/28/16.
//  Copyright © 2016 tufts. All rights reserved.
//

#ifndef misc_h
#define misc_h

#include <stdio.h>
#include "particles.h"


#ifndef TRUE
#define TRUE -1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef PI
#define PI 3.141592653589793
#endif

#ifndef SQRT2
#define SQRT2 1.414213562
#endif

#ifndef SQRT3d3
#define SQRT3d3 0.57735026919
#endif

#ifndef MACHINE_EPSILON
#define MACHINE_EPSILON 1e-14
#endif



#ifndef DEBUG0
#define getIndex3d(a,b,c,a0) ((a)*a0[1]*a0[2]+(b)*a0[2]+(c))
#define getIndex3dShort(a,a0) (a[0]*a0[1]*a0[2] + a[1]*a0[2] + a[2])
#define getIndex4d(a,b,c,d,a0) ((a)*a0[1]*a0[2]*a0[3] + (b)*a0[2]*a0[3] + (c)*a0[3] + d)
#define getIndex4dShort(a,a0) (a[0]*a0[1]*a0[2]*a0[3] + a[1]*a0[2]*a0[3] + a[2]*a0[3] + a[3])

#define vectorNorm(v1) (sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]) )
#define vectorNormSquare(v1) (v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2])
#define vectorDotProduct(v1,v2) (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])
#endif

void vectorAdd(const double v1[3], const double v2[3], double *v3);
void vectorSubtract(double v1[3],double v2[3],double *v3);
void vectorScale(const double v1[3], double scale, double v2[3]);
void vectorCopy(const double *x1, double *x2);
#ifdef DEBUG0
double vectorNorm(double v1[3]);
double vectorNormSquare(double v1[3]);
double vectorDotProduct(double v1[3],double v2[3]);
#endif


//partition array alloc backgroun
double* arrayAlloc1d(double *x, int n);
int* arrayAlloc1dInt(int *x, int n);
particle** arrayAlloc1dParticlePointer(particle **x, int n);
double** arrayAlloc2d(double **x, int n1, int n2);

int compareDoubles(const void *a, const void *b);


void randinitialize(void);
int randInteger(int max);
double randNormal(void);




#endif /* misc_h */
