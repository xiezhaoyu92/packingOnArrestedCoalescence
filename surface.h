//
//  surface.h
//  ellipsoidRelaxation
//
//  Created by Chris Burke on 2/28/16.
//  Copyright © 2016 tufts. All rights reserved.
//

#ifndef surface_h
#define surface_h

#include <stdio.h>


#ifndef MACHINE_EPSILON
#define MACHINE_EPSILON 1e-14
#endif

#ifndef DEBUG0
#define constraintFunction(a,b,x) ((x[0]/(a))*(x[0]/(a))+(x[1]/(a))*(x[1]/(a))+(x[2]/(b))*(x[2]/(b))-1)
#endif


//surface projection
int closestEllipsePoint(double a, double b, double rho, double z, double *rhoNew, double *zNew);
int closestCirclePoint(double a, double rho, double z, double *rhoNew, double *zNew);
int project(double a, double b, const double *xin, double *xout);
int moveWithSurface(double a1, double b1, double a2, double b2, const double *xin, double *xout);

int tangent1(double a, double b, const double *xin, double *t1) ;
int tangent2(double a, double b, const double *xin, double *t2);


double ellipsoidSurfaceArea(double a, double b);
double areaAtTime(double a0, double b0, double f, double t);
double timeAtArea(double a0, double b0, double f, double area);
double surfaceAreaAtEvoParameter(double a0, double b0, double s);
double evoParameterAtSurfaceArea(double a0, double b0, double area);
void ellipsoidUpdate(double t, double a0, double b0, double *a, double *b);

double conformalFactor(double a, double b, double u);
double conformalInvert(double a, double b, double u2);

#ifdef DEBUG0
double constraintFunction(double a, double b, double *x);
#endif
void constraintGradient(double a, double b, double *x, double *grad);
void vectorFromSurface(double a, double b, double *x, double *dx);



#endif /* surface_h */
