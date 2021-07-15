//
//  DropletCoalescence.h
//  ellipsoidRelaxation
//
//  Created by Zhaoyu Xie on 7/11/16.
//  Copyright © 2016 tufts. All rights reserved.
//

#ifndef DropletCoalescence_h
#define DropletCoalescence_h

#ifndef PI
#define PI 3.141592653589793
#endif

#include "particles.h"

double DCconstraintFunction(double z0, double a, double b, double *x);
void DCconstraintGradient(double z0, double a, double b, double *x, double *grad);

void DCRelaxationProject(double z0, double aOld, double bOld, double a, double b, particle *p);

int DCtangent1(double z0, double a, double b, const double *x, double *t1);
int DCtangent2(double z0, double a, double b, const double *x, double *t2);

double DCconformalFactor(double z0, double a, double b, double u);
double DCconformalInvert(double z0, double a, double b, double u2);


void SetupPartitions(double rad, double a, double b, double z0, int nPart[], double pPosArray[]);

void DCgradDescentConstrainedIntegrate(double z0, double a, double b, particle *p, double dt);
void DCundoOverlaps(double z0, double a, double b, double maxRadius, particle p[],  int nPart[], int pCount[], particle *pp[], double pPosArray[],
                    int maxUndoSteps, double *dtOverlap, int *nearJammed);
int DCstep(double z0, double a, double b, double delta, double *xin, double *xout);

void DCsetPartitionNumber(int np, double rad, double z0, double a, double b, int nPart[], double pFactor);
particle** DCallocatePP(int np, double rad, double z0, double a, double b, int nPart[], particle *pp[], double pFactor);
int* DCallocatePCount(double rad, double z0, double a, double b, int nPart[], int pCount[], double pFactor);
void DCmakePartitions(double z0, double a, double b, double maxRadius, int nPart[], double pPosArray[]);

double rhoSquare(double z0, double a, double b, double z);
double volumeDC (double z0, double a, double b);
void findNewParameters(double *a, double *b, double volume, double lambda);
double findA ( double b, double volume, double machep, double t);

#endif /* DropletCoalescence_h */
