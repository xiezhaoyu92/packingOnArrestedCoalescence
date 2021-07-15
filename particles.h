//
//  particles.h
//  ellipsoidRelaxation
//
//  Created by Chris Burke on 2/28/16.
//  Copyright © 2016 tufts. All rights reserved.
//

#ifndef particles_h
#define particles_h

#include <stdio.h>



typedef struct {
    double position[3];
    double oldPosition[3];
    double preRelaxPosition[3];
    double postPreviousRelaxPosition[3];
    double rad;
    double forceDet[3]; // deterministic force
    double forceStoc[3]; // stochastic force
    int    coord[4];
} particle;

double closestDistance(int np, particle *p);
int overlapQ(particle *p1, particle *p2);
int anyOverlapQ(int nPart[], double maxRadius, particle *p, int pCount[], particle *pp[],  double pPosArray[] );

void undoOverlaps(double *aa[], double maxRadius, particle p[], int nPart[], int pCount[], particle *pp[], double pPosArray[], int maxUndoSteps, double *dtOverlap, int *nearJammed);
void undoOverlapsCompare(double *aa[], double maxRadius, particle q[] ,particle p[],  int nPart[], int pCount[], particle *pp[], double pPosArray[],
                         int maxUndoSteps, double *dtOverlap, int *nearJammed);
void undoOverlapsOld(int np, double a, double b, double rad, particle *p,
                     int maxUndoSteps, double *dtOverlap, int *nearJammed);
void undoOverlapsOldProfile(int np, double a, double b, double rad, particle *p, int maxUndoSteps, double *dtOverlap, int *nearJammed, FILE* undoFile);

void verletIntegrate(particle *p, double dt);
void gradDescentIntegrate(particle *p, double dt);
void gradDescentConstrainedIntegrate(double a, double b, particle *p, double dt);

double findTypicalSpacing(int np, double rad, particle *p, int *nPart, int *pCount, particle **pp, double *pPosArray);
double findTypicalSpacingOld(int np, double rad, particle *p);
int distancesToNeighbors(int nPart[], double maxRadius, particle *p, int pCount[], particle *pp[],  double pPosArray[], int *nNeighbors, double *distanceList);

int step(double a, double b, double delta, const double *xin, double *xout);double surfaceAreaAtEvoParameter(double a0, double b0, double s);

#endif /* particles_h */
