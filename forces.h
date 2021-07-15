//
//  forces.h
//  ellipsoidRelaxation
//
//  Created by Chris Burke on 2/28/16.
//  Copyright © 2016 tufts. All rights reserved.
//
#ifndef FORCES_H
#define FORCES_H



int addAttractiveForce(particle *p1, particle *p2, double attrStrength, double* newForce);
int addAttractiveLJForce(particle *p1, particle *p2, double attrStrength, double* newForce);
int addRepulsiveForce(particle *p1, particle *p2, double attrStrength, double newForce[3]);
double repulsiveEnergy(particle *p1, particle *p2, double repStrength);
int addScreenedAttractiveForce(particle *p1, particle *p2, double screeningLength, double attrStrength);
int addScreenedAttractiveLJForce(particle *p1, particle *p2, double screeningLength, double attrStrength);
#ifdef USE_GSL
double screenedAttractiveEnergy(particle *p1, particle *p2, double screeningLength, double attrStrength);
#endif
int addGravitationalForce(particle *p, double gravStrength);
int addStochasticForce(double a, double b, double D, double rad, particle *p);
int addGenericForce(particle *p, double f[3], double strength);
int addOverlapForce(particle *p1, particle *p2, double rad, double strength);
int addHertzianForce(particle *p1, particle *p2, double rad, double strength);
double hertzianEnergy(particle *p1, particle *p2, double rad, double strength);
int addHarmonicSurfaceForce(double a, double b, particle *p1, double strength);
double harmonicSurfaceEnergy(double a, double b, particle *p1, double strength);


int applyToPairs(double searchRadius, int centerCenterCutoff, double rad, double strength, particle p[],  int nPart[], int pCount[], particle *pp[], double pPosArray[],
                 int (*pairInteraction)(particle*, particle*, double, double));

double totalPairEnergy(double searchRadius, int centerCenterCutoff, double rad, double strength, particle p[], int nPart[], int pCount[], particle *pp[], double pPosArray[],
                       double (*pairEnergy)(particle*, particle*, double, double));

#endif