//
//  partition.h
//  ellipsoidRelaxation
//
//  Created by Chris Burke on 2/28/16.
//  Copyright © 2016 tufts. All rights reserved.
//
#ifndef PARTITION_H
#define PARTITION_H

#include "particles.h"


//partition setup
void setPartitionNumber(int np, double rad, double *aa[], int nPart[], double pFactor);
void getMaxDimensions(double *aa[], int nPart[], double max[]);
particle** allocatePP(int np, double rad, double *aa[], int nPart[], particle *pp[], double pFactor);
int* allocatePCount(double rad, double *aa[], int nPart[], int pCount[], double pFactor);
void makePartitions(double *aa[], double maxRadius, int nPart[], double pPosArray[] );

//partition main
void rePartitionParticle(particle *p, int nPart[], int pCount[], particle *pp[], double pPosArray[] );
void partitionParticles(particle p[], int nPart[], int pCount[], particle *pp[], double pPosArray[]);
void findPartition(double p[], int nPart[], double pPosArray[] ,int c[] );
int findPartitionPos(double position, int nPart, double pPos[]);
int greaterPartition(int dx, int dy, int dz);

//partition indexing
/*int getIndex3d(int a, int b, int c, int a0[4]);
 int getIndex4d(int a, int b, int c, int d, int a0[4]);
 int getIndex4dShort(int a[4], int a0[4]);
 int getIndex3dShort(int a[4], int a0[4]);*/

//partition debugging
int checkpCount(int np, int nPart[], int pCount[]);
int checkpCountP(int np, int nPart[], int pCount[], particle *p);
int checkpCountPP(int np, int nPart[], int pCount[], particle *pp[]);
int checkPvsPP(int np, int nPart[], int pCount[], particle *p, particle *pp[]);
int checkOverlapFull(particle *p, int np);
int runChecks(int np, int nPart[], int pCount[], particle *p, particle *pp[]);
int runChecksMinusOverlap(int np, int nPart[], int pCount[], particle *p, particle *pp[]);

#endif