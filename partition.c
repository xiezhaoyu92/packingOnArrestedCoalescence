//
//  partition.c
//  ellipsoidRelaxation
//
//  Created by Chris Burke on 2/28/16.
//  Copyright © 2016 tufts. All rights reserved.
//
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <stdio.h> //for debugging
#include <float.h>
#include <time.h>
#include "particles.h"
#include "partition.h"
#include "mt19937ar.h"
#include "misc.h"

/*********************
 * partition setup    *
 *********************/

/*
 * Input:   rad             - max particle radius
 *          aa              - length of axes
 *          nPart           - number of partitions
 *          pFactor         - a factor controlling how big the partitions are = nPart goes like 1/pFactor
 * Output:  nPart
 * Returns: new pointer to pp
 */
void setPartitionNumber(int np, double rad, double *aa[], int nPart[], double pFactor) {
    for (int i=0; i<3; i++) {
        if( (int) *aa[i]/(pFactor*rad) > 1 ){
            nPart[i] = (int) *aa[i]/(pFactor*rad);
        } else {
            nPart[i] = 1;
        }
    }
    nPart[3] = np;
}

/*
 * Input:   aa[3]           - length of axis
 *          nPart[3]        -
 * Output:  max[3]          - max dimensions in each direction for any time
 */
void getMaxDimensions(double *aa[], int nPart[], double max[]){
    double finalRadius = pow(*aa[0]* *aa[1]* *aa[2],(1./3.));
    for(int i=0; i<3; i++){
        if( *aa[i] > finalRadius){
            max[i] = *aa[i];
        } else {
            max[i] = finalRadius;
        }
    }
}


/*
 * Input:   np              - particle number
 *          rad             - max particle radius
 *          aa              - length of axes
 *          nPart           - number of partitions
 *          pFactor         - a factor controlling how big the partitions are = nPart goes like 1/pFactor
 * Output:  nPart
 * Returns: new pointer to pp
 */
particle** allocatePP(int np, double rad, double *aa[], int nPart[], particle *pp[], double pFactor) {
    if (pp) {
        free(pp);
    }
    double maxDimensions[3];
    getMaxDimensions(aa,nPart,maxDimensions);
    int maxSize[3];
    for(int i=0; i<3; i++){
        if( maxDimensions[i]/(pFactor*rad) > 1 ){
            maxSize[i] = (int) (maxDimensions[i]/(pFactor*rad));
        } else {
            maxSize[i] = 1;
        }
        //    printf("maxsize: %d\n", maxSize[i]);
    }
    pp = arrayAlloc1dParticlePointer(pp, maxSize[0]*maxSize[1]*maxSize[2]*np); //we can now allocate less than this
    for (int i=0; i<maxSize[0]*maxSize[1]*maxSize[2]*np; i++) {
        pp[i]=NULL;
    }
    
    return pp;
}


int* allocatePCount(double rad, double *aa[], int nPart[], int pCount[], double pFactor) {
    if (pCount) {
        free(pCount);
    }
    double maxDimensions[3];
    getMaxDimensions(aa,nPart,maxDimensions);
    int maxSize[3];
    for(int i=0; i<3; i++){
        if( maxDimensions[i]/(pFactor*rad) > 1 ){
            maxSize[i] = (int) (maxDimensions[i]/(pFactor*rad));
        } else {
            maxSize[i] = 1;
        }
    }
    
    pCount = arrayAlloc1dInt(pCount, maxSize[0]*maxSize[1]*maxSize[2]); //we can now allocate less than this
    for (int i=0; i<maxSize[0]*maxSize[1]*maxSize[2]; i++) {
        pCount[i]=0;
    }
    
    return pCount;
}

//sets pPosArray
void makePartitions(double *aa[], double maxRadius, int nPart[], double pPosArray[] ){
    //vary these functions for non-uniform partitioning
    for(int i=0; i<nPart[0]; i++){
        pPosArray[i] = (i+1) * 2 * (*aa[0] + 2*maxRadius) / nPart[0] - (*aa[0] + 2*maxRadius);
    }
    for(int i=0; i<nPart[1]; i++){
        pPosArray[nPart[0] + i] = (i+1) * 2 * (*aa[1] + 2*maxRadius) / nPart[1] - (*aa[1] + 2*maxRadius);
    }
    for(int i=0; i<nPart[2]; i++){
        pPosArray[nPart[0] + nPart[1] + i] = (i+1) * 2 * (*aa[2] + 2*maxRadius) / nPart[2] - (*aa[2] + 2*maxRadius);
    }
    
}




/*********************
 * partition main    *
 *********************/


void rePartitionParticle(particle *p, int nPart[], int pCount[], particle *pp[], double pPosArray[] ){
    int c[4] = {p->coord[0], p->coord[1], p->coord[2], p->coord[3]}; // c holds p's partition... which might be out of date
    int d[3];
    //findPartitionOld(pp[getIndex4d(c[0],c[1],c[2],c[3],nPart)], nPart, aa, cc); //why pass pp[...]? wouldn't p work? wait, what is the difference between c and cc?
    findPartition(p->position, nPart, pPosArray, d); // d is the updated partition
    // if particle has changed partitions, then repartition
    if( d[0] != c[0] || d[1] != c[1] || d[2] != c[2] ){
        //printf("move from %i %i %i %i\n",c[0],c[1],c[2],c[3]);
        //printf("       to %i %i %i %i\n",d[0],d[1],d[2],pCount[getIndex3dShort(d, nPart)]);
        //printf("repartitioned\n");
        // set particle's partition to new partition
        p->coord[0] = d[0];
        p->coord[1] = d[1];
        p->coord[2] = d[2];
        
        //     if(pCount[getIndex3d(d[0],d[1],d[2],nPart)] == NULL){ printf("trying to read zero value: %i\n", pCount[getIndex3d(d[0],d[1],d[2],nPart)]); exit(1);}
        
        p->coord[3] = pCount[getIndex3dShort(d,nPart)]; // update particle's number in partition
        //printf("coord[3] before and after: %i %i\n",c[3],p->coord[3]);
        pp[getIndex4dShort(p->coord,nPart)] = p; // move particle to end of new partition
        pp[getIndex4dShort(c,nPart)]=NULL; // clear p's old spot in pp
        pCount[getIndex3dShort(d,nPart)] += 1; // update number of particles in new partition
        
        pCount[getIndex3dShort(c,nPart)] -= 1; // update number of particles in old partition
        // if the old partition is not empty now, move last particle to moved particle's old spot
        /*for (int i=0; i<nPart[0]*nPart[1]*nPart[2]; i++) {
         printf("pCount[%i] %i\n",i,pCount[i]);
         }*/
        if( pCount[getIndex3dShort(c,nPart)] != 0 && c[3] != pCount[getIndex3dShort(c, nPart)]) {
            // move particle
            pp[getIndex4dShort(c,nPart)]  = pp[getIndex4d(c[0],c[1],c[2],pCount[getIndex3dShort(c,nPart)],nPart)]; //move last particle to position of old particle.
            // update particle's partition and number in partition
            pp[getIndex4dShort(c,nPart)]->coord[3] = c[3];
            pp[getIndex4d(c[0],c[1],c[2],pCount[getIndex3dShort(c,nPart)],nPart)] = NULL; //clear partition of last particle
            /*for (int i=0; i<4; i++) {
             pp[getIndex4dShort(c,nPart)]->coord[i] = c[i];
             }*/
        }
        
    }
}

void partitionParticles(particle p[], int nPart[], int pCount[], particle *pp[], double pPosArray[]){
    int np = nPart[3];
    int c[4];
    
    //initialize pp, nPart;
    for(int i=0; i<nPart[0]; i++){
        for(int j=0; j<nPart[1]; j++){
            for(int k=0; k<nPart[2]; k++){
                pCount[getIndex3d(i,j,k,nPart)] = 0;
                for(int n=0; n<np; n++){
                    pp[getIndex4d(i,j,k,n,nPart)] = NULL;
                }
            }
        }
    }
    
    for(int n=0; n<np; n++){
        findPartition(p[n].position, nPart, pPosArray, c);
        p[n].coord[0] = c[0];
        p[n].coord[1] = c[1];
        p[n].coord[2] = c[2];
        /*  if( getIndex3dShort(c,nPart) != 0){
         printf("gI3DS(c): %i\n", getIndex3dShort(c,nPart) ); //sometimes returns -1
         printf("c[]: %i %i %i %i\n", c[0], c[1], c[2], c[3]);
         } */
        p[n].coord[3] = pCount[getIndex3dShort(c,nPart)];
        pp[getIndex4dShort(p[n].coord,nPart)] = &p[n]; //check if the pointer stuff is right
        pCount[getIndex3dShort(c,nPart)] += 1;
    }
    
}

//pPosArray  = [x_0, ..., x_n, y_0, ..., y_m, z_0, ..., z_l]
void findPartition(double p[], int nPart[], double pPosArray[] ,int c[] ){
    c[0] = findPartitionPos(p[0], nPart[0], &pPosArray[ 0 ] );
    c[1] = findPartitionPos(p[1], nPart[1], &pPosArray[ nPart[0] ] );
    c[2] = findPartitionPos(p[2], nPart[2], &pPosArray[ nPart[0] + nPart[1] ] );
    //    printf("%i, %i, %i\n", c[0], c[1], c[2]);
}



//this should be remade into a search tree type thing maybe. Depends on average size of nPart I guess.
int findPartitionPos(double position, int nPart, double pPos[]){
    for(int i=0; i<nPart-1; i++){ // check through nPart-2 because nPart-1 is the alternate returned value
        if( position <= pPos[i] ){
            //           printf("pos: %f, nPart: %i, pPos[nPart-1]: %f\n", position, nPart, pPos[nPart-1]);
            return i;
        }
    }
#ifdef DEBUG0
    //note this may occur when finding delta+/delta-
    printf("error!: pos: %f, nPart: %i, pPos[nPart-1]: %f\n", position, nPart, pPos[nPart-1]);
#endif
    return nPart-1; //particle exceeds uppper boundary
}


//returns true if particle is in higher partition than other particle
int greaterPartition(int dx, int dy, int dz){
    return (dx > 0 || ( dx == 0 && dy > 0) || (dx == 0 && dy == 0 && dz >= 0) );
}
//if we restirct dx,dy,dz<n
//return  (n+1)^2*dx + (n+1)*dy + dz >= 0; ?? is this right, is this faster?



/*********************
 * testing functions *
 *********************/

// check all pairs of particle for overlap. print the particle number for overlapping
// pairs and return the number of overlaps
int checkOverlapFull(particle *p, int np) {
    int anyOverlap = 0;
    int overlapCheck;
    for (int i = 0; i < np; i++) {
        for (int j=0; j<i; j++) {
            overlapCheck = overlapQ(&p[i], &p[j]);
            if (overlapCheck) {
                printf("overlap %i %i\n",i,j);
                anyOverlap+=overlapCheck;
                //double dx[3];
                //vectorSubtract(p[i].position, p[j].position, dx);
                //printf("%e\n",vectorNorm(dx)-p[i].rad-p[j].rad);
            }
        }
    }
    return anyOverlap;
}


//determine if pCount contains invalid counts
// make sure it sums to np
// make sure no elements are negative
int checkpCount(int np, int nPart[], int pCount[]) {
    int pSum = 0;
    int errorQ = FALSE;
    for (int i = 0; i < nPart[0]; i++){
        for (int j = 0; j < nPart[1]; j++){
            for (int k = 0; k < nPart[2]; k++){
                pSum+=pCount[getIndex3d(i,j,k,nPart)];
                if(pCount[getIndex3d(i,j,k,nPart)]<0){
                    errorQ = TRUE;
                }
            }
        }
    }
    if (pSum != np) {
        errorQ = TRUE;
    }
    if( errorQ ){
        printf("pCount contains invalid counts!\n");
        printf("pSum= %i\n", pSum);
        printf("i j k count:\n");
        for (int i = 0; i < nPart[0]; i++){
            for (int j = 0; j < nPart[1]; j++){
                for (int k = 0; k < nPart[2]; k++){
                    printf("%i %i %i %i\n",i,j,k,pCount[getIndex3d(i,j,k,nPart)]);
                }
            }
        }
        return TRUE;
    }
    return FALSE;
}

// compare to pCount to p
// look at each p[i].coord and make a tally
int checkpCountP(int np, int nPart[], int pCount[], particle *p){
    int errorQ = FALSE;
    int qCount[nPart[0]*nPart[1]*nPart[2]];
    
    for (int i = 0; i < nPart[0]; i++){
        for (int j = 0; j < nPart[1]; j++){
            for (int k = 0; k < nPart[2]; k++){
                int index[3] = {i,j,k};
                qCount[getIndex3dShort(index, nPart)]=0;
            }
        }
    }
    for(int i=0; i<np; i++){
        qCount[getIndex3dShort(p[i].coord,nPart)] +=1;
    }
    for (int i = 0; i < nPart[0]; i++){
        for (int j = 0; j < nPart[1]; j++){
            for (int k = 0; k < nPart[2]; k++){
                if(pCount[getIndex3d(i,j,k,nPart)] != qCount[getIndex3d(i,j,k,nPart)]){
                    errorQ = TRUE;
                }
            }
        }
    }
    if( errorQ ){
        printf("pCount does not match actual counts based on p!\n");
        printf("i j k count:\n");
        printf("pCount:\n");
        for (int i = 0; i < nPart[0]; i++){
            for (int j = 0; j < nPart[1]; j++){
                for (int k = 0; k < nPart[2]; k++){
                    printf("%i %i %i %i\n",i,j,k,pCount[getIndex3d(i,j,k,nPart)]);
                }
            }
        }
        printf("actual counts:\n");
        printf("i j k count:\n");
        for (int i = 0; i < nPart[0]; i++){
            for (int j = 0; j < nPart[1]; j++){
                for (int k = 0; k < nPart[2]; k++){
                    printf("%i %i %i %i\n",i,j,k,qCount[getIndex3d(i,j,k,nPart)]);
                }
            }
        }
        return TRUE;
    }
    return FALSE;
}


// compare to pCount to pp
// look at each element in pp and make a tally
// look at the full pp array, counting all non-null entries
int checkpCountPP(int np, int nPart[], int pCount[], particle *pp[]){
    int errorQ = FALSE;
    int qCount[nPart[0]*nPart[1]*nPart[2]];
    
    for (int i = 0; i < nPart[0]; i++){
        for (int j = 0; j < nPart[1]; j++){
            for (int k = 0; k < nPart[2]; k++){
                //int index[3] = {i,j,k};
                qCount[getIndex3d(i,j,k, nPart)]=0;
            }
        }
    }
    
    for (int i = 0; i < nPart[0]; i++){
        for (int j = 0; j < nPart[1]; j++){
            for (int k = 0; k < nPart[2]; k++){
                for(int n=0; n<np; n++){ //this assumes no straggling particles //what does this mean?
                    int c[4] = {i,j,k,n};
                    if( pp[getIndex4dShort(c,nPart)] != NULL ){
                        qCount[getIndex3dShort(c,nPart)] +=1;
                    }
                }
            }
        }
    }
    for (int i = 0; i < nPart[0]; i++){
        for (int j = 0; j < nPart[1]; j++){
            for (int k = 0; k < nPart[2]; k++){
                if(pCount[getIndex3d(i,j,k,nPart)] != qCount[getIndex3d(i,j,k,nPart)]){
                    errorQ = TRUE;
                }
            }
        }
    }
    if( errorQ ){
        printf("pCount does not match actual counts based on pp!\n");
        
        printf("i j k count:\n");
        printf("pCount:\n");
        for (int i = 0; i < nPart[0]; i++){
            for (int j = 0; j < nPart[1]; j++){
                for (int k = 0; k < nPart[2]; k++){
                    printf("%i %i %i %i\n",i,j,k,pCount[getIndex3d(i,j,k,nPart)]);
                }
            }
        }
        printf("actual counts:\n");
        printf("i j k count:\n");
        for (int i = 0; i < nPart[0]; i++){
            for (int j = 0; j < nPart[1]; j++){
                for (int k = 0; k < nPart[2]; k++){
                    printf("%i %i %i %i\n",i,j,k,qCount[getIndex3d(i,j,k,nPart)]);
                    /*for (int n=0; n<qCount[getIndex3d(i,j,k,nPart)]; n++) {
                     printf("%i\n",pp[getIndex4d(i, j, k, n, nPart)]);
                     }*/
                }
            }
        }
        return TRUE;
    }
    return FALSE;
}

// check each particle in p to make sure coord points to the correct position in pp
// also check each particle in pp and make sure it has the right coord
int checkPvsPP(int np, int nPart[], int pCount[], particle *p, particle *pp[]) {
    int error=FALSE;
    int c[4];
    
    // for each particle in p, check that coord points to the correct spot in pp
    for (int n=0; n<np; n++) {
        for (int i=0; i<4; i++) {
            c[i]=p[n].coord[i];
        }
        if (pp[getIndex4dShort(c, nPart)] != p+n) {
            printf("p[%i].coord is wrong\n",n);
            printf("coord %i %i %i %i\n",c[0],c[1],c[2],c[3]);
            printf("points to pp containing p[%i]\n",pp[getIndex4dShort(c, nPart)]-p);
            error = TRUE;
        }
        
    }
    
    // for each particle pointer in pp, check that it points to the right particle in p (are these errors symmetric with above?)
    int ppError = FALSE;
    for (int i=0; i<nPart[0]; i++) {
        for (int j=0; j<nPart[1]; j++) {
            for (int k=0; k<nPart[2]; k++) {
                for (int n=0; n<pCount[getIndex3d(i, j, k, nPart)]; n++) {
                    ppError = FALSE;
                    int index[4] = {i,j,k,n};
                    for (int m=0; m<4; m++) {
                        c[m] = pp[getIndex4dShort(index, nPart)]->coord[m];
                        if (c[m] != index[m]) {
                            ppError = TRUE;
                        }
                    }
                    if (ppError) {
                        printf("pp[%i,%i,%i,%i] is wrong\n",i,j,k,n);
                        printf("points to particle with coord %i %i %i %i\n",c[0],c[1],c[2],c[3]);
                        printf("p[%i]\n",pp[getIndex4dShort(c, nPart)]-p);
                        error = TRUE;
                    }
                    
                }
            }
        }
    }
    if (error) {
        printf("pCount:\n");
        for (int i = 0; i < nPart[0]; i++){
            for (int j = 0; j < nPart[1]; j++){
                for (int k = 0; k < nPart[2]; k++){
                    printf("%i %i %i %i\n",i,j,k,pCount[getIndex3d(i,j,k,nPart)]);
                }
            }
        }
    }
    
    return error;
}

// run data structure checks
int runChecks(int np, int nPart[], int pCount[], particle *p, particle *pp[]) {
    int error = 0;
    error += checkpCount(np, nPart, pCount);
    error += checkpCountP(np, nPart, pCount, p);
    error += checkpCountPP(np, nPart, pCount, pp);
    error += checkPvsPP(np, nPart, pCount, p, pp);
    error += checkOverlapFull(p, np);
    if (error) {
        printf("error...\n");
        //exit(1);
    }
    
    return error;
}

// run data structure checks
int runChecksMinusOverlap(int np, int nPart[], int pCount[], particle *p, particle *pp[]) {
    int error = 0;
    error += checkpCount(np, nPart, pCount);
    error += checkpCountP(np, nPart, pCount, p);
    error += checkpCountPP(np, nPart, pCount, pp);
    error += checkPvsPP(np, nPart, pCount, p, pp);
    if (error) {
        printf("error...\n");
        //exit(1);
    }
    
    return error;
}



