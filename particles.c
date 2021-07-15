//
//  particles.c
//  ellipsoidRelaxation
//
//  Created by Chris Burke on 2/28/16.
//  Copyright © 2016 tufts. All rights reserved.
//

#include <float.h>
#include <stdlib.h>
#include <math.h>
#include "surface.h"
#include "particles.h"
#include "partition.h"
#include "forces.h"
#include "misc.h"


/********************
 * particle spacing *
 ********************/

/* find the closest center-to-center distance among all particles
 * Input:       np             - number of particles
 p              - list of particles
 * Output:      None
 * Returns:	closest distance
 */
// note - this function isn't being used. passing of p may be done wrong
// very slow, doesn't use partitioning
double closestDistance(int np, particle *p) {
    double closest = DBL_MAX;
    double sepVector[3];
    double seperation = 0.0;
    
    
    for( int i = 0; i < np; i++) {
        for( int j = 0; j < i; j++) {
            vectorSubtract(p[i].position, p[j].position, sepVector);
            seperation = vectorNorm(sepVector)-(p[i].rad+p[j].rad);
            if (seperation < closest) closest = seperation;
        }
    }
    
    return closest;
}


/* determine whether particles are overlapping
 * Input:	p1		- first particle
 *	       	p2		- second particle
 * Output:      None
 * Returns:	TRUE if overlapping, FALSE otherwise
 */
int overlapQ(particle *p1, particle *p2) {
    double sepVector[3];
    //double seperation;
    double seperationSqr;
    
    
    vectorSubtract(p1->position, p2->position, sepVector);
    //seperation = vectorNorm(sepVector);
    seperationSqr = vectorNormSquare(sepVector);
    
    if (seperationSqr < (p1->rad+p2->rad)*(p1->rad+p2->rad)) {
        return TRUE;
    } else {
        return FALSE;
    }
    
}


/* Undo any overlap between particles. If it takes too long to undo
 * all overlap, signal that the simulation is jammed.
 * Input:   np	            - number of particles
 *          aa               - length of axis
 *          p               - list of particles
 *          pp              - partitioned pointers to p
 *          nCount          -
 *          pCount          - number of particles in a given partition
 *          maxUndoSteps    - jamming criteria - jam if more steps than this are taken
 *          dtOverlap       - timestep
 * Output:  jammed          - whether the simulation is jammed or not
 *          p               - particles will be moved
 */
void undoOverlaps(double *aa[], double maxRadius, particle p[],  int nPart[], int pCount[], particle *pp[], double pPosArray[],
                  int maxUndoSteps, double *dtOverlap, int *nearJammed) {
    int np = nPart[3];
    int totalRelaxOverlap;
    //int particleRelaxOverlap;
    //int currentOverlap;
    int undoSteps = 0;
    *nearJammed = FALSE;
    
    //particle *x[np][np]; // big array, may want to malloc
    
    for (int j = 0; j < np; j++) {
        vectorCopy(p[j].position,p[j].oldPosition);
    }
    
    do {
        undoSteps++;
        // remake pp.
        partitionParticles(p, nPart, pCount, pp, pPosArray);
        
        // reset forces
        for (int k = 0; k < np; k++) {
            for (int l = 0; l < 3; l++) {
                p[k].forceDet[l] = 0;
                p[k].forceStoc[l] = 0;
            }
        }
        
        // calculate forces
        totalRelaxOverlap = applyToPairs(2*maxRadius, TRUE, maxRadius, 0, p, nPart, pCount, pp, pPosArray, &addOverlapForce);
        
        // integrate
        for (int k = 0; k < np; k++) {
            //verletIntegrate(&p[k]);
            gradDescentConstrainedIntegrate(*aa[0],*aa[2],&p[k],*dtOverlap);
            //project(*aa[0], *aa[2], p[k].position, p[k].position);
        }
        
        if (undoSteps%maxUndoSteps==0) {
            //*dtOverlap *= 0.5;
            *nearJammed = TRUE;
            printf("smaller dtOverlap %e\n",*dtOverlap);
        }
        
    } while (totalRelaxOverlap && !*nearJammed);
    //printf("undoSteps %i\n",undoSteps);
    
    partitionParticles(p, nPart, pCount, pp, pPosArray);
    
}
/*
 * Input:   np	            - number of particles
 *          aa              - length of axis
 *          maxRadius       - search radius
 *          p               - particle in question
 *          pCount
 *          pp              -
 *          pPosArray              - partition of particles
 * Output:  TRUE if overlap, FALSE if no overlap
 */
int anyOverlapQ(int nPart[], double maxRadius, particle *p, int pCount[], particle *pp[],  double pPosArray[] ) {
    //int np = nPart[3];
    //int index = 0;
    //    rePartitionParticle(p, nPart, pCount, pp, pPosArray );
    
    int i = p->coord[0];
    int j = p->coord[1];
    int k = p->coord[2];
    int n = p->coord[3];
    
    int c[4] = {i,j,k,n};
    
    int z = 5;
    int deltaPlus[3]={z,z,z};
    int deltaMinus[3]={-z,-z,-z};
    
    double posPlus[3] = {p->position[0] + 1.05*maxRadius,
        p->position[1] + 1.05*maxRadius,
        p->position[2] + 1.05*maxRadius};
    
    double posMinus[3] = {p->position[0] - 1.05*maxRadius,
        p->position[1] - 1.05*maxRadius,
        p->position[2] - 1.05*maxRadius};
    
    // double posPlus[3] = {0,0,0};
    // double posMinus[3] = {0,0,0};
    
    findPartition(posPlus , nPart, pPosArray, deltaPlus);
    findPartition(posMinus , nPart, pPosArray, deltaMinus);
    
    //    printf("%i, %i, %i\n", deltaPlus[0],deltaPlus[1],deltaPlus[2]);
    for(int l=0; l<3; l++){
        deltaPlus[l] -= c[l];
        deltaMinus[l] -= c[l];
        //  printf("d+: %d, d-: %d\n", deltaPlus[l], deltaMinus[l]);
    }
    //    printf("d+: %i, %i, %i\n", deltaPlus[0],deltaPlus[1],deltaPlus[2]);
    //    printf("d-: %i, %i, %i\n", deltaMinus[0],deltaMinus[1],deltaMinus[2]);
    //look in partitions i+di, j+dj, k+dk
    for (int di = deltaMinus[0]; di <= deltaPlus[0]; di++){
        for (int dj = deltaMinus[1]; dj <= deltaPlus[1]; dj++){
            for (int dk = deltaMinus[2]; dk <= deltaPlus[2]; dk++){
                
                //if particle overlaps with any in nearby partitions but isn't itself
                /*				if( (di == 0) && (dj == 0) && (dk == 0) ){
                 for (int m=0; m<pCount[getIndex3d(i,j,k,nPart)]; m++){ // [i+j*nPart[0]+k*nPart[0]*nPart[1]]
                 if( m != pp[getIndex4dShort(c,nPart)]->coord[3]){
                 if( overlapQ(p, pp[getIndex4d(i+di,j+dj,k+dk,m,nPart)]) ) return TRUE;
                 }
                 }
                 } else { //for rest of higher particles in adjascent partitions
                 for (int m=0; m<pCount[getIndex3d(i+di,j+dj,k+dk,nPart)]; m++){
                 if( overlapQ(p, pp[getIndex4d(i+di,j+dj,k+dk,m,nPart)]) ) return TRUE;
                 }
                 }  */
                
                for (int m=0; m<pCount[getIndex3d(i+di,j+dj,k+dk,nPart)]; m++){
                    if( overlapQ(p, pp[getIndex4d(i+di,j+dj,k+dk,m,nPart)]) ){
                        //if( di != 0 || dj != 0 || dk != 0 || m != pp[getIndex4dShort(c,nPart)]->coord[3] ){ // if p != pp[]  doesn't seem to work?
                        if( p != pp[getIndex4d(i+di,j+dj,k+dk,m,nPart)] ){ // if p != pp[]  doesn't seem to work?
                            return TRUE;
                        }
                    }
               	}
            }
        }
    }
    return FALSE;
}

/* Undo any overlap between particles. If it takes too long to undo
 * all overlap, signal that the simulation is jammed.
 * Input:   np	            - number of particles
 *          a               - semi-major/minor axis in z-direction
 *          b               - semi-major/minor axis in xy-plane
 *          p               - list of particles
 *          maxUndoSteps    - jamming criteria - jam if more steps than this are taken
 *          dtOverlap       - timestep
 * Output:  jammed          - whether the simulation is jammed or not
 *          p               - particles will be moved
 */
void undoOverlapsOld(int np, double a, double b, double rad, particle *p,
                     int maxUndoSteps, double *dtOverlap, int *nearJammed) {
    int totalRelaxOverlap;
    int particleRelaxOverlap;
    int currentOverlap;
    int undoSteps = 0;
    *nearJammed = FALSE;
    
    for (int j = 0; j < np; j++) {
        vectorCopy(p[j].position,p[j].oldPosition);
    }
    do {
        undoSteps++;
        // reset forces
        for (int k = 0; k < np; k++) {
            for (int l = 0; l < 3; l++) {
                p[k].forceDet[l] = 0;
                p[k].forceStoc[l] = 0;
            }
        }
        // calculate forces
        totalRelaxOverlap = 0;
        //for (int k = 0; k < np; k++) particleRelaxOverlapTable[k] = 0;
        for (int k = 0; k < np; k++) {
            particleRelaxOverlap = 0;
            for (int l = k+1; l < np; l++) {
                currentOverlap = addOverlapForce(&p[k], &p[l], rad, 0);
                particleRelaxOverlap += currentOverlap;
                //particleRelaxOverlapTable[k] += currentOverlap; //
                //particleRelaxOverlapTable[l] += currentOverlap; // add both of these to account for them later
            }
            totalRelaxOverlap += particleRelaxOverlap;
        }
        // integrate
        for (int k = 0; k < np; k++) {
            //verletIntegrate(&p[k]);
            /*gradDescentIntegrate(&p[k],*dtOverlap);
             project(a, b, p[k].position, p[k].position);*/
            gradDescentConstrainedIntegrate(a,b,&p[k],*dtOverlap);
        }
        //printf("%i\n",relaxOverlap);
        
        
        // debugging - output info to determine jamming criteria
        
        //if (undoSteps % 100 == 0) {
        //printf("undosteps: %i   overlaps: %i\n",undoSteps,totalRelaxOverlap);
        /*if (overlapFile) {
         fprintf(overlapFile, "%i ", totalRelaxOverlap);
         } else {
         printf("overlapFile pointer is null");
         }*/
        //}
        
        if (undoSteps%maxUndoSteps==0) {
            *dtOverlap *= 0.5;
            *nearJammed = TRUE;
            printf("smaller dtOverlap %e\n",*dtOverlap);
        }
        
        
    } while (totalRelaxOverlap && !*nearJammed);
    
    
    //printf("undosteps: %i\n",undoSteps);
    
}

/* Undo any overlap between particles. If it takes too long to undo
 * all overlap, signal that the simulation is jammed.
 * Input:   np	            - number of particles
 *          a               - semi-major/minor axis in z-direction
 *          b               - semi-major/minor axis in xy-plane
 *          p               - list of particles
 *          maxUndoSteps    - jamming criteria - jam if more steps than this are taken
 *          dtOverlap       - timestep
 * Output:  jammed          - whether the simulation is jammed or not
 *          p               - particles will be moved
 */
void undoOverlapsOldProfile(int np, double a, double b, double rad, particle *p,
                            int maxUndoSteps, double *dtOverlap, int *nearJammed, FILE* undoFile) {
    int totalRelaxOverlap;
    int particleRelaxOverlap;
    int currentOverlap;
    int undoSteps = 0;
    *nearJammed = FALSE;
    
    if (undoFile) {
        fprintf(undoFile, "%.16e ", *dtOverlap);
    } else {
        printf("undoFile pointer is null\n");
        exit(1);
    }
    
    for (int j = 0; j < np; j++) {
        vectorCopy(p[j].position,p[j].oldPosition);
    }
    do {
        undoSteps++;
        // reset forces
        for (int k = 0; k < np; k++) {
            for (int l = 0; l < 3; l++) {
                p[k].forceDet[l] = 0;
                p[k].forceStoc[l] = 0;
            }
        }
        // calculate forces
        totalRelaxOverlap = 0;
        //for (int k = 0; k < np; k++) particleRelaxOverlapTable[k] = 0;
        for (int k = 0; k < np; k++) {
            particleRelaxOverlap = 0;
            for (int l = k+1; l < np; l++) {
                currentOverlap = addOverlapForce(&p[k], &p[l], rad, 0);
                particleRelaxOverlap += currentOverlap;
                //particleRelaxOverlapTable[k] += currentOverlap; //
                //particleRelaxOverlapTable[l] += currentOverlap; // add both of these to account for them later
            }
            totalRelaxOverlap += particleRelaxOverlap;
        }
        // integrate
        for (int k = 0; k < np; k++) {
            //verletIntegrate(&p[k]);
            /*gradDescentIntegrate(&p[k],*dtOverlap);
             project(a, b, p[k].position, p[k].position);*/
            gradDescentConstrainedIntegrate(a,b,&p[k],*dtOverlap);
        }
        //printf("%i\n",relaxOverlap);
        
        
        // debugging - output info to determine jamming criteria
        
        //if (undoSteps % 100 == 0) {
        //printf("undosteps: %i   overlaps: %i\n",undoSteps,totalRelaxOverlap);
        /*if (overlapFile) {
         fprintf(overlapFile, "%i ", totalRelaxOverlap);
         } else {
         printf("overlapFile pointer is null");
         }*/
        //}
        
        if (undoSteps%maxUndoSteps==0) {
            *dtOverlap *= 0.5;
            *nearJammed = TRUE;
            printf("smaller dtOverlap %e\n",*dtOverlap);
        }
        
        
    } while (totalRelaxOverlap && !*nearJammed);
    
    if (undoFile) {
        fprintf(undoFile, "%i\n", undoSteps);
    } else {
        printf("undoFile pointer is null\n");
        exit(1);
    }
    
    //printf("undosteps: %i\n",undoSteps);
    
}


/* Use the force on a particle to update its position
 * Input:   p	            - particle to be integrated
 * Output:  p.position	    - particle's new position
 * Returns:	Nothing
 */
void verletIntegrate(particle *p, double dt) {
    double newPosition[3];
    
    for (int i=0; i<3; i++) {
        newPosition[i] = 2*(p->position[i]) - p->oldPosition[i] + (p->forceDet[i])*dt +(p->forceStoc[i])*sqrt(dt);
    }
    vectorCopy(p->position,p->oldPosition);
    vectorCopy(newPosition,p->position);
    
}

/* Use the force on a particle to update its position
 * Input:   p	            - particle to be integrated
 * Output:  p.position	    - particle's new position
 * Returns:	Nothing
 */
void gradDescentIntegrate(particle *p, double dt) {
    
    for (int i=0; i<3; i++) {
        p->position[i] = p->position[i] + (p->forceDet[i])*dt + (p->forceStoc[i])*sqrt(dt);
    }
    
}


/*****************************************************
 * Constraint enforcement using Lagrange multipliers *
 ****************************************************/

/* Use the force on a particle to update its position, subject to the surface constraint
 * (enforced by Lagrange multiplier)
 * Input:  - a              - ellipsoid semi major/minor axis length
 *         - b              - ellipsoid semi major/minor axis length
 *         - p	            - particle to be integrated
 * Output:  p.position	    - particle's new position
 * Returns:	Nothing
 */
void gradDescentConstrainedIntegrate(double a, double b, particle *p, double dt) {
    double newPosition[3];
    double sigmaNew; // constraint function value
    double lambda; // Lagrange muliplier
    double gradOld[3]; // gradient of constraint function at old position
    double gradNew[3]; // gradient of constrain function at new position
    double updateForce[3];
    
    for (int i=0; i<3; i++) {
        newPosition[i] = p->position[i]+ (p->forceDet[i])*dt + (p->forceStoc[i])*sqrt(dt);
    }
    
    sigmaNew = constraintFunction(a, b, newPosition);
    
    do {
        constraintGradient(a,b,newPosition,gradNew);
        constraintGradient(a,b,p->position,gradOld);
        lambda = sigmaNew/(vectorDotProduct(gradNew, gradOld));
        vectorScale(gradOld, -lambda, updateForce);
        vectorAdd(newPosition, updateForce, newPosition);
        sigmaNew = constraintFunction(a, b, newPosition);
    } while (sigmaNew > MACHINE_EPSILON);
    vectorCopy(newPosition,p->position);
    
    
}



/* Find typical interparticle spacing
 * In       - np
 *          - rad
 *          - p
 * Out      - none
 * Returns  - typical interaprticle spacing
 * slow, not using partitions
 */
double findTypicalSpacing(int np, double rad, particle *p, int *nPart, int *pCount, particle **pp, double *pPosArray) {
    
    // get 3 smallest spacings for each particle
    double *spacings;
    spacings = arrayAlloc1d(spacings, 3*np);
    double *distanceList;
    distanceList = arrayAlloc1d(distanceList, np);
    int nNeighbors=0;
    //int deficit=0; // account for any particles with < 3 neighbors
    for (int i=0; i<np; i++) {
        for (int j=0; j<3; j++) {
            spacings[3*i+j] = 100*rad;
        }
        
        distancesToNeighbors(nPart, 4*rad, &p[i], pCount, pp, pPosArray, &nNeighbors, distanceList);
        for (int j=0; j<nNeighbors; j++) {
            for (int k=0; k<3; k++) {
                if (distanceList[j] < spacings[3*i+k]) {
                    for (int m=2; m>k; m--) {
                        spacings[3*i+m]=spacings[3*i+m-1];
                    }
                    spacings[3*i+k]=distanceList[j];
                    break;
                }
            }
        }
        // account for particles with less than three neighbors within the search distance - give search distance as distance
        if (nNeighbors<3) {
            for (int j=3; j>nNeighbors; j--) {
                spacings[3*i+j-1]=4*rad;
                //deficit++;
            }
        }
        
    }
    
    // geometric mean
    double logSum=0;
    for (int i=0; i<3*np; i++) {
        logSum+=log(spacings[i]);
    }
    //double typicalSpacing = exp(logSum/(3*np-deficit));
    double typicalSpacing = exp(logSum/(3*np));
    
    free(spacings);
    free(distanceList);
    
    return typicalSpacing;
}


/* Find typical interparticle spacing
 * In       - np
 *          - rad
 *          - p
 * Out      - none
 * Returns  - typical interaprticle spacing
 * slow, not using partitions
 */
double findTypicalSpacingOld(int np, double rad, particle *p) {
    
    double *spacings;
    spacings = arrayAlloc1d(spacings, 3*np);
    for (int i=0; i<np; i++) {
        for (int j=0; j<3; j++) {
            spacings[3*i+j] = 100*rad;
        }
        
        
        for (int j=0; j<np; j++) {
            if (j!=i) {
                double deltaX[3];
                vectorSubtract(p[i].position, p[j].position, deltaX);
                double distance = vectorNorm(deltaX) - p[i].rad - p[j].rad;
                for (int k=0; k<3; k++) {
                    if (distance < spacings[3*i+k]) {
                        for (int m=2; m>k; m--) {
                            spacings[3*i+m]=spacings[3*i+m-1];
                        }
                        spacings[3*i+k]=distance;
                        break;
                    }
                }
            }
        }
        
    }
    
    // geometric mean
    double logSum=0;
    for (int i=0; i<3*np; i++) {
        logSum+=log(spacings[i]);
    }
    double typicalSpacing = exp(logSum/(3*np));
    
    free(spacings);
    
    return typicalSpacing;
}


/* For the given particle, look for neighbors within some cutoff distance and return
 * the interparticle distances to them
 * In      - nPart, pCount, pp:     - partitioning data structures
 *           maxRadius              - center-to-center particle distance to search
 * Out     - nNeighbors             - number of neighbors within the search distance
 *         - distances              - interparticle distances to neighbors, should be length np
 * Returns - TRUE on success, FALSE otherwise
 * NOTE: Currently being used in startup diffusion, not in actual simulation
 */
int distancesToNeighbors(int nPart[], double maxRadius, particle *p, int pCount[], particle *pp[],  double pPosArray[], int *nNeighbors, double *distanceList) {
    
    int i = p->coord[0];
    int j = p->coord[1];
    int k = p->coord[2];
    int n = p->coord[3];
    
    int c[4] = {i,j,k,n};
    
    int z = 5;
    int deltaPlus[3]={z,z,z};
    int deltaMinus[3]={-z,-z,-z};
    
    double posPlus[3] = {p->position[0] + 1.05*maxRadius,
        p->position[1] + 1.05*maxRadius,
        p->position[2] + 1.05*maxRadius};
    
    double posMinus[3] = {p->position[0] - 1.05*maxRadius,
        p->position[1] - 1.05*maxRadius,
        p->position[2] - 1.05*maxRadius};
    
    findPartition(posPlus , nPart, pPosArray, deltaPlus);
    findPartition(posMinus , nPart, pPosArray, deltaMinus);
    
    for(int l=0; l<3; l++){
        deltaPlus[l] -= c[l];
        deltaMinus[l] -= c[l];
    }
    
    *nNeighbors = 0;
    double dx[3];
    double distance;
    //look in partitions i+di, j+dj, k+dk
    for (int di = deltaMinus[0]; di <= deltaPlus[0]; di++){
        for (int dj = deltaMinus[1]; dj <= deltaPlus[1]; dj++){
            for (int dk = deltaMinus[2]; dk <= deltaPlus[2]; dk++){
                for (int m=0; m<pCount[getIndex3d(i+di,j+dj,k+dk,nPart)]; m++){
                    vectorSubtract(p->position, pp[getIndex4d(i+di,j+dj,k+dk,m,nPart)]->position, dx);
                    distance = vectorNorm(dx);
                    if( distance < maxRadius ){
                        if( p != pp[getIndex4d(i+di,j+dj,k+dk,m,nPart)] ){
                            distanceList[*nNeighbors] = distance - p->rad - pp[getIndex4d(i+di,j+dj,k+dk,m,nPart)]->rad;
                            //*nNeighbors = *nNeighbors + 1;
                            (*nNeighbors)++;
                        }
                    }
               	}
            }
        }
    }
    return FALSE;

    
}


/* For the given particle, take a random step in the tangent plane
 * In      - a      - semi-major/minor axis of the ellipsoid
 *           b      -    "                            "
 *         - delta  - average size of the step
 *           xin    - position of the particle to move
 * Out     - xout   - new position
 * Returns - TRUE on success, FALSE otherwise
 * NOTE: Currently being used in startup diffusion, not in actual simulation
 */
int step(double a, double b, double delta, const double *xin, double *xout) {
    
    double d1;
    double d2;
    double t1[3];
    double t2[3];
    double t1Scaled[3];
    double t2Scaled[3];
    double stepVector[3];
    
    // unit tangent vectors
    tangent1(a, b, xin, t1);
    tangent2(a, b, xin, t2);
    
    // step size along each tangent
    d1 = delta*randNormal();
    d2 = delta*randNormal();
    
    // scale tangent vectors accordingly
    vectorScale(t1, d1, t1Scaled);
    vectorScale(t2, d2, t2Scaled);
    
    // compose step vector
    vectorAdd(t1Scaled, t2Scaled, stepVector);
    
    // take step
    vectorAdd(xin, stepVector, xout);
    
    return TRUE;
}

