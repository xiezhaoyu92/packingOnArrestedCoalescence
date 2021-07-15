//
//  forces.c
//  ellipsoidRelaxation
//
//  Created by Chris Burke on 2/28/16.
//  Copyright © 2016 tufts. All rights reserved.
//

#include <stdlib.h>
#include <math.h>
#include "surface.h"
#include "particles.h"
#include "partition.h"
#include "misc.h"


/* Add an attractive force between two particles.
 * The force magnitude goes like 1/seperation and is along the
 * particle seperation vector
 * Input:   p1	            - first particle
 *          p2              - second particle
 *          attrStrength    - factor that determines strength of attractive force
 * Output:  p1.force	    - add force to p1
 *          p2.force        - add force to p2
 *          newForce        - components of the added force - for debugging
 * to do - specify strength based on system's physical parameters
 */
int addAttractiveForce(particle *p1, particle *p2, double attrStrength, double newForce[3]) {
    double deltaX[3]; //seperation vector
    double distance; // magnitude of seperation
    
    vectorSubtract(p1->position, p2->position, deltaX);
    distance = vectorNorm(deltaX);
    
    
    for (int i=0; i<3; i++){
        newForce[i] = attrStrength*deltaX[i]/(distance*distance);
        p1->forceDet[i] += -newForce[i];
        p2->forceDet[i] += newForce[i];
    }
    
    return 0;
}

/* Add the attractive component of a Lennard-Jones force between two particles.
 * The force magnitude goes like 1/r^7 and is along the
 * particle seperation vector
 * Input:   p1	            - first particle
 *          p2              - second particle
 *          attrStrength    - factor that determines strength of attractive force
 * Output:  p1.force	    - add force to p1
 *          p2.force        - add force to p2
 *          newForce        - components of the added force - for debugging
 * to do - specify strength based on system's physical parameters
 */
int addAttractiveLJForce(particle *p1, particle *p2, double attrStrength, double newForce[3]) {
    double deltaX[3]; //seperation vector
    double distance; // magnitude of seperation
    
    vectorSubtract(p1->position, p2->position, deltaX);
    distance = vectorNorm(deltaX);
    
    
    for (int i=0; i<3; i++){
        newForce[i] = attrStrength*6*deltaX[i]/pow(distance,8);
        p1->forceDet[i] += -newForce[i];
        p2->forceDet[i] += newForce[i];
    }
    
    return 0;
}


/* Add a repulsive force between two particles.
 * The force magnitude goes like 1/seperation and is along the
 * particle seperation vector
 * Input:   p1	            - first particle
 *          p2              - second particle
 *          attrStrength    - factor that determines strength of attractive force
 * Output:  p1.force	    - add force to p1
 *          p2.force        - add force to p2
 *          newForce        - components of the added force - for debugging
 * to do - specify strength based on system's physical parameters
 */
int addRepulsiveForce(particle *p1, particle *p2, double attrStrength, double newForce[3]) {
    double deltaX[3]; //seperation vector
    double distance; // magnitude of seperation
    
    vectorSubtract(p1->position, p2->position, deltaX);
    distance = vectorNorm(deltaX);
    
    
    for (int i=0; i<3; i++){
        newForce[i] = -attrStrength*deltaX[i]/(distance*distance);
        p1->forceDet[i] += -newForce[i];
        p2->forceDet[i] += newForce[i];
    }
    
    return 0;
}

/* Calculate the pair-potential of two particles with a repulsive interaction.
 * This is the potential corresponding to the force in addRepulsiveForce()
 * Input:   p1	            - first particle
 *          p2              - second particle
 *          repStrength     - factor that determines strength of repulsive force
 * Output:  none
 * Returns: particle pair potential
 * to do - specify strength based on system's physical parameters
 */
double repulsiveEnergy(particle *p1, particle *p2, double repStrength) {
    double deltaX[3]; //seperation vector
    double distance; // magnitude of seperation
    
    vectorSubtract(p1->position, p2->position, deltaX);
    distance = vectorNorm(deltaX);
    
    return -repStrength*log(distance);
}


/* Add a screened attractive force between two particles.
 * The unscreened force magnitude goes like 1/seperation and is along the
 * particle seperation vector.  The screening is a quadratic dropoff to zero at the given cutoff
 * Input:   p1	            - first particle
 *          p2              - second particle
 *          attrStrength    - factor that determines strength of attractive force
 *          screeningLength - determines the lengthscale over which the attractive force decays -- in terms
 *                              of distance between particle surfaces, not including radius
 * Output:  p1.force	    - add force to p1
 *          p2.force        - add force to p2
 * to do - specify strength based on system's physical parameters
 */
int addScreenedAttractiveForce(particle *p1, particle *p2, double screeningLength, double attrStrength) {
    double deltaX[3]; //seperation vector
    double distance; // magnitude of seperation
    double newForce[3];
    
    vectorSubtract(p1->position, p2->position, deltaX);
    distance = vectorNorm(deltaX);
    double sep = distance - p1->rad - p2->rad;
    
    if (sep < screeningLength) {
        double magnitude = -attrStrength/distance*(1-sep/screeningLength);
        if (fabs(magnitude)>1) {
            printf("big mag %e\n",magnitude);
        }
        for (int i=0; i<3; i++){
            newForce[i] = magnitude*deltaX[i]/distance;
            p1->forceDet[i] += newForce[i];
            p2->forceDet[i] += -newForce[i];
        }
    }
    
    return 0;
}


/* Add a screened attractive Lennard-Jones force between two particles.
 * The unscreened force magnitude goes like 1/seperation and is along the
 * particle seperation vector. The screening is a quadratic dropoff to zero at the given cutoff
 * Input:   p1	            - first particle
 *          p2              - second particle
 *          attrStrength    - factor that determines strength of attractive force
 *          screeningLength - determines the lengthscale over which the attractive force decays -- in terms
 *                              of distance between particle surfaces, not including radius
 * Output:  p1.force	    - add force to p1
 *          p2.force        - add force to p2
 *          newForce        - components of the added force - for debugging
 * to do - specify strength based on system's physical parameters
 */
int addScreenedAttractiveLJForce(particle *p1, particle *p2, double screeningLength, double attrStrength) {
    double deltaX[3]; //seperation vector
    double distance; // magnitude of seperation
    double newForce[3]; 
    
    vectorSubtract(p1->position, p2->position, deltaX);
    distance = vectorNorm(deltaX);
    double sep = distance - p1->rad - p2->rad;
    
    if (sep < screeningLength) {
        double magnitude = attrStrength*6/pow(distance,7)*(1-sep*sep/(screeningLength*screeningLength));
        for (int i=0; i<3; i++){
            newForce[i] = magnitude*deltaX[i]/distance;
            p1->forceDet[i] += -newForce[i];
            p2->forceDet[i] += newForce[i];
        }
    }
    
    
    return 0;
}

#ifdef USE_GSL
/* Calculate the pair-potential of two particles with a screened attractive interaction.
 * This is the potential corresponding to the force in addScreenedAttractiveForce()
 * Input:   p1	            - first particle
 *          p2              - second particle
 *          attrStrength    - factor that determines strength of attractive force
 *          screeningLength - determines the lengthscale over which the attractive force decays
 * Output:  none
 * Returns: particle pair potential
 * to do - specify strength based on system's physical parameters
 */
double screenedAttractiveEnergy(particle *p1, particle *p2, double screeningLength, double attrStrength) {
    double deltaX[3]; //seperation vector
    double distance; // magnitude of seperation
    
    vectorSubtract(p1->position, p2->position, deltaX);
    distance = vectorNorm(deltaX);
    
    return attrStrength*gsl_sf_expint_Ei(-distance/screeningLength);
}
#endif

/* Add gravitational force to a particle
 * Input:   p               - particle
 *          gravStrength    - constant determining magnitude of force
 * Output:  p.force         - add force to p
 */
int addGravitationalForce(particle *p, double gravStrength) {
    p->forceDet[1] -= gravStrength;
    
    return 0;
}

/* Add a stochastic force to a particle
 * Should be proportional to sqrt(dt*2D)
 * Input:   a          - semi-major/minor axis of the ellipsoid
 *          b          -    "                            "
 *          D          - diffusion constant for a particle with radius rad (below)
 *          rad        - radius by which the scale of the diffusion constant is set
 *          dt         - timestep
 *          p          - particle
 * Output:  p.force    - add force to p
 */
int addStochasticForce(double a, double b, double D, double rad, particle *p) {
    // old method: add force along tangent plane
    /*
     double d1;
     double d2;
     double t1[3];
     double t2[3];
     double t1Scaled[3];
     double t2Scaled[3];
     double forceVector[3];
     
     // unit tangent vectors
     tangent1(a, b, p->position, t1);
     tangent2(a, b, p->position, t2);
     
     // step size along each tangent
     d1 = sqrt(D*rad/(p->rad))*randNormal();
     d2 = sqrt(D*rad/(p->rad))*randNormal();
     
     // scale tangent vectors accordingly
     vectorScale(t1, d1, t1Scaled);
     vectorScale(t2, d2, t2Scaled);
     
     // compose step vector
     vectorAdd(t1Scaled, t2Scaled, forceVector);
     */
    
    // make move in 3D
    for (int i=0; i<3; i++){
        p->forceStoc[i] += sqrt(D*rad/(p->rad))*randNormal();
    }
    
    return 0;
}

/* Add a specified force to a particle
 * Input:   p               - particle
 *          f               - a force to apply
 *          strength        - a scaling factor for the force
 * Output:  p.force         - add force to p
 */
int addGenericForce(particle *p, double f[3], double strength) {
    for (int i=0; i<3; i++) {
        p->forceDet[i] += strength*f[i];
    }
    
    return 0;
}

/* Add a force to two particles if they are overlapping
 * The force magnitude is fixed and is a repulsive force along the
 * particle seperation vector
 * Input:   p1	            - first particle
 *          p2              - second particle
 * Output:  p1.force	    - add force to p1
 *          p2.force        - add force to p2
 * Returns:	-1 if particle are overlapping, 0 otherwise
 */
int addOverlapForce(particle *p1, particle *p2, double rad, double strength) {
    double deltaX[3]; //seperation vector
    double distance; // magnitude of seperation
    int isOverlap = 0;
    
    if (overlapQ(p1, p2)) {
        isOverlap = -1;
        vectorSubtract(p1->position, p2->position, deltaX);
        distance = vectorNorm(deltaX);
        
        for (int i=0; i<3; i++){
            p1->forceDet[i] += rad*deltaX[i]/distance;
            p2->forceDet[i] += -rad*deltaX[i]/distance;
        }
        
    }
    
    return isOverlap;
}

/* Add a Hertzian contact force to two particles if they are overlapping
 * Force goes like d^(5/2) where d is the overlap divided by particle radius
 * Input:   p1	            - first particle
 *          p2              - second particle
 *          strength        - coefficient determing force strength
 *          rad             - max particle radius, to set length scale
 * Output:  p1.force	    - add force to p1
 *          p2.force        - add force to p2
 * Returns:	-1 if particle are overlapping, 0 otherwise
 */
int addHertzianForce(particle *p1, particle *p2, double rad, double strength) {
    double deltaX[3]; //seperation vector
    // should the magnitude of this force be based on the particle radius? set here, or externally when determining strength?
    int isOverlap = FALSE;
    
    if (overlapQ(p1, p2)) {
        isOverlap = TRUE;
        vectorSubtract(p1->position, p2->position, deltaX);
        double distance = vectorNorm(deltaX);
        double sigma = p1->rad + p2->rad;
        double d = 1 - distance/sigma;
        
        double mag = strength/sigma*pow(d, 3./2.);
        
        /*if (mag>0.1) {
            printf("large force: %e\n",mag);
            printf("   distance: %e\n",distance);
            printf("   overlap:  %e\n",d);
            //blah
        }*/
        
        for (int i=0; i<3; i++){
            p1->forceDet[i] += deltaX[i]/distance*mag;
            p2->forceDet[i] += -deltaX[i]/distance*mag;
        }
        
    }
    
    return isOverlap;
}

/* Calculate Hertzian contact potential to two particles if they are overlapping
 * Energy goes like d^(5/2) where d is the overlap divided by particle radius
 * Input:   p1	            - first particle
 *          p2              - second particle
 *          strength        - coefficient determing force strength
 *          rad             - max particle radius, to set length scale
 * Output:  p1.force	    - add force to p1
 *          p2.force        - add force to p2
 * Returns:	-1 if particle are overlapping, 0 otherwise
 */
double hertzianEnergy(particle *p1, particle *p2, double rad, double strength) {
    double deltaX[3]; //seperation vector
    // should the magnitude of this force be based on the particle radius? set here, or externally when determining strength?
    //int isOverlap = FALSE;
    double U=0;
    
    if (overlapQ(p1, p2)) {
        //isOverlap = TRUE;
        vectorSubtract(p1->position, p2->position, deltaX);
        double distance = vectorNorm(deltaX);
        double d = 1 - distance/(p1->rad + p2->rad);
        
        U = (2./5.)*strength*pow(d, 5./2.);
        
    }
    
    return U;
}

/* Add a harmonic surface constraint force
 * Force goes like the distance to the surface
 * Input:   a, b            - surface axes lengths
 *          p1	            - particle
 *          strength        - coefficient determing force strength\
 * Output:  p1.forceDet	    - add force to p1
 * Returns:	-1 if particle are overlapping, 0 otherwise
 */
int addHarmonicSurfaceForce(double a, double b, particle *p1, double strength) {
    double dx[3];
    vectorFromSurface(a, b, p1->position, dx);
    
    for (int i=0; i<3; i++) {
        p1->forceDet[i] += -strength*dx[i];
    }
    
    return 0;
}

/* Calculate a harmonic surface constraint potential
 * Force goes like the distance to the surface
 * Input:   a, b            - surface axes lengths
 *          p1	            - particle
 *          strength        - coefficient determing force strength\
 * Output:  p1.forceDet	    - add force to p1
 */
double harmonicSurfaceEnergy(double a, double b, particle *p1, double strength) {
    double dx[3];
    vectorFromSurface(a, b, p1->position, dx);
    
    return -0.5*strength*vectorDotProduct(dx, dx);
}


//applies generic function to all particle pairs in partitions overlapping a specified search radius
//returns number of successful overlaps
//generic function should return "int -1" if it was applied to the pair, "int 0" otherwise
//searchRadius is the force cutoff distance -- this is fed to pairInteraction
//centerCenterCutoff is a boolean specifying wehther searchRadius includes the particle radii (i.e. is a center-to-center distance) or whether it is a surface-to-surface distance
// rad is the maximum particle radius
int applyToPairs(double searchRadius, int centerCenterCutoff, double rad, double strength, particle p[],  int nPart[], int pCount[], particle *pp[], double pPosArray[],
                 int (*pairInteraction)(particle*, particle*, double, double)) {
    double searchRadiusCC; // how far to look, including the maximum possible particle radii
    if (centerCenterCutoff) {
        searchRadiusCC = searchRadius;
    } else {
        searchRadiusCC = searchRadius+2*rad;
    }
    int particleApplications;
    int totalApplications = 0;
    for (int i = 0; i<nPart[0]; i++){
        for (int j = 0; j<nPart[1]; j++){
            for (int k = 0; k<nPart[2]; k++){
                for (int n=0; n<pCount[getIndex3d(i,j,k,nPart)]; n++){ //for every particle
                    particleApplications = 0;
                    int c[4] = {i,j,k,n};
                    //      printf("particle: %i%i%i%i\n",i,j,k,n);
                    if(pp[getIndex4dShort(c,nPart)] == NULL){ printf("trying to read null value!!\n"); exit(1);}
                    
                    int x = 5;
                    int deltaPlus[3]={x,x,x};
                    int deltaMinus[3]={-x,-x,-x};
                    
                    double posPlus[3] = {pp[getIndex4dShort(c,nPart)]->position[0] + 1.05*searchRadiusCC,
                        pp[getIndex4dShort(c,nPart)]->position[1] + 1.05*searchRadiusCC,
                        pp[getIndex4dShort(c,nPart)]->position[2] + 1.05*searchRadiusCC};
                    
                    double posMinus[3] = {pp[getIndex4dShort(c,nPart)]->position[0] - 1.05*searchRadiusCC,
                        pp[getIndex4dShort(c,nPart)]->position[1] - 1.05*searchRadiusCC,
                        pp[getIndex4dShort(c,nPart)]->position[2] - 1.05*searchRadiusCC};
                    
                    findPartition(posPlus , nPart, pPosArray, deltaPlus);
                    findPartition(posMinus , nPart, pPosArray, deltaMinus);
                    
                    for(int l=0; l<3; l++){
                        deltaPlus[l] -= c[l];
                        deltaMinus[l] -= c[l];
                    }
                    
                    //look in partitions i+di, j+dj, k+dk
                    for (int di = deltaMinus[0]; di <= deltaPlus[0]; di++){
                        for (int dj = deltaMinus[1]; dj <= deltaPlus[1]; dj++){
                            for (int dk = deltaMinus[2]; dk <= deltaPlus[2]; dk++){
                                if( greaterPartition(di, dj, dk) ){ // does this really need to be another function?
                                    //#ifdef  DEBUG0
                                    if( (c[0]+di < 0) || (c[1]+dj < 0) || (c[2]+dk < 0) ||
                                       (c[0]+di >= nPart[0]) || (c[1]+dj >= nPart[1]) || (c[2]+dk >= nPart[2])  ){
                                        printf("delta out of range in undoOverlaps\n");
                                        printf("di: %i, dj: %i, dk: %i\n", di,dj,dk);
                                        printf("ci: %i, cj: %i, ck: %i\n", c[0],c[1],c[2]);
                                        printf("pos[0]: %f, pos[1]: %f pos[3]: %f\n", (pp[getIndex4dShort(c,nPart)]->position[0]),
                                               (pp[getIndex4dShort(c,nPart)]->position[1]),
                                               (pp[getIndex4dShort(c,nPart)]->position[2]) );
                                        
                                    } else {
                                        //#endif
                                        if( (di == 0) && (dj == 0) && (dk == 0) ){
                                            for (int m=n+1; m<pCount[getIndex3d(i,j,k,nPart)]; m++){ // [i+j*nPart[0]+k*nPart[0]*nPart[1]]
                                                particleApplications += pairInteraction(pp[getIndex4d(i,j,k,n,nPart)], pp[getIndex4d(i,j,k,m,nPart)], searchRadius,strength);
                                                
                                            }
                                        } else { //for rest of higher particles in adjascent partitions
                                            for (int m=0; m<pCount[getIndex3d(i+di,j+dj,k+dk,nPart)]; m++){
                                                particleApplications += pairInteraction(pp[getIndex4d(i,j,k,n,nPart)], pp[getIndex4d(i+di,j+dj,k+dk,m,nPart)], searchRadius,strength);
                                            }
                                        }
                                        //#ifdef  DEBUG0
                                    }
                                    //#endif
                                }
                            }
                        }
                        totalApplications += particleApplications;
                    }
                }
            }
        }
    }
    return totalApplications;
}

//calculates energy from generic pair interaction function for all particle pairs in partitions overlapping a specified search radius
//returns number of successful overlaps
//generic function should return "int -1" if it was applied to the pair, "int 0" otherwise
//searchRadius is the force cutoff distance -- this is fed to pairInteraction
//centerCenterCutoff is a boolean specifying wehther searchRadius includes the particle radii (i.e. is a center-to-center distance) or whether it is a surface-to-surface distance
// rad is the maximum particle radius
double totalPairEnergy(double searchRadius, int centerCenterCutoff, double rad, double strength, particle p[], int nPart[], int pCount[], particle *pp[], double pPosArray[],
                       double (*pairEnergy)(particle*, particle*, double, double)) {
    double energyTotal = 0.;
    double searchRadiusCC; // how far to look, including the maximum possible particle radii
    if (centerCenterCutoff) {
        searchRadiusCC = searchRadius;
    } else {
        searchRadiusCC = searchRadius+2*rad;
    }
    int particleApplications;
    int totalApplications = 0;
    for (int i = 0; i<nPart[0]; i++){
        for (int j = 0; j<nPart[1]; j++){
            for (int k = 0; k<nPart[2]; k++){
                for (int n=0; n<pCount[getIndex3d(i,j,k,nPart)]; n++){ //for every particle
                    particleApplications = 0;
                    int c[4] = {i,j,k,n};
                    //      printf("particle: %i%i%i%i\n",i,j,k,n);
                    if(pp[getIndex4dShort(c,nPart)] == NULL){ printf("trying to read null value!!\n"); exit(1);}
                    
                    int x = 5;
                    int deltaPlus[3]={x,x,x};
                    int deltaMinus[3]={-x,-x,-x};
                    
                    double posPlus[3] = {pp[getIndex4dShort(c,nPart)]->position[0] + 1.05*searchRadiusCC,
                        pp[getIndex4dShort(c,nPart)]->position[1] + 1.05*searchRadiusCC,
                        pp[getIndex4dShort(c,nPart)]->position[2] + 1.05*searchRadiusCC};
                    
                    double posMinus[3] = {pp[getIndex4dShort(c,nPart)]->position[0] - 1.05*searchRadiusCC,
                        pp[getIndex4dShort(c,nPart)]->position[1] - 1.05*searchRadiusCC,
                        pp[getIndex4dShort(c,nPart)]->position[2] - 1.05*searchRadiusCC};
                    
                    findPartition(posPlus , nPart, pPosArray, deltaPlus);
                    findPartition(posMinus , nPart, pPosArray, deltaMinus);
                    
                    for(int l=0; l<3; l++){
                        deltaPlus[l] -= c[l];
                        deltaMinus[l] -= c[l];
                    }
                    
                    //look in partitions i+di, j+dj, k+dk
                    for (int di = deltaMinus[0]; di <= deltaPlus[0]; di++){
                        for (int dj = deltaMinus[1]; dj <= deltaPlus[1]; dj++){
                            for (int dk = deltaMinus[2]; dk <= deltaPlus[2]; dk++){
                                if( greaterPartition(di, dj, dk) ){ // does this really need to be another function?
                                    //#ifdef  DEBUG0
                                    if( (c[0]+di < 0) || (c[1]+dj < 0) || (c[2]+dk < 0) ||
                                       (c[0]+di >= nPart[0]) || (c[1]+dj >= nPart[1]) || (c[2]+dk >= nPart[2])  ){
                                        printf("delta out of range in undoOverlaps\n");
                                        printf("di: %i, dj: %i, dk: %i\n", di,dj,dk);
                                        printf("ci: %i, cj: %i, ck: %i\n", c[0],c[1],c[2]);
                                        printf("pos[0]: %f, pos[1]: %f pos[3]: %f\n", (pp[getIndex4dShort(c,nPart)]->position[0]),
                                               (pp[getIndex4dShort(c,nPart)]->position[1]),
                                               (pp[getIndex4dShort(c,nPart)]->position[2]) );
                                        
                                    } else {
                                        //#endif
                                        if( (di == 0) && (dj == 0) && (dk == 0) ){
                                            for (int m=n+1; m<pCount[getIndex3d(i,j,k,nPart)]; m++){ // [i+j*nPart[0]+k*nPart[0]*nPart[1]]
                                                energyTotal += pairEnergy(pp[getIndex4d(i,j,k,n,nPart)], pp[getIndex4d(i,j,k,m,nPart)], searchRadius, strength);
                                                
                                            }
                                        } else { //for rest of higher particles in adjascent partitions
                                            for (int m=0; m<pCount[getIndex3d(i+di,j+dj,k+dk,nPart)]; m++){
                                                energyTotal += pairEnergy(pp[getIndex4d(i,j,k,n,nPart)], pp[getIndex4d(i+di,j+dj,k+dk,m,nPart)], searchRadius,strength);
                                            }
                                        }
                                        //#ifdef  DEBUG0
                                    }
                                    //#endif
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return energyTotal;
}

