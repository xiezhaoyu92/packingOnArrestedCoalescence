//
//  DropletCoalescence.c
//  ellipsoidRelaxation
//
//  Created by Zhaoyu Xie on 7/11/16.
//  Copyright © 2016 tufts. All rights reserved.
//

#include "DropletCoalescence.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "surface.h"
#include "partition.h"
#include "particles.h"
#include "forces.h"
#include "misc.h"
#include "mt19937ar.h"

#ifndef MACHINE_EPSILON
#define MACHINE_EPSILON 1e-14
#endif

double DCconstraintFunction(double z0, double a, double b, double *x){
    return (a*a*(-x[2]+z0)*(-x[2]+z0)+b*b*(x[0]*x[0]+x[1]*x[1]))/pow(-x[2]*x[2]+x[2]*z0-(x[0]*x[0]+x[1]*x[1]), 2)-1;
}

void DCconstraintGradient(double z0, double a, double b, double *x, double *grad){
    grad[0]=(4*b*b*x[0]*(x[0]*x[0]+x[1]*x[1])+4*a*a*x[0]*(-x[2]+z0)*(-x[2]+z0))/pow(-x[0]*x[0]-x[1]*x[1]-x[2]*x[2]+x[2]*z0,3)+2*b*b*x[0]/pow(-x[0]*x[0]-x[1]*x[1]-x[2]*x[2]+x[2]*z0,2);
    grad[1]=(4*b*b*x[1]*(x[0]*x[0]+x[1]*x[1])+4*a*a*x[1]*(-x[2]+z0)*(-x[2]+z0))/pow(-x[0]*x[0]-x[1]*x[1]-x[2]*x[2]+x[2]*z0,3)+2*b*b*x[1]/pow(-x[0]*x[0]-x[1]*x[1]-x[2]*x[2]+x[2]*z0,2);
    grad[2]=-(2*b*b*(x[0]*x[0]+x[1]*x[1])*(-2*x[2]+z0)+2*a*a*(-2*x[2]+z0)*(-x[2]+z0)*(-x[2]+z0))/pow(-x[0]*x[0]-x[1]*x[1]-x[2]*x[2]+x[2]*z0,3)-2*a*a*(-x[2]+z0)/pow(-x[0]*x[0]-x[1]*x[1]-x[2]*x[2]+x[2]*z0,2);
}

void DCRelaxationProject(double z0, double aOld, double bOld, double a, double b, particle *p){
    double newPosition[3];
    double gradOld[3];
    double gradNew[3];
    double lambda;
    double normalDirection[3];
    double updateForce[3];
    double sigma;
    int i;
    double normalization;
    double oldSigma;

/*    DCconstraintGradient(z0, aOld, bOld, p->position, normalDirection);
    normalization=vectorNorm(normalDirection);
    for(i=0;i<3;i++){
        normalDirection[i]=normalDirection[i]/normalization;
    }
    for(i=0;i<3;i++){
        newPosition[i]=p->position[i]+1/2000000*normalDirection[i];
    }*/
    for(i=0;i<3;i++){
        newPosition[i]=p->position[i];
    }
    sigma=DCconstraintFunction(z0, a, b, newPosition);
    oldSigma=sigma;
    DCconstraintGradient(z0, aOld, bOld, p->position, gradOld);
    while(fabs(sigma) > MACHINE_EPSILON){
        DCconstraintGradient(z0, a, b, newPosition, gradNew);
        lambda = sigma/(vectorDotProduct(gradNew, gradOld));
        vectorScale(gradOld, -lambda, updateForce);
        vectorAdd(newPosition, updateForce, newPosition);
        sigma = DCconstraintFunction(z0, a, b, newPosition);
        /*if (fabs(fabs(oldSigma)-fabs(sigma)) < 1e-15)
            break;*/
        oldSigma=sigma;
    }
    vectorCopy(newPosition,p->position);
}

int DCtangent1(double z0, double a, double b, const double *x, double *t1) {
    
    double rho = sqrt(x[0]*x[0] + x[1]*x[1]);
    double u1=b*b+2*x[2]*(-x[2]+z0);
    double u2=b*b*b*b+4*a*a*(x[2]-z0)*(x[2]-z0)+4*b*b*x[2]*(-x[2]+z0);
    if (fabs(rho) < MACHINE_EPSILON) {
        // special case at poles
        t1[0] = 1;
        t1[1] = 0;
        t1[2] = 0;
    } else {
        t1[0] = -x[1]*sqrt(u1+sqrt(u2))/(sqrt(2)*rho);
        t1[1] = x[0]*sqrt(u1+sqrt(u2))/(sqrt(2)*rho);
        t1[2] = 0;
        double normalization=vectorNorm(t1);
        t1[0]=t1[0]/normalization;
        t1[1]=t1[1]/normalization;
        t1[2]=t1[2]/normalization;
    }
    
    return TRUE;
}

int DCtangent2(double z0, double a, double b, const double *x, double *t2){
    
    double rho = sqrt(x[0]*x[0] + x[1]*x[1]);
    double u1=-2*x[2]+z0;
    double u2=2*a*a*(x[2]-z0)+b*b*(-2*x[2]+z0);
    double u3=b*b*b*b+4*a*a*(x[2]-z0)*(x[2]-z0)+4*b*b*x[2]*(-x[2]+z0);
    double u4=b*b+2*x[2]*(-x[2]+z0);
    
    if (fabs(rho) < MACHINE_EPSILON) {
        // special case at poles
        t2[0] = 0;
        t2[1] = 1;
        t2[2] = 0;
    } else {
        t2[0]=a*x[0]*(u1+u2/sqrt(u3))/(sqrt(2)*rho*sqrt(u4+sqrt(u3)));
        t2[1]=a*x[1]*(u1+u2/sqrt(u3))/(sqrt(2)*rho*sqrt(u4+sqrt(u3)));
        t2[2]=a;
        double normalization=vectorNorm(t2);
        t2[0]=t2[0]/normalization;
        t2[1]=t2[1]/normalization;
        t2[2]=t2[2]/normalization;
    }
    
    return TRUE;
}


//VectorFromSurface isn't defined right now, cause it only be used concerning harmonic surface energy and force.


double DCconformalFactor(double z0, double a, double b, double u) {
    double u1=b*b+2*a*u*(-a*u+z0);
    double u2=b*b*b*b+4*a*b*b*u*(-a*u+z0)+4*a*a*(-a*u+z0)*(-a*u+z0);
    double u3=2*a*a*a*u-2*a*b*b*u-2*a*a*z0+b*b*z0;
    
    return (double)1/2*sqrt(a*a*(u1+sqrt(u2))*(2+pow((-2*a*u+z0+u3/sqrt(u2)),2)/(u1+sqrt(u2))));
}

double DCconformalInvert(double z0, double a, double b, double u2) {
    double du=0.00001;
    double sum=0;
    double intMax=0;
    double result=1.0; // will return 1 in case the integral never reaches its max for numerical reasons
    
    if (u2<0-MACHINE_EPSILON || u2>1+MACHINE_EPSILON) {
        printf("conformal coordinate out of range - must be between 0 and 1\n");
        return(0);
    }
        
   
    for (double i=0+du/2; i<1.0; i+=du) {
        intMax += DCconformalFactor(z0,a,b,i)*du;
    }
    for (double i=0+du/2; i<1.0; i+=du) {
        sum += DCconformalFactor(z0,a,b,i)*du;
        if (sum > u2*intMax) {
            result = i;
            break;
        }
    }
    
    
    return(result);
    
}

void SetupPartitions(double rad, double a, double b, double z0, int nPart[], double pPosArray[]){
    
/*    if (*pp) {
        free(*pp);
    }
    if (*pCount) {
        free(*pCount);
    }
    
    //pp = arrayAlloc1dParticlePointer(pp, nPart[0]*nPart[1]*nPart[2]*np);
    *pp=malloc(sizeof(particle*)*nPart[0]*nPart[1]*nPart[2]*np);
    for (int i=0; i<nPart[0]*nPart[1]*nPart[2]*np;i++) {
        (*pp)[i]=NULL;
    }
    //pCount = arrayAlloc1dInt(pCount, nPart[0]*nPart[1]*nPart[2]);
    *pCount=malloc(sizeof(int)*nPart[0]*nPart[1]*nPart[2]);
    for (int i=0; i<nPart[0]*nPart[1]*nPart[2]; i++) {
        (*pCount)[i]=0;
    }*/
    
    for(int i=0; i<nPart[0]; i++){
        pPosArray[i] = (i+1) * 2 * (a + 2*rad) / nPart[0] - (a + 2*rad);
    }
    for(int i=0; i<nPart[1]; i++){
        pPosArray[nPart[0] + i] = (i+1) * 2 * (a + 2*rad) / nPart[1] - (a + 2*rad);
    }
    for(int i=0; i<nPart[2]; i++){
        pPosArray[nPart[0] + nPart[1] + i] = (i+1) * 2 * (a + 2*rad) / nPart[2] - (a + 2*rad);
    }


}//Merge functions of initialing partition, including setPartitionNumber, getMaxDimensions, allocatePP, allocate PCount and makePartitions. Keep partition constant during the relaxation.


//Other functions in partition.c don't need to change.


void DCgradDescentConstrainedIntegrate(double z0, double a, double b, particle *p, double dt) {
    double newPosition[3];
    double sigmaNew; // constraint function value
    double lambda; // Lagrange muliplier
    double gradOld[3]; // gradient of constraint function at old position
    double gradNew[3]; // gradient of constrain function at new position
    double updateForce[3];
    double sigmaOld;
    
    for (int i=0; i<3; i++) {
        newPosition[i] = p->position[i]+ (p->forceDet[i])*dt + (p->forceStoc[i])*sqrt(dt);
    }
    
    sigmaNew = DCconstraintFunction(z0, a, b, newPosition);
    sigmaOld=sigmaNew;
    
    do {
        DCconstraintGradient(z0,a,b,newPosition,gradNew);
        DCconstraintGradient(z0,a,b,p->position,gradOld);
        lambda = sigmaNew/(vectorDotProduct(gradNew, gradOld));
        vectorScale(gradOld, -lambda, updateForce);
        vectorAdd(newPosition, updateForce, newPosition);
        sigmaNew = DCconstraintFunction(z0, a, b, newPosition);
        /*if (fabs(fabs(sigmaNew)-fabs(sigmaOld)) < 1e-15)
            break;*/
        sigmaOld=sigmaNew;
    } while (fabs(sigmaNew) > MACHINE_EPSILON);
    vectorCopy(newPosition,p->position);
    
}

void DCundoOverlaps(double z0, double a, double b, double maxRadius, particle p[],  int nPart[], int pCount[], particle *pp[], double pPosArray[],
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
            DCgradDescentConstrainedIntegrate(z0,a,b,&p[k],*dtOverlap);
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


/*For the given particle, take a random step in the tangent plan and project to the suface.*/
int DCstep(double z0, double a, double b, double delta, double *xin, double *xout) {
    
    double d1;
    double d2;
    double t1[3];
    double t2[3];
    double t1Scaled[3];
    double t2Scaled[3];
    double stepVector[3];
    double newPosition[3];
    
    double sigma;
    double oldSigma;
    double lambda;
    double gradOld[3];
    double gradNew[3];
    double updateForce[3];
    
    // unit tangent vectors
    DCtangent1(z0, a, b, xin, t1);
    DCtangent2(z0, a, b, xin, t2);
    
    // step size along each tangent
    d1 = delta*randNormal();
    d2 = delta*randNormal();
    
    // scale tangent vectors accordingly
    vectorScale(t1, d1, t1Scaled);
    vectorScale(t2, d2, t2Scaled);
    
    // compose step vector
    vectorAdd(t1Scaled, t2Scaled, stepVector);
    
    // take step
    vectorAdd(xin, stepVector, newPosition);
    
    DCconstraintGradient(z0, a, b, xin, gradOld);
    sigma=DCconstraintFunction(z0, a, b, newPosition);
    oldSigma=sigma;
    while(fabs(sigma)>MACHINE_EPSILON){
        DCconstraintGradient(z0, a, b, newPosition, gradNew);
        lambda=sigma/vectorDotProduct(gradNew, gradOld);
        vectorScale(gradOld, -lambda, updateForce);
        vectorAdd(newPosition, updateForce, newPosition);
        sigma=DCconstraintFunction(z0, a, b, newPosition);
        /*if (fabs(fabs(oldSigma)-fabs(sigma)) < 1e-15)
            break;*/
        oldSigma=sigma;
    }
    vectorCopy(newPosition, xout);
    return TRUE;
}

//addStochasticForce, addHarmonicSurfaceForce, harmonicSurfaceEnergy should also be changed. But I don't need it for now.
void DCsetPartitionNumber(int np, double rad, double z0, double a, double b, int nPart[], double pFactor) {
    for (int i=0; i<3; i++) {
        if( (int) a/(pFactor*rad) > 1 ){
            nPart[i] = (int) a/(pFactor*rad);
        } else {
            nPart[i] = 1;
        }
    }
    nPart[3] = np;
}

particle** DCallocatePP(int np, double rad, double z0, double a, double b, int nPart[], particle *pp[], double pFactor) {
    if (pp) {
        free(pp);
    }
    pp = arrayAlloc1dParticlePointer(pp, nPart[0]*nPart[1]*nPart[2]*np); //we can now allocate less than this
    for (int i=0; i<nPart[0]*nPart[1]*nPart[2]*np; i++) {
        pp[i]=NULL;
    }
    
    return pp;
}

int* DCallocatePCount(double rad, double z0, double a, double b, int nPart[], int pCount[], double pFactor) {
    if (pCount) {
        free(pCount);
    }
    pCount = arrayAlloc1dInt(pCount, nPart[0]*nPart[1]*nPart[2]); //we can now allocate less than this
    for (int i=0; i<nPart[0]*nPart[1]*nPart[2]; i++) {
        pCount[i]=0;
    }
    
    return pCount;
}

void DCmakePartitions(double z0, double a, double b, double maxRadius, int nPart[], double pPosArray[] ){
    //vary these functions for non-uniform partitioning
    for(int i=0; i<nPart[0]; i++){
        pPosArray[i] = (i+1) * 2 * (a + 2*maxRadius) / nPart[0] - (a + 2*maxRadius);
    }
    for(int i=0; i<nPart[1]; i++){
        pPosArray[nPart[0] + i] = (i+1) * 2 * (a + 2*maxRadius) / nPart[1] - (a + 2*maxRadius);
    }
    for(int i=0; i<nPart[2]; i++){
        pPosArray[nPart[0] + nPart[1] + i] = (i+1) * 2 * (a + 2*maxRadius) / nPart[2] - (a + 2*maxRadius);
    }
    
}

//integrand of volume
double rhoSquare(double z0, double a, double b, double z){
    return 1.0/2*(b*b+2*z*(-z+z0)+sqrt(b*b*b*b+4*a*a*(z-z0)*(z-z0)+4*b*b*z*(-z+z0)));
}

//Romberg Integration to calculate the volume.
double volumeDC (double z0, double a, double b){
    double R[20] = {0};
    unsigned int m = 0;
    double h = a;
    R[0] = 1.0/2*h*(rhoSquare(z0, a, b, 0)+rhoSquare(z0, a, b, a));
    h = h/2;
    m++;
    R[1] = 1.0/2*R[0]+h*rhoSquare(z0, a, b, a/2);
    R[0] = 4.0/3*R[1]-1.0/3*R[0];
    while(fabs(R[0]-R[1]) > MACHINE_EPSILON){
        m++;
        if(m>=20)
            return -1;
        h = h/2;
        R[m] = 1.0/2*R[m-1];
        for (int i=1; i<=(int)pow(2, m-1); i++)
            R[m] = R[m]+h*rhoSquare(z0, a, b, (2*i-1)*h);
        for (int j=m-1, i=1; j>=0; j--,i++)
            R[j] = 1.0/(pow(4,i)-1)*(pow(4,i)*R[j+1]-R[j]);
    };
    //printf("iteration: %d\n",m);
    return 2*PI*R[0];
}

void findNewParameters(double *a, double *b, double volume, double lambda){
    double initialVolume = volumeDC(0, 1, lambda);
    *a = pow(volume/initialVolume,1.0/3);
    *b = lambda * (*a);
}

double findA ( double b, double volume, double machep, double t)
// findA seeks the a in an interval [0,1] to match the volume, given the value of b.

/*Parameters:
 
 Input, double b, the value of b.
 
 Input, double MACHEP, an estimate for the relative machine
 precision.
 
 Input, double T, a positive error tolerance.
 
 Input, double volume, target volume.
 
 Output, the estimated value of a that matches the volue*/
{
    double c;
    double d;
    double e;
    double fa;
    double fb;
    double fc;
    double m;
    double p;
    double q;
    double r;
    double s;
    double sa;
    double sb;
    double tol;
    /*
     Make local copies of A and B.
     */
    sa = 0;
    sb = 1;
    fa = volumeDC(0, sa, b)-volume;
    fb = volumeDC(0, sb, b)-volume;
    
    c = sa;
    fc = fa;
    e = sb - sa;
    d = e;
    
    for ( ; ; )
    {
        if ( fabs ( fc ) < fabs ( fb ) )
        {
            sa = sb;
            sb = c;
            c = sa;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        
        tol = 2.0 * machep * fabs ( sb ) + t;
        m = 0.5 * ( c - sb );
        
        if ( fabs ( m ) <= tol || fb == 0.0 )
        {
            break;
        }
        
        if ( fabs ( e ) < tol || fabs ( fa ) <= fabs ( fb ) )
        {
            e = m;
            d = e;
        }
        else
        {
            s = fb / fa;
            
            if ( sa == c )
            {
                p = 2.0 * m * s;
                q = 1.0 - s;
            }
            else
            {
                q = fa / fc;
                r = fb / fc;
                p = s * ( 2.0 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0 ) );
                q = ( q - 1.0 ) * ( r - 1.0 ) * ( s - 1.0 );
            }
            
            if ( 0.0 < p )
            {
                q = - q;
            }
            else
            {
                p = - p;
            }
            
            s = e;
            e = d;
            
            if ( 2.0 * p < 3.0 * m * q - fabs ( tol * q ) &&
                p < fabs ( 0.5 * s * q ) )
            {
                d = p / q;
            }
            else
            {
                e = m;
                d = e;
            }
        }
        sa = sb;
        fa = fb;
        
        if ( tol < fabs ( d ) )
        {
            sb = sb + d;
        }
        else if ( 0.0 < m )
        {
            sb = sb + tol;
        }
        else
        {
            sb = sb - tol;
        }
        
        fb = volumeDC (0, sb, b)-volume;
        
        if ( ( 0.0 < fb && 0.0 < fc ) || ( fb <= 0.0 && fc <= 0.0 ) )
        {
            c = sa;
            fc = fa;
            e = sb - sa;
            d = e;
        }
    }
    return sb;
}
