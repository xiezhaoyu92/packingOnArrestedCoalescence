//
//  surface.c
//  ellipsoidRelaxation
//
//  Created by Chris Burke on 2/28/16.
//  Copyright © 2016 tufts. All rights reserved.
//

#include <complex.h>
#include <float.h>
#include <math.h>
#include "surface.h"
#include "misc.h"
#include "mt19937ar.h"



/**********************
 * surface projection *
 **********************/

/* Projects a point onto the surface of an ellipsoid
 * In      - a       - semi-major/minor axis of the ellipsoid
 *           b       -    "                            "
 *           rho     - rho component of position in cylindrical coordinates
 *           z       - z component of position in cylindrical coordinates
 * Out     - rhoNew  - rho component of the new point on the surface - can be the same as rho0
 *         - zNew    - z component of the new point on the surface - can be the same as z0
 * Returns - none
 */
int closestEllipsePoint(double a, double b, double rho, double z, double *rhoNew, double *zNew)
{
    
    double complex aC;
    double complex bC;
    double complex rhoC;
    double complex zC;
    double complex rhoNewC;
    double complex zNewC;
    int zSign;
    int flipped;
    
    // quartic root expressions
    double complex c0;
    double complex c1;
    double complex c2;
    double complex d0;
    double complex d1;
    double complex d2;
    double complex d3;
    double complex d4;
    double complex p;
    double complex q;
    double complex delta0;
    double complex delta1;
    double complex S;
    double complex Q;
    
    double flipCutoff = 1e-5; // cutoff for the distance from the z axis, below which the axes are flipped
    
    
    if (fabs(rho) < flipCutoff) {
        flipped = TRUE;
        aC = b;
        bC = a;
        rhoC = fabs(z);
        zC = rho;
    } else {
        flipped = FALSE;
        aC = a;
        bC = b;
        rhoC = rho;
        zC = fabs(z);
    }
    
    if (fabs(z) < MACHINE_EPSILON) {
        zSign = randInteger(2)*2-3; // from 2 or 4 to -1 or 1
        //printf("zSign: %i\n", zSign);
    } else if (z > 0) {
        zSign = 1;
    } else {
        zSign = -1;
    }
    
    c2 = aC*aC - bC*bC;
    c0 = bC*bC*zC;
    c1 = -aC*aC*rhoC;
    d4 = bC*bC*c2*c2;
    d3 = 2*bC*bC*c1*c2;
    d2 = bC*bC*c1*c1 + aC*aC*(c0*c0 - bC*bC*c2*c2);
    d1 = -2*aC*aC*bC*bC*c1*c2;
    d0 = -aC*aC*bC*bC*c1*c1;
    p = (-3*d3*d3 + 8*d2*d4)/(8.*d4*d4);
    q = (d3*d3*d3 - 4*d2*d3*d4 + 8*d1*d4*d4)/(8.*d4*d4*d4);
    delta0 = d2*d2 - 3*d1*d3 + 12*d0*d4;
    delta1 = 2*d2*d2*d2 - 9*d1*d2*d3 + 27*d0*d3*d3 + 27*d1*d1*d4 - 72*d0*d2*d4;
    Q = cpow((delta1 + csqrt(-4*delta0*delta0*delta0 + delta1*delta1))/2.,0.3333333333333333);
    S = 0.5*csqrt(-((2*p)/3.) + (delta0/Q + Q)/(3.*d4));
    if ( ((a > b) && flipped) || ((a < b) && !flipped)){
        rhoNewC = -(d3/(4.*d4)) + S + 0.5*csqrt(-2*p - q/S - 4*S*S); // "x3" root in mathematica notebook
    } else {
        rhoNewC = -(d3/(4.*d4)) - S + 0.5*csqrt(-2*p + q/S - 4*S*S); // "x1" root in mathematica notebook
    }
    zNewC = bC*csqrt(1-rhoNewC*rhoNewC/(aC*aC));
    
    if (flipped) {
        *rhoNew = creal(zNewC);
        *zNew = zSign*creal(rhoNewC);
    } else {
        //printf("z: %f+i*%f\n",creal(zNewC),cimag(zNewC));
        //printf("rho: %f+i*%f\n",creal(rhoNewC),cimag(rhoNewC));
        *rhoNew = creal(rhoNewC);
        *zNew = zSign*creal(zNewC);
    }
    //printf("rho, z: %f %f\n", rho, z);
    //printf("rhoNew, zNew: %f %f\n", *rhoNew, *zNew);
    
    return TRUE;
}

/* Projects a point onto the surface of an circle
 * In      - a      - radius of the circle
 *           rho    - rho component of position in cylindrical coordinates
 *           z      - z component of position in cylindrical coordinates
 * Out     - rhoNew  - rho component of the new point on the surface - can be the same as rho0
 *         - zNew    - z component of the new point on the surface - can be the same as z0
 * Returns - none
 */
int closestCirclePoint(double a, double rho, double z, double *rhoNew, double *zNew)
{
    double r;
    
    r = sqrt(rho*rho + z*z);
    
    *rhoNew = rho*a/r;
    *zNew = z*a/r;
    
    return TRUE;
}


/* Projects a point onto the surface of an ellipsoid
 * In      - a      - semi-major/minor axis of the ellipsoid
 *           b      -    "                            "
 *           xin    - a point to project
 * Out     - xout   - the projected point (this can be the same as xin)
 * Returns - TRUE on success, FALSE otherwise
 */
int project(double a, double b, const double *xin, double *xout) {
    
    double phi;
    double rho = sqrt(xin[0]*xin[0]+xin[1]*xin[1]);
    double z = xin[2];
    
    // if on the z axis, choose a random phi
    if ( (fabs(xin[0]) < MACHINE_EPSILON) && (fabs(xin[1]) < MACHINE_EPSILON) ) {
        phi = genrand_real2()*2*PI;
        //printf("phi: %f\n", phi);
    } else {
        phi = atan2(xin[1], xin[0]);
    }
    
    
    //printf("%f %f %f\n", phi, rho, z);
    if (fabs(a-b) < MACHINE_EPSILON) {
        // sphere
        if ( fabs(rho) < MACHINE_EPSILON && fabs(z) < MACHINE_EPSILON ) {
            // if at origin, move to random point (all points have equal probability)
            z = a*(genrand_real1()*2 - 1);
            rho = sqrt(a*a-z*z);
        } else {
            closestCirclePoint(a, rho, z, &rho, &z);
        }
    } else if (a > b) {
        // oblate
        if ( fabs(rho) < MACHINE_EPSILON && fabs(z) < MACHINE_EPSILON ) {
            // if at origin, move to top or bottom pole
            z = b*(randInteger(2)*2-3); // b or -b
            rho = 0;
        } else {
            closestEllipsePoint(a, b, rho, z, &rho, &z);
        }
    } else {
        // prolate
        if ( fabs(rho) < MACHINE_EPSILON && fabs(z) < MACHINE_EPSILON ) {
            // if at origin, move to random point on equator
            z = 0;
            rho = a;
        } else if (fabs(rho) < MACHINE_EPSILON && fabs(z) >= b) {
            // if on z axis outside of surface, go to closest pole
            if (z > 0) {
                z = b;
            } else {
                z = -b;
            }
        } else {
            closestEllipsePoint(a, b, rho, z, &rho, &z);
        }
    }
    //printf("phi, rho, z: %f %f %f\n", phi, rho, z);
    xout[0] = cos(phi)*rho;
    xout[1] = sin(phi)*rho;
    xout[2] = z;
    //printf("x, y, z:  %f %f %f\n", xout[0], xout[1], xout[2]);
    
    return TRUE;
}

#ifdef DropletCoalescense

#endif

/* moves a point with surface of an ellipsoid - scale particle coordinates with ellispoid axes
 * In      - a1     - initial semi-major/minor axis of the ellipsoid
 *           b1     -    "                            "
 *         - a2     - final semi-major/minor axis of the ellipsoid
 *           b2     -    "                            "
 *           xin    - a point to move
 * Out     - xout   - the moved point (this can be the same as xin)
 * Returns - TRUE on success, FALSE otherwise
 */
int moveWithSurface(double a1, double b1, double a2, double b2, const double *xin, double *xout) {
    xout[0] = xin[0]*a2/a1;
    xout[1] = xin[1]*a2/a1;
    xout[2] = xin[2]*b2/b1;
    
    return TRUE;
}


/* Returns the tangent vector along the phi direction at a given point
 * In      - a      - semi-major/minor axis of the ellipsoid
 *           b      -    "                            "
 *           xin    - a point to project
 * Out     - t1     - the tangent vector at xin)
 * Returns - TRUE on success, FALSE otherwise
 */
int tangent1(double a, double b, const double *xin, double *t1) {
    
    double rho = sqrt(xin[0]*xin[0] + xin[1]*xin[1]);
    
    if (fabs(rho) < MACHINE_EPSILON) {
        // special case at poles
        t1[0] = 1;
        t1[1] = 0;
        t1[2] = 0;
    } else {
        t1[0] = xin[1]/rho;
        t1[1] = -xin[0]/rho;
        t1[2] = 0;
    }
    
    return TRUE;
}

/* Returns the tangent vector orthogonal to the phi direction at a given point
 * In      - a      - semi-major/minor axis of the ellipsoid
 *           b      -    "                            "
 *           xin    - a point to project
 * Out     - t2     - the tangent vector at xin)
 * Returns - TRUE on success, FALSE otherwise
 */
int tangent2(double a, double b, const double *xin, double *t2) {
    
    double rho = sqrt(xin[0]*xin[0] + xin[1]*xin[1]);
    double sigma = sqrt( (xin[0]*xin[0] + xin[1]*xin[1])/(a*a*a*a) + xin[2]*xin[2]/(b*b*b*b) );
    
    if (fabs(rho) < MACHINE_EPSILON) {
        // special case at poles
        t2[0] = 0;
        t2[1] = 1;
        t2[2] = 0;
    } else {
        t2[0] = -(xin[0]*xin[2])/(b*b*rho*sigma);
        t2[1] = -(xin[1]*xin[2])/(b*b*rho*sigma);
        t2[2] = rho/(a*a*sigma);
    }
    
    return TRUE;
}

/*********************
 * surface evolution *
 *********************/




/* Calculates the surface area of an ellipsoid given the semi-major/minor axis lengths
 * In      - a      - semi-major/minor axis of the ellipsoid
 *           b      -    "                            "
 * Out     - none
 * Returns - surface area
 */
double ellipsoidSurfaceArea(double a, double b) {
    double e;
    if (fabs(a-b) < MACHINE_EPSILON) {
        // sphere
        return 4*PI*a*a;
    } else if (b>a) {
        // prolate
        e = sqrt(b*b-a*a)/b;
        return 2*PI*a*(a+b/e*asin(e));
    } else {
        // oblate
        e = sqrt(a*a-b*b)/a;
        return PI*(2*a*a+b*b/e*log((1+e)/(1-e)));
    }
}


/* Give the surface area at a given time such that the surface area
 * decays exponentially. At t=1, the area will have decreased by a fraction f
 * of the total ammount it would have decreased by at t=infinity
 * Input:	a0      - inital semi-major/minor axis length
 *          b0      - ''
 *          f       - fraction of the t=inf decay that will occur by t=1
 *          t       - time, between 0 and 1
 * Output:  Nothing
 * Returns:	surface area at the given time
 */
double areaAtTime(double a0, double b0, double f, double t) {
    double area0 = ellipsoidSurfaceArea(a0,b0);
    double area1 = surfaceAreaAtEvoParameter(a0,b0,1.0);
    double areaInf = area0-(area0-area1)/f;
    double tau = 1/log(1/(1-f));
    
    return areaInf+(area0-areaInf)*exp(-t/tau);
    
}


/* Give the time at which the surface has a given surface area. Assume surface area
 * decays exponentially. At t=1, the area will have decreased by a fraction f
 * of the total ammount it would have decreased by at t=infinity
 * Input:	a0      - inital semi-major/minor axis length
 *          b0      - ''
 *          f       - fraction of the t=inf decay that will occur by t=1
 *          area    - area at desired time
 * Output:  Nothing
 * Returns:	time at which the given surface area is reached
 */
double timeAtArea(double a0, double b0, double f, double area) {
    double area0 = ellipsoidSurfaceArea(a0,b0);
    double area1 = surfaceAreaAtEvoParameter(a0,b0,1.0);
    double areaInf = area0-(area0-area1)/f;
    double tau = 1/log(1/(1-f));
    
    return -tau*log((area-areaInf)/(area0-areaInf));
    
}


/* Give the surface area for a shape evolution parameter s given an initial shape
 * Input:	a0      - inital semi-major/minor axis length
 *          b0      - ''
 *          s       - shape evolution parameter, between 0 and 1
 * Output:  Nothing
 * Returns:	surface area at the given evolution parameter
 */
double surfaceAreaAtEvoParameter(double a0, double b0, double s) {
    double a;
    double b;
    ellipsoidUpdate(s, a0, b0, &a, &b);
    double area = ellipsoidSurfaceArea(a,b);
    
    return area;
}

/* Given a desired area and an initial shape, return the shape evolution parameter
 * at which that area occurs. Uses bisection method to invert surfaceAreaAtEvoParameter.
 * Input:	a0      - inital semi-major/minor axis length
 *          b0      - ''
 *          area    - the surface area corresponding to the desired evolution parameter
 * Output:  Nothing
 * Returns:	evolution parameter corresponding to the given area
 */
double evoParameterAtSurfaceArea(double a0, double b0, double area) {
    double eps = 1e-14; // precision goal
    double x1 = 0.0; //brackets --- we know the function is monotonically decreasing between these values
    double x2 = 1.0;
    int maxIter = 50; // maximum iterations
    double dx;
    double xMid;
    double root;
    
    double f = surfaceAreaAtEvoParameter(a0,b0,x1)-area;
    double fMid = surfaceAreaAtEvoParameter(a0,b0,x2)-area;
    
    
    
    if (f*fMid>0.0) {
        printf("not bracketed correctly: %f %f\n",f,fMid);
    }
    
    //orient search from lower f value
    if (f<=0.0) {
        root = x1;
        dx = x2-x1;
    } else {
        root = x2;
        dx = x1-x2;
    }
    
    for (int i=0; i<maxIter; i++) {
        dx*=0.5;
        xMid = root+dx;
        fMid = surfaceAreaAtEvoParameter(a0,b0,xMid)-area;
        if (fMid<=0.0) {
            root = xMid;
        }
        if (fabs(dx)<eps || fMid==0.0) {
            break;
        }
        if (i==maxIter-1) {
            printf("Too many bisections. dx: %f\n",dx);
        }
    }
    
    return root;
}


/* Give the ellipsoid axis parameters for a given shape evolution parameter
 * Input:	s		- evolution parameter, runs from 0 to 1
 *          a0      - inital semi-major/minor axis length
 *          b0      - ''
 * Output:  a		- semi-major/minor axis length at evolution parameter s
 *          b       - ''
 * Returns:	Nothing
 */
void ellipsoidUpdate(double s, double a0, double b0, double *a, double *b) {
    *a = a0+(pow(a0*a0*b0,1.0/3.0)-a0)*s;
    *b = b0*(a0*a0)/((*a)*(*a));
}


#ifdef  DEBUG0
/* return the constraint function value at a given position. equals 0 when constraint is satisfied
 * In      - a      - ellipsoid semi major/minor axis length
 *         - b      - ellipsoid semi major/minor axis length
 *         - x      - position (3d vector)
 * Out     - none
 * Returns - constraint function value
 */
double constraintFunction(double a, double b, double *x) {
    return (x[0]/a)*(x[0]/a)+(x[1]/a)*(x[1]/a)+(x[2]/b)*(x[2]/b)-1;
}
#endif

/* return the gradient of the constraint function at a given position
 * In      - a      - ellipsoid semi major/minor axis length
 *         - b      - ellipsoid semi major/minor axis length
 *         - x      - position (length 3 array)
 * Out     - grad   - constraint function gradient
 * Returns - none
 */
void constraintGradient(double a, double b, double *x, double *grad) {
    grad[0]=2*x[0]/(a*a);
    grad[1]=2*x[1]/(a*a);
    grad[2]=2*x[2]/(b*b);
}


/* Returns the distance vector to a given point from the surface.
 * First order approximation, valid for points near the surface.
 * In      - a      - ellipsoid semi major/minor axis length
 *         - b      - ellipsoid semi major/minor axis length
 *         - x      - position (length 3 array)
 * Out     - dx     - distance vector from nearest point on the surface
 * Returns - none
 */
void vectorFromSurface(double a, double b, double *x, double *dx) {
    double grad[3];
    constraintGradient(a, b, x, grad);
    double c = constraintFunction(a, b, x);
    
    // derived from delta c = grad c . delta x
    vectorScale(grad, c/vectorDotProduct(grad, grad), dx);
}


/* Calculate the conformal factor at a point on the surface
 * In      - a      - semi-major/minor axis of the ellipsoid
 *           b      -    "                            "
 *           u      - scaled z-axis coordinate
 * Out     - none
 * Returns - conformal factor
 */
double conformalFactor(double a, double b, double u) {
    return sqrt(a*a*(b*b+(a*a-b*b)*u*u));
}

/* Given a conformal coordinate, give the corresponding scaled z-axis coordinate.
 * Used for generating an even surface distribution.
 * This is done by inverting an integral equation numerically - this is good enough
 * to generate an evenly distributed starting configuration, which will then be allowed to diffuse,
 * but does not have enough precision for general usage.
 * Includes exact solution for the special case of a sphere.
 * In      - a      - semi-major/minor axis of the ellipsoid
 *           b      -    "                            "
 *           u2     - conformal coordinate - must be in range [0,1]
 * Out     - none
 * Returns - scaled z-axis coordinate
 */
double conformalInvert(double a, double b, double u2) {
    double du=0.00001;
    double sum=0;
    double intMax=0;
    double result=1.0; // will return 1 in case the integral never reaches its max for numerical reasons
    
    if (u2<0-MACHINE_EPSILON || u2>1+MACHINE_EPSILON) {
        printf("conformal coordinate out of range - must be between 0 and 1\n");
        return(0);
    }
    
    
    if (fabs(a-b) < MACHINE_EPSILON) { //special case: sphere
        result = 2*u2-1;
    } else { // general ellipsoid case
        for (double i=-1.0+du/2; i<1.0; i+=du) {
            intMax += conformalFactor(a,b,i)*du;
        }
        for (double i=-1.0+du/2; i<1.0; i+=du) {
            sum += conformalFactor(a,b,i)*du;
            if (sum > u2*intMax) {
                result = i;
                break;
            }
        }
    }
    
    
    return(result);
    
}




/*Zhaoyu's Relaxation*/

void gradDescentConstrainedIntegrateProject(double aOld, double bOld, double a, double b, particle *p) {
    double newPosition[3];
    double sigmaNew; // constraint function value
    double lambda; // Lagrange muliplier
    double gradOld[3]; // gradient of constraint function at old position
    double gradNew[3]; // gradient of constrain function at new position
    double updateForce[3];
    double normalization;
    double normalDirection[3];
    constraintGradient(aOld, bOld, p->position, normalDirection);
    normalization=vectorNorm(normalDirection);
    for (int i=0; i<3; i++) {
        newPosition[i] = p->position[i]+ 0.05*normalDirection[i]/normalization;
    }
    
    sigmaNew = constraintFunction(a, b, newPosition);
    constraintGradient(aOld, bOld, p->position, gradOld);
    
    do {
        constraintGradient(a,b,newPosition,gradNew);
        lambda = sigmaNew/(vectorDotProduct(gradNew, gradOld));
        vectorScale(gradOld, -lambda, updateForce);
        vectorAdd(newPosition, updateForce, newPosition);
        sigmaNew = constraintFunction(a, b, newPosition);
    } while (fabs(sigmaNew) > MACHINE_EPSILON);
    vectorCopy(newPosition,p->position);
}
