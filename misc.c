//
//  misc.c
//  ellipsoidRelaxation
//
//  Created by Chris Burke on 2/28/16.
//  Copyright © 2016 tufts. All rights reserved.
//

#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "particles.h"
#include "misc.h"
#include "mt19937ar.h"

/***************
 * vector math *
 ***************/

/* Vector addition
 * Input:	v1		- A point.
 *			v2		- A second point
 * Output:  v3		- v1 + v2
 * Returns:	Nothing
 */
void vectorAdd(const double v1[3], const double v2[3], double *v3) {
    v3[0] = v1[0] + v2[0];
    v3[1] = v1[1] + v2[1];
    v3[2] = v1[2] + v2[2];
    return;
}

/* Vector subtraction
 * Input:	v1		- A point.
 *			v2		- A second point
 * Output:  v3		- v1 - v2
 * Returns:	Nothing
 */
void vectorSubtract(double v1[3],double v2[3],double *v3) {
    v3[0] = v1[0] - v2[0];
    v3[1] = v1[1] - v2[1];
    v3[2] = v1[2] - v2[2];
    return;
}

/* Vector scaling
 * Input:	v1		- A point.
 *			scale	- a scaling factor
 * Output:  v2		- scale*v1
 * Returns:	Nothing
 */
void vectorScale(const double v1[3], double scale, double v2[3]) {
    v2[0] = v1[0]*scale;
    v2[1] = v1[1]*scale;
    v2[2] = v1[2]*scale;
    return;
}

/* copy the value of a vector
 * Input:	x1		- the vector to be copied
 * Output:  x2		- array the vector value is copied to
 * Returns:	Nothing
 */
void vectorCopy(const double *x1, double *x2) {
    x2[0] = x1[0];
    x2[1] = x1[1];
    x2[2] = x1[2];
    return;
}

#ifdef  DEBUG0
/* Vector norm
 * Input:	v1		- A vector.
 * Output:  None
 * Returns:	norm(v1)
 */
double vectorNorm(double v1[3]) {
    return sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
}

/* Square of a vector's norm - to avoid taking square roots during comparisons
 * Input:	v1		- A vector.
 * Output:  None
 * Returns:	norm(v1)
 */
double vectorNormSquare(double v1[3]) {
    return v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
}

/* Vector dot product
 * Input:	v1		- A vector.
 *			v2		- A second vector.
 * Output:  None
 * Returns:	v1 . v2
 */
double vectorDotProduct(double v1[3],double v2[3]) {
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}
#endif



/*********************
 * memory allocation *
 *********************/

/* Allocates a 1D array
 * In      - x      - pointer of array to be allocated
 *           n      - size of array
 * Out     - none
 * Returns - x
 */
double* arrayAlloc1d(double *x, int n) {
    
    x=malloc(sizeof(double)*n);
    if (!x) {
        printf("array allocation failed");
    }
    
    return(x);
}

/* Allocates a 1D array of ints
 * In      - x      - pointer of array to be allocated
 *           n      - size of array
 * Out     - none
 * Returns - x
 */
int* arrayAlloc1dInt(int *x, int n) {
    
    x=malloc(sizeof(int)*n);
    if (!x) {
        printf("array allocation failed");
    }
    
    return(x);
}

/* Allocates a 1D array of particle pointers
 * In      - x      - pointer of array to be allocated
 *           n      - size of array
 * Out     - none
 * Returns - x
 */
particle** arrayAlloc1dParticlePointer(particle **x, int n) {
    
    x=malloc(sizeof(particle*)*n);
    if (!x) {
        printf("array allocation failed");
    }
    
    return(x);
}


/* Allocates a 2D array
 * In      - x      - pointer of array to be allocated
 *           n1     - size of each dimensions
 n2
 * Out     - none
 * Returns - x      - pointer of array to be allocated
 */
double** arrayAlloc2d(double **x, int n1, int n2) {
    int success=TRUE;
    
    x=malloc(sizeof(double *)*n1);
    if (x) {
        for(int i=0; i<n1; i++) x[i]=NULL;
        
        for(int i=0; i<n1; i++) {
            x[i]=malloc(sizeof(double)*n2);
            if (!x[i]) {
                success=FALSE;
                break;
            }
        }
        
        if (!success) {
            for (int i=0; i<n1; i++) if (x[i]) free(x[i]);
            free(x);
            printf("array allocation failed");
        }
    } else {
        printf("array allocation failed");
    }
    
    return(x);
}



/* a comparison function for sorting
 * In      - a, b       - pointers to doubles to be compared
 * Out     - none
 * Returns - negative value if *a < *b, 0 if equal, positive if *a > *b
 * From GNU documentation
 */
int compareDoubles(const void *a, const void *b)
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;
    
    return (*da > *db) - (*da < *db);
}



/***************************
 * random number functions *
 ***************************/



/* initialize random number generator
 * Input:       None
 * Output:      None
 * Returns:     None
 */
void randinitialize(void) {
    FILE *urandom;
    unsigned long seed=(unsigned long) time(NULL);
    char bytes[sizeof(unsigned long)];
    
    /* Initialize the random number generator */
    urandom=fopen("/dev/urandom", "r");
    if (urandom) {
        for(int i=0;i<sizeof(unsigned long);i++) bytes[i]=(char) fgetc(urandom);
        fclose(urandom);
        
        seed = *((unsigned long *) bytes);
    } else printf("Warning: initializing random number generator using time-not recommended for production runs.\n");
    
    init_genrand(seed);
}

/* Generate a random integer from 1 to max
 * Input:	max	      - largest possible number to return
 * Output:      None
 * Returns:	a random integer from 1 to max
 */
int randInteger(int max) {
    
    return (int) (genrand_real1()*max+1.0);
}

/* Generate a random number with a gaussian distribution with sigma=1
 * Input:	None
 * Output:      None
 * Returns:	a random double
 */
double randNormal(void) {
    double x,y,r;
    
    do {
        x=2.0*genrand_real1()-1.0;
        y=2.0*genrand_real1()-1.0;
        
        r=x*x+y*y;
    } while (r>=1.0);
    
    return x*sqrt((-2.0*log(r))/r);
}



