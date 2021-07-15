//
//  main.c
//  Droplet Coalescence
//
//  Created by Zhaoyu Xie on 8/1/16.
//  Copyright © 2016 Zhaoyu Xie. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "surface.h"
#include "DropletCoalescence.h"
#include "partition.h"
#include "particles.h"
#include "forces.h"
#include "misc.h"
#include "mt19937ar.h"
#include <getopt.h>


#ifndef TRUE
#define TRUE -1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef PI
#define PI 3.141592653589793
#endif

#ifndef MACHINE_EPSILON
#define MACHINE_EPSILON = 1e-14
#endif


int main(int argc, char* argv[]) {
    double z0=0;
    int np = 1000;
    double a0 = 1.0;
    double b0 = 0.1;
    double rad = 0.03;
    double pFactor = 3.0;
    int outputSteps = 10000;
    double ratio=1.0;//radius ratio of smaller particles to larger particles.
    double mixRatio=0.5;//quantity ratio of large particles.
    double scalingFactor = 0.5; // determine how b scales with time
    double dtDiffusion0 = 1;
    double dtOverlap0 = 1e-3;
    double dtOverlap = dtOverlap0;
    double dtOverlapStartup = 1e-3;
    double expDecayFrac = 0.9;
    double timestepScaleFactor = 100;
    
    int maxUndoStepsStart = 1e6;
    int startupDiffusionSteps = 20000;
    int maxUndoSteps = 1e4;
    int relaxationSteps = 1e5;
    double equilibrateTime = 1e6;
    
    double diffusionTimescale = 0.25; //default dtDiffusionScale, unless defined later.
    double deltaPhi = 1.0e-6;
    
    int readUnfinishedQ = FALSE;
    int animationQ = FALSE;
    int readJammedQ = FALSE;
    int dynamicsQ = FALSE;
    int moveWithSurfaceQ = FALSE; //if true, scaling particles to new surface for static case
    
    int animationRate = 100;

    FILE *animationFile=NULL;

    
    static struct option longOptions[]={
        {"readUnfinished",      no_argument,        NULL, 'u'},
        {"aspectRatio",         required_argument,  NULL, 'b'},
        {"particleNumber",      required_argument,  NULL, 'n'},
        {"radius",              required_argument,  NULL, 'r'},
        {"radRatio",            required_argument,  NULL, 't'},
        {"intermixRatio",       required_argument,  NULL, 'i'},
        {"animate",             no_argument,        NULL, 'm'},
        {"readJammed",          no_argument,        NULL, 'j'},
        {"dynamics",            no_argument,        NULL, 'd'}, //if set -d, no need to set -b
        {"dtDiffusionScale",    required_argument,  NULL, 's'},
        {"scalingFactor",       required_argument,  NULL, 'f'}, //only used for dynamic case
        {"moveWithSurface",     no_argument,        NULL, 'e'}, //only used for static case
        {0,                 0,                  0,     0 }
    };
    
    int optIndex = 0;
    int opt;
    while ((opt = getopt_long(argc, argv,"ub:n:r:t:i:mjds:f:e",
                              longOptions, &optIndex )) != -1) {
        switch (opt) {
            case 0:
                break;
            case 'b' :
                b0 = atof(optarg);
                break;
            case 'u' :
                readUnfinishedQ = TRUE;
                break;
            case 'n' :
                np = atoi(optarg);
                break;
            case 'r' :
                rad = atof(optarg);
                break;
            case 't' :
                ratio = atof(optarg);
                break;
            case 'i' :
                mixRatio = atof(optarg);
                break;
            case 'm' :
                animationQ = TRUE;
                break;
            case 'j' :
                readJammedQ = TRUE;
                break;
            case 'd' :
                dynamicsQ = TRUE;
                break;
            case 's':
                diffusionTimescale = atof(optarg);
                break;
            case 'f':
                scalingFactor = atof(optarg);
                break;
            case 'e':
                moveWithSurfaceQ = TRUE;
                break;
            default:
                exit(1);
        }
    }
    
    if (readUnfinishedQ && readJammedQ) {
        printf("readJammmed and readUnfinished are mutually exclusive options. Please specify only one.\n");
        exit(1);
    }
    if (readUnfinishedQ || readJammedQ) {
        printf("note: previous simulation state is being read in. Any new simulation parameters will be ignored\n");
    }

    if (readUnfinishedQ || readJammedQ) {
        FILE *npInFile=NULL;
        npInFile = fopen("npts.dat","r");
        if (npInFile) {
            fscanf(npInFile, "%i", &np);
            fclose(npInFile);
        } else {
            printf("npInFile pointer is null\n");
            exit(1);
        }
    }
    
    particle p[np];
    
    double simTimeStart = 0;
    int simStep = 0;
    double nextRelaxationTime = 0;
    double a;
    double b;
    double aOld;
    double bOld;
    
    if (readUnfinishedQ || readJammedQ) {
        FILE *radInFile=NULL;
        radInFile = fopen("radius.dat","r");
        if (radInFile) {
            fscanf(radInFile, "%lf", &rad);
            fclose(radInFile);
        } else {
            printf("radInFile pointer is null\n");
            exit(1);
        }
        FILE *radiiInFile=NULL;
        radiiInFile = fopen("radii.dat","r");
        if (radiiInFile) {
            for (int i=0; i<np; i++) {
                fscanf(radiiInFile, "%lf", &p[i].rad);
            }
            fclose(radiiInFile);
        } else {
            printf("radiiInFile pointer is null\n");
            exit(1);
        }
        FILE *arInFile=NULL;
        arInFile = fopen("DCParamsInit.dat","r");
        if (arInFile) {
            fscanf(arInFile, "%lf %lf %lf", &z0, &a0, &b0);
            fclose(arInFile);
        } else {
            printf("arInFile pointer is null\n");
            exit(1);
        }

        
        FILE *diffTimeFileIn=NULL;
        diffTimeFileIn = fopen("diffusionTimescale.dat","r");
        if (diffTimeFileIn) {
            fscanf(diffTimeFileIn, "%lf", &diffusionTimescale);
            fclose(diffTimeFileIn);
        } else {
            printf("diffTimeFileIn pointer is null\n");
            exit(1);
        }
        
        
        FILE *scalingFactorFileIn=NULL;
        scalingFactorFileIn = fopen("scalingFactor.dat","r");
        if (scalingFactorFileIn) {
            fscanf(scalingFactorFileIn, "%lf", &scalingFactor);
            fclose(scalingFactorFileIn);
        } else {
            printf("scalingFactorFileIn pointer is null\n");
            exit(1);
        }

        if (readUnfinishedQ){
            
            FILE *timeInFile=NULL;
            timeInFile = fopen("simTimeTemp.dat","r");
            if (timeInFile) {
                fscanf(timeInFile, "%lf", &simTimeStart);
                fclose(timeInFile);
            } else {
                printf("timeInFile pointer is null\n");
                exit(1);
            }
            
            FILE *nextRelaxationFile=NULL;
            nextRelaxationFile = fopen("nextRelaxationTime.dat","r");
            if (nextRelaxationFile) {
                fscanf(nextRelaxationFile, "%lf", &nextRelaxationTime);
                fclose(nextRelaxationFile);
            } else {
                printf("nextRelaxationFile pointer is null\n");
                exit(1);
            }
            
            FILE *arInFile=NULL;
            arInFile = fopen("ParamsTemp.dat","r");
            if (arInFile) {
                fscanf(arInFile, "%lf %lf %lf", &z0, &a, &b);
                fclose(arInFile);
            } else {
                printf("arInFile pointer is null\n");
                exit(1);
            }
            
            FILE *configInFile=NULL;
            configInFile = fopen("configurationTemp.asc","r");
            if (configInFile) {
                for (int n=0; n<np; n++) {
                    int pCount = fscanf(configInFile, "%lf %lf %lf", &p[n].position[0], &p[n].position[1], &p[n].position[2]);
                    if (pCount != 3) {
                        printf("not 3 particles\n");
                        exit(1);
                    }
                }
                fclose(configInFile);
            } else {
                printf("configInFile pointer is null\n");
                exit(1);
            }
        }
        
        if (readJammedQ) {
            
            FILE *timeInFile=NULL;
            timeInFile = fopen("simTimeArrested.dat","r");
            if (timeInFile) {
                fscanf(timeInFile, "%lf", &simTimeStart);
                fclose(timeInFile);
            } else {
                printf("timeInFile pointer is null\n");
                exit(1);
            }
            
            FILE *arInFile=NULL;
            arInFile = fopen("DCParams.dat","r");
            if (arInFile) {
                fscanf(arInFile, "%lf %lf %lf", &z0, &a, &b);
                fclose(arInFile);
            } else {
                printf("arInFile pointer is null\n");
                exit(1);
            }
            
            FILE *configInFile=NULL;
            configInFile = fopen("configuration.asc","r");
            if (configInFile) {
                for (int n=0; n<np; n++) {
                    int pCount = fscanf(configInFile, "%lf %lf %lf", &p[n].position[0], &p[n].position[1], &p[n].position[2]);
                    if (pCount != 3) {
                        printf("not 3 particles\n");
                        exit(1);
                    }
                }
                fclose(configInFile);
            } else {
                printf("configInFile pointer is null\n");
                exit(1);
            }
        }
        
        for (int i=0; i<np; i++) {
            vectorCopy(p[i].position, p[i].postPreviousRelaxPosition);
        }
    }
    
    if(!(readUnfinishedQ||readJammedQ)){
        a = a0;
        b = b0;
        if (dynamicsQ){  //To use dynamics, one must set b0=0 initially
            b0 = 0;
            b = b0;
        }
    }
    
    double volume = volumeDC(0, a0, b0);
    double startupDelta = 2*rad/100;
    double delta = 2.*rad/1000;
    double diffusionCoeff = delta*delta;
    double simTimeFinal = 2*rad*rad/(diffusionCoeff*diffusionTimescale);
    double scalingConst = pow(0.25, 1.0/3)*a0/pow(simTimeFinal,scalingFactor);// the scaling constant to let b evolve as t^scalingFactor
    
    double simTime;
    double dtRelaxation0;
    double dtRelaxation;
    double diffusionSteps;
    
    int nPart[4];
    particle **pp = NULL;
    int *pCount = NULL;
    
    for (int i=0; i<3; i++) {
        if( (int) a/(pFactor*rad) > 1 ){
            nPart[i] = (int) a/(pFactor*rad);
        } else {
            nPart[i] = 1;
        }
    }
    nPart[3] = np;
    
    if (pp) {
        free(pp);
    }
    pp=malloc(sizeof(particle*)*nPart[0]*nPart[1]*nPart[2]*np);
    for (int i=0; i<nPart[0]*nPart[1]*nPart[2]*np;i++) {
        pp[i]=NULL;
    }
    
    if (pCount) {
        free(pCount);
    }
    pCount=malloc(sizeof(int)*nPart[0]*nPart[1]*nPart[2]);
    for (int i=0; i<nPart[0]*nPart[1]*nPart[2]; i++) {
        pCount[i]=0;
    }

    
    double pPosArray[nPart[0]+nPart[1]+nPart[2]];
    SetupPartitions(rad, a, b, z0, nPart, pPosArray);
    
    dtRelaxation0 = simTimeFinal/relaxationSteps;
    dtRelaxation = dtRelaxation0;
    diffusionSteps = round(simTimeFinal/dtDiffusion0);
    
    if (diffusionSteps<relaxationSteps) {
        printf("Fast relaxation rate: using smaller timestep.\n");
        
        diffusionSteps=relaxationSteps;
        dtDiffusion0=dtRelaxation;
        printf("%e\n",dtDiffusion0);
    }
    
    double avgTotalDiffusionDistance = sqrt(2*diffusionCoeff*simTimeFinal)/(2*rad);
    
    double dtDiffusion = dtDiffusion0;
    double dtTol = 1e-5*dtRelaxation;
    double minSpacing = 1e-7;
    double dtDiffusionMin = timestepScaleFactor*minSpacing*minSpacing/(2*diffusionCoeff);
    
    int jammed = FALSE;
    
    int initialRelaxationStep;
    
    clock_t begin, end;
    double timeSpent;
    
    randinitialize();
    begin = clock();
    
    if(!(readUnfinishedQ||readJammedQ)){
        for (int i=0; i< (int)(np*mixRatio); i++) {
            p[i].rad=rad;
        }
        for (int i=(int)(np*mixRatio); i<np; i++) {
            p[i].rad=rad*ratio;
        }

        FILE *npFile=NULL;
        npFile = fopen("npts.dat","w");
        if (npFile) {
            fprintf(npFile, "%i\n", np);
            fclose(npFile);
        } else {
            printf("npFile pointer is null\n");
        }
        
        FILE *radFile=NULL;
        radFile = fopen("radius.dat","w");
        if (radFile) {
            fprintf(radFile, "%.15lf\n", rad);
            fclose(radFile);
        } else {
            printf("radFile pointer is null\n");
        }
        
        FILE *radiiFile=NULL;
        radiiFile = fopen("radii.dat","w");
        if (radiiFile) {
            for (int i=0; i<np; i++) {
                fprintf(radiiFile, "%.15lf\n", p[i].rad);
            }
            fclose(radiiFile);
        } else {
            printf("radiiFile pointer is null\n");
        }
        FILE *diffTimeFile=NULL;
        diffTimeFile = fopen("diffusionTimescale.dat","w");
        if (diffTimeFile) {
            fprintf(diffTimeFile, "%.15lf\n", diffusionTimescale);
            fclose(diffTimeFile);
        } else {
            printf("diffTimeFile pointer is null\n");
        }
        FILE *abInitFile=NULL;
        abInitFile = fopen("DCParamsInit.dat","w");
        if (abInitFile) {
            fprintf(abInitFile, "%.15lf %.15lf %.15lf\n", z0, a0, b0);
            fclose(abInitFile);
        } else {
            printf("abInitFile pointer is null\n");
        }
        FILE *scalingFactorFile=NULL;
        scalingFactorFile = fopen("scalingFactor.dat","w");
        if (scalingFactorFile) {
            fprintf(scalingFactorFile, "%.15lf\n", scalingFactor);
            fclose(scalingFactorFile);
        } else {
            printf("scalingFactorFile pointer is null\n");
        }

    }
    
    printf("np: %i\n",np);
    printf("rad: %f\n",rad);
    printf("a: %f\n",a);
    printf("b: %f\n",b);
    printf("ratio: %f\n",ratio);
    printf("mixRatio: %f\n",mixRatio);
    
    if (!readUnfinishedQ) {
        printf("maxUndoStepsStart: %i\n",maxUndoStepsStart);
        printf("startupDiffusionSteps: %i\n",startupDiffusionSteps);
        printf("startupDelta: %e\n",startupDelta);
    }
    
    printf("simTimeStart: %e\n",simTimeStart);
    printf("simTimeFinal: %e\n",simTimeFinal);
    
    printf("dtRelaxation: %e\n",dtRelaxation);
    printf("relaxationSteps: %i\n",relaxationSteps);
    printf("maxUndoSteps: %i\n",maxUndoSteps);
    
    printf("diffusionTimescale: %e\n", diffusionTimescale);
    printf("avgTotalDiffusionDistance: %e\n",avgTotalDiffusionDistance);
    printf("dtDiffusion: %e\n",dtDiffusion);
    printf("diffusionSteps: %lf\n",diffusionSteps);
    printf("delta: %e\n",delta);
    printf("diffusionCoeff: %e\n",diffusionCoeff);
    printf("scalingFactor: %lf\n", scalingFactor);
    printf("scalingConst: %e\n",scalingConst);

    if(!(readUnfinishedQ||readJammedQ)){
        FILE* init1 = NULL;
        init1 = fopen("init1.asc","w");
        if(!init1){
            printf("file init1.asc failed to open");
            exit(1);
        }
        for (int i = 0; i < np; i++) {
            // place particles in surface parameter space, based on what would give an even distribution on a sphere
            double u = DCconformalInvert(z0, a, b, genrand_real1());
            double v = 2*PI*genrand_real2();
            
            p[i].position[0]=sqrt((double)1/2*(b*b-2*a*a*u*u+2*a*u*z0+sqrt(b*b*b*b+4*a*b*b*u*(-a*u+z0)+4*a*a*(-a*u+z0)*(-a*u+z0))))*cos(v);
            p[i].position[1]=sqrt((double)1/2*(b*b-2*a*a*u*u+2*a*u*z0+sqrt(b*b*b*b+4*a*b*b*u*(-a*u+z0)+4*a*a*(-a*u+z0)*(-a*u+z0))))*sin(v);
            if (genrand_real2()>0.5)
                p[i].position[2]=a*u;
            else
                p[i].position[2]=-a*u;
            fprintf(init1, "%.15lf %.15lf %.15lf\n", p[i].position[0], p[i].position[1], p[i].position[2]);
        }
        if (init1) fclose(init1);

        DCundoOverlaps(z0, a, b, rad, p, nPart, pCount, pp, pPosArray, maxUndoStepsStart, &dtOverlapStartup, &jammed);
        
        if (jammed) {
            printf("overcrowded starting point - use fewer particles or smaller radius\n");
            exit(1);
        }
        
        FILE *init2=NULL; // after relaxation
        init2 = fopen("init2.asc","w");
        if (!init2) {
            printf("file init2.asc failed to open");
            exit(1);
        }
        for (int i = 0; i < np; i++) {
            fprintf(init2, "%.15lf %.15lf %.15lf\n", p[i].position[0], p[i].position[1], p[i].position[2]);
        }
        fclose(init2);

        for (int i = 0; i < startupDiffusionSteps; i++) {
            
            if(i%1000==0) printf("startup shuffle: %i\n",i);
            
            // one set of np moves
            for (int j0 = 0; j0 < np; j0++) {
                int j = randInteger(np)-1;
                
                vectorCopy(p[j].position, p[j].oldPosition);
                DCstep(z0, a, b, startupDelta, p[j].position, p[j].position);
                
                if ((DCconstraintFunction(z0,a,b,p[j].position))){
                    vectorCopy(p[j].oldPosition, p[j].position);
                    continue;
                }

                //addStochasticForce(a, b, startupDiffusionCoeff, rad, &p[j]);
                //gradDescentConstrainedIntegrate(a,b,&p[j],dtDiffusion);
                rePartitionParticle(&p[j], nPart, pCount, pp, pPosArray );
                //    partitionParticles(p, nPart, pCount, pp, pPosArray);
                
                
                if ( anyOverlapQ(nPart, 2*rad, &p[j], pCount, pp, pPosArray )  ) {
                    //printf("overlap! %i %i\n", j, overlapCheck);
                    vectorCopy(p[j].oldPosition, p[j].position);
                    rePartitionParticle(&p[j], nPart, pCount, pp, pPosArray );
                    //               partitionParticles(p, nPart, pCount, pp, pPosArray);
                }
            }
        }
        
        FILE *init3=NULL; // after initial diffusion
        init3 = fopen("init3.asc","w");
        if (!init3) {
            printf("file init3.asc failed to open");
            exit(1);
        }
        for (int i = 0; i < np; i++) {
            fprintf(init3, "%.15lf %.15lf %.15lf\n", p[i].position[0], p[i].position[1], p[i].position[2]);
            //printf("%3d: %.15lf %.15lf %.15lf\n", i, p[i].position[0], p[i].position[1], p[i].position[2]);
        }
        fclose(init3);
    }
    
    if (animationQ) {
        
        // record a sequence of configurations throughout the simulation
        animationFile = fopen("animation.dat", "a");
        if (animationFile) {
            fprintf(animationFile, "%.15lf\n", simTimeStart/simTimeFinal);
            fprintf(animationFile, "%.15lf %.15lf %.15lf\n", z0, a, b);
            for (int i=0; i<np; i++) {
                fprintf(animationFile, "%.15lf %.15lf %.15lf\n", p[i].position[0], p[i].position[1], p[i].position[2]);
            }
            fclose(animationFile);
        } else {
            printf("animationFile pointer is null\n");
        }
    }
    
    int overlapCount = 0;
    int stepCount = 0;
    
    double simTimeLastRelax = simTimeStart;
    double nextRelaxationTimeLastRelax = nextRelaxationTime;
    int relaxationStep = 0;
    
    partitionParticles(p, nPart, pCount, pp, pPosArray);
    
    double typicalSpacing = findTypicalSpacing(np, rad, p, nPart, pCount, pp, pPosArray);
    
    for(simTime=simTimeStart;simTime<=simTimeFinal && !jammed;simTime+=dtDiffusion){
        if(simTime > nextRelaxationTime){
            //printf("simTime/simTimeFinal %.16f\n", simTime/simTimeFinal);
            aOld = a;
            bOld = b;
            for (int j=0; j<np; j++)
                vectorCopy(p[j].position, p[j].preRelaxPosition); // not being used... is this still useful?
            if (dynamicsQ) {
                b = scalingConst*pow(simTime, scalingFactor);
                a = findA(b, volume, DBL_EPSILON, MACHINE_EPSILON);
            } else {
                a=a0*(1-simTime/simTimeFinal);
                b=b0*(1-simTime/simTimeFinal);
            }
            for(int j=0;j<np;j++){
                vectorCopy(p[j].position, p[j].oldPosition);
                if (moveWithSurfaceQ)
                    moveWithSurface(aOld, bOld, a, b, p[j].position, p[j].position);
                else
                    DCRelaxationProject(z0, aOld, bOld, a, b, &p[j]);
            }
            
           /* if(  (int) ( (simTime/simTimeFinal) * 100000 ) % 10000 == 0 ){ //this should be changed to be a function of pFactor and aa
                
                for (int i=0; i<3; i++) {
                    if( (int) a/(pFactor*rad) > 1 ){
                        nPart[i] = (int) a/(pFactor*rad);
                    } else {
                        nPart[i] = 1;
                    }
                }
                nPart[3] = np;
                
                SetupPartitions(rad, a, b, z0, nPart, pPosArray);
                
                partitionParticles(p, nPart, pCount, pp, pPosArray);
            }*/
            
            int nearJammed = FALSE;
            DCundoOverlaps(z0, a, b, rad, p, nPart, pCount, pp, pPosArray, maxUndoSteps, &dtOverlap, &nearJammed);
            
            if (nearJammed) {
                for (int j=0; j<np; j++) {
                    vectorCopy(p[j].postPreviousRelaxPosition, p[j].position);
                }
                a=aOld;
                b=bOld;
                
                simTime = simTimeLastRelax;
                printf("simTime %e\n",simTime);
                nextRelaxationTime = nextRelaxationTimeLastRelax;
                dtRelaxation *= 0.5;
                dtOverlap *= 0.5;
                
                printf("reducing timesteps %e %e\n",dtRelaxation,dtOverlap);
                
                if (dtRelaxation < dtTol) jammed = TRUE;
            } else {
                simTimeLastRelax = simTime;
                nextRelaxationTimeLastRelax = nextRelaxationTime; // save old time, which will be updated using new timestep on reset
                for (int j=0; j<np; j++)
                    vectorCopy(p[j].position, p[j].postPreviousRelaxPosition);
                // move aOld here for consistency?
                
                typicalSpacing = findTypicalSpacing(np, rad, p, nPart, pCount, pp, pPosArray);
                {
                    dtDiffusion = timestepScaleFactor*typicalSpacing*typicalSpacing/(2*diffusionCoeff);
                    if (dtDiffusion>dtDiffusion0) dtDiffusion = dtDiffusion0;
                    if (dtDiffusion<dtDiffusionMin) dtDiffusion = dtDiffusionMin;
                }
                
                if(  (int) ( (simTime/simTimeFinal) * 100000 ) % 10000 == 0 ){
                    for (int i=0; i<3; i++) {
                        if( (int) a/(pFactor*rad) > 1 ){
                            nPart[i] = (int) a/(pFactor*rad);
                        } else {
                            nPart[i] = 1;
                        }
                    }
                    nPart[3] = np;
                
                    SetupPartitions(rad, a, b, z0, nPart, pPosArray);
                }
                
                
            }
            
            if (animationQ && relaxationStep % animationRate == 0) {
                animationFile = fopen("animation.dat", "a");
                if (animationFile) {
                    fprintf(animationFile, "%.15lf\n", simTime/simTimeFinal);
                    fprintf(animationFile, "%.15lf %.15lf %.15lf\n", z0, a, b);
                    for (int j=0; j<np; j++) {
                        fprintf(animationFile, "%.15lf %.15lf %.15lf\n", p[j].position[0], p[j].position[1], p[j].position[2]);
                    }
                    fclose(animationFile);
                } else {
                    printf("animationFile pointer is null\n");
                }
            }
            
            partitionParticles(p, nPart, pCount, pp, pPosArray);
            nextRelaxationTime += dtRelaxation;
            relaxationStep++;

        }
        
        //One set of particle moves
        
        // reset forces
        for (int k=0; k<np; k++)
            for (int l=0; l<3; l++) {
                p[k].forceDet[l]=0;
                p[k].forceStoc[l]=0;
            }
            
        // calculate forces
        double netForcePreSum[np][3];
        for (int k = 0; k < np; k++)
            for (int m=0; m<3; m++)
                netForcePreSum[k][m]=0;
            
        for (int k = 0; k < np; k++) {
            double newForce[3];
            for (int m=0; m<3; m++)
                newForce[m]=0;
/*
            if (gravityQ && !relaxationOnlyQ) {
                addGravitationalForce(&p[k],gravStrength);
            }
            if (unjamBiasQ && !relaxationOnlyQ) {
                addGenericForce(&p[k], biasForce[k], biasStrength);
            }
            if (harmonicSurfaceQ) {
                addHarmonicSurfaceForce(a, b, &p[k], surfaceStrength);
                energy+=harmonicSurfaceEnergy(a, b, &p[k], surfaceStrength);
            }
 */
            addStochasticForce(a, b, diffusionCoeff, rad, &p[k]);
        
                // repulsion
                /*for (int m=0; m < k; m++) {
                 addRepulsiveForce(&p[k], &p[m], attrStrength, newForce);
                 configEnergy += repulsiveEnergy(&p[k], &p[m], attrStrength);
                 }*/
                
        }
            
        // going to loop through particles in a random order - set up shuffled array of indices
        int randIndices[np];
        for (int k=0; k<np; k++)
            randIndices[k]=k;
        for (int k=np-1; k>0; k--) {
            int swapWith = randInteger(k-1);
            int hold = randIndices[k];
            randIndices[k] = randIndices[swapWith];
            randIndices[swapWith] = hold;
        }
            
            
        // integrate forces
        for (int k0 = 0; k0 < np; k0++) {
            int k = randIndices[k0];
                
            // apply forces to make move
            vectorCopy(p[k].position, p[k].oldPosition);
                
            // move and project
            //gradDescentIntegrate(&p[k],dtDiffusion);
            //project(a, b, p[k].position, p[k].position);
                
            DCgradDescentConstrainedIntegrate(z0, a, b, &p[k], dtDiffusion);
            
            rePartitionParticle(&p[k], nPart, pCount, pp, pPosArray );
            
            stepCount++;
            if ( anyOverlapQ(nPart, 2*rad, &p[k], pCount, pp, pPosArray )||isnan(DCconstraintFunction(z0,a,b,p[k].position)) ) {
                //printf("overlap! %i %i\n", k, overlapCheck);
                vectorCopy(p[k].oldPosition, p[k].position);
                rePartitionParticle(&p[k], nPart, pCount, pp, pPosArray );
                overlapCount++;
            }
                
        }
        
        if (simStep%outputSteps==0) {
            FILE *configurationTemp=NULL;
            configurationTemp = fopen("configurationTemp.asc.tmp","w");
            if (configurationTemp) {
                for (int i=0; i<np; i++) {
                    fprintf(configurationTemp, "%.15lf %.15lf %.15lf\n", p[i].position[0], p[i].position[1], p[i].position[2]);
                }
                fclose(configurationTemp);
            } else {
                printf("configurationTemp pointer is null\n");
            }
            
            FILE *paramsTempFile=NULL;
            paramsTempFile = fopen("ParamsTemp.dat.tmp","w");
            if (paramsTempFile) {
                fprintf(paramsTempFile, "%.15lf %.15lf %.15lf\n", z0, a, b);
                fclose(paramsTempFile);
            } else {
                printf("paramsTempFile pointer is null\n");
            }
            
            FILE *simTimeTempFile=NULL;
            simTimeTempFile = fopen("simTimeTemp.dat.tmp","w");
            if (simTimeTempFile) {
                fprintf(simTimeTempFile, "%lf\n", simTime);
                fclose(simTimeTempFile);
            } else {
                printf("simTimeTempFile pointer is null\n");
            }
            
            
            FILE *nextRelaxationFile=NULL;
            nextRelaxationFile = fopen("nextRelaxationTime.dat.tmp","w");
            if (nextRelaxationFile) {
                fprintf(nextRelaxationFile, "%lf\n", nextRelaxationTime);
                fclose(nextRelaxationFile);
            } else {
                printf("nextRelaxationFile pointer is null\n");
                exit(1);
            }
            
            system("mv -f configurationTemp.asc.tmp configurationTemp.asc");
            system("mv -f ParamsTemp.dat.tmp ParamsTemp.dat");
            system("mv -f simTimeTemp.dat.tmp simTimeTemp.dat");
            system("mv -f nextRelaxationTime.dat.tmp nextRelaxationTime.dat");
        }
        

        //printf("simTime/simTimeFinal: %.16lf\n", simTime/simTimeFinal);
        simStep++;

    }
    
    end = clock();
    timeSpent = (double)(end - begin) / CLOCKS_PER_SEC;
    
    FILE *timeFile=NULL;
    timeFile = fopen("time.dat","w");
    if (timeFile) {
        fprintf(timeFile, "%lf\n", timeSpent);
        fclose(timeFile);
    } else {
        printf("timeFile pointer is null\n");
    }
    
    FILE *paramFile=NULL;
    paramFile = fopen("DCParams.dat","w");
    if (paramFile) {
        fprintf(paramFile, "%.15lf %.15lf %.15lf\n", z0, a, b);
        fclose(paramFile);
    } else {
        printf("paramFile pointer is null\n");
    }
    
    FILE *configuration=NULL;
    configuration = fopen("configuration.asc","w");
    if (configuration) {
        for (int i=0; i<np; i++) {
            fprintf(configuration, "%.15lf %.15lf %.15lf\n", p[i].position[0], p[i].position[1], p[i].position[2]);
        }
        fclose(configuration);
    } else {
        printf("configuration pointer is null\n");
    }
    
    
    FILE *simTimeArrestedFile=NULL;
    simTimeArrestedFile = fopen("simTimeArrested.dat","w");
    if (simTimeArrestedFile) {
        fprintf(simTimeArrestedFile, "%lf\n", simTime);
        fclose(simTimeArrestedFile);
    } else {
        printf("simTimeArrestedFile pointer is null\n");
    }
    
    
    FILE *nextRelaxationArrestedFile=NULL;
    nextRelaxationArrestedFile = fopen("nextRelaxationTimeArrested.dat","w");
    if (nextRelaxationArrestedFile) {
        fprintf(nextRelaxationArrestedFile, "%lf\n", nextRelaxationTime);
        fclose(nextRelaxationArrestedFile);
    } else {
        printf("nextRelaxationArrestedFile pointer is null\n");
        exit(1);
    }

    printf("Hello, World!\n");
    return 0;
}
