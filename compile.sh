#!/bin/bash
gcc -lm -std=c99 -O3 -g main.c partition.c particles.c surface.c forces.c misc.c mt19937ar.c DropletCoalescence.c
