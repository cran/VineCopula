#ifndef VINE_H
#define VINE_H

/*
** vine.c - C code of the package VineCopula
**
** with contributions from Carlos Almeida, Aleksey Min,
** Ulf Schepsmeier, Jakob Stoeber and Eike Brechmann
**
** A first version was based on code
** from Daniel Berg <daniel at danielberg.no>
** provided by personal communication.
**
*/

#ifndef   	VINE_H_
# define   	VINE_H_
#endif 	    /* !VINE_H_ */



#include <R.h>
#include <Rinternals.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <setjmp.h>
#include <Rdefines.h>
#include <Rmath.h>
#include "VineCopula/vine.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#define XINFMAX DBL_MAX

#define STRICT_R_HEADERS 1

/* define boolean type for C */
typedef unsigned int boolean;

#define ep 1e-6

#define UMAX  1-1e-12
#define UMIN  1e-12

#define TOL 1e-4


#endif
