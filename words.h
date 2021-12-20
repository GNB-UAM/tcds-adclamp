////////////////////////////////////////////////////////////////////////////////
/*!
   \file words.h
   \date 29-01-2012
   \author Grupo Neurocomputacion Biologica (GNB). Spain

   \brief File with words module definition

   File with data structures, constant and functions for hr3 model

   \version 1.0
*/
////////////////////////////////////////////////////////////////////////////////


#ifndef __HR3_H__
#define __HR3_H__

#include "model.h"
#include "common.h"
#include "int.h"
#include "wordsBuffer.h"
#include <linux/types.h>
#include <linux/random.h>

#define BIT_DETECTED_OUT 1
#define BIT_NOT_DETECTED_OUT 0
#define END_WINDOW_OUT 2
#define WORD_DETECTED_OUT 3
#define END_HISTOGRAM_TIME -2

#define FALSE 0
#define TRUE 1

double model(double currentV, double time);
double histogramModel(double currentV, double time);
double detectorModel(double currentV, double time);
double aleatModel(double currentV, double time);
int detect_spike (double voltage, double time);
void nm_func2(void);
void init_neurons (NEURON_MODEL *neuronModel);
void asign_neurons_names (NEURON_MODEL *neuronModel);
NEURON_MODEL *init_nm ( unsigned long n_neurons );
void free_nm (void);
void int_function(double time, double *vars, double * dvars, double *p);

static inline char *ftoa(double x, char *str) {
    char s = (x >= 0 ? ' ' : '-');
    int tmp1, tmp2;

    x = fabs(x);
    tmp1 = floor(x);
    tmp2 = 1e6*(x - tmp1);
    sprintf(str, "%c%d.%06d", s, tmp1, tmp2);
    return str;
}

static inline double vAbs (double x){
    return fabs(x);
}
//!< Number of state variables
#define DIM 3

//!< State variables
#define X 0
#define Y 1
#define Z 2

//!< Number of model parameters
#define DIM_PARAM 3

//!< Parameters
#define I_EXT 0
#define E     1
#define K     2

//!< Initial Values
#define X_0	-1.638766
#define Y_0	-12.181940
#define Z_0	  2.791712

#define I_0   0.0
#define E_0   3.3
#define K_0   1.0

//!< HR Constants, that could be future parameters in the future
#define a_cte         1.0
#define b_cte         3.0
#define c_cte         1.0
#define d_cte         1.0
#define e_cte         1.0
#define f_cte         5.0
#define g_cte         0.0021
#define h1_cte        1.6
#define Svalue_cte    4.0

#define FIRST_CALL 0
#define NEW_TIME_WINDOW 1

#endif
