/*******************************************************************************
* Copyright (c) 2015-2017
* School of Electrical, Computer and Energy Engineering, Arizona State University
* PI: Prof. Shimeng Yu
* All rights reserved.
*   
* This source code is part of NeuroSim - a device-circuit-algorithm framework to benchmark 
* neuro-inspired architectures with synaptic devices(e.g., SRAM and emerging non-volatile memory). 
* Copyright of the model is maintained by the developers, and the model is distributed under 
* the terms of the Creative Commons Attribution-NonCommercial 4.0 International Public License 
* http://creativecommons.org/licenses/by-nc/4.0/legalcode.
* The source code is free and you can redistribute and/or modify it
* by providing that the following conditions are met:
*   
*  1) Redistributions of source code must retain the above copyright notice,
*     this list of conditions and the following disclaimer. 
*   
*  2) Redistributions in binary form must reproduce the above copyright notice,
*     this list of conditions and the following disclaimer in the documentation
*     and/or other materials provided with the distribution.
*   
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
* 
* Developer list: 
*   Pai-Yu Chen     Email: pchen72 at asu dot edu 
*                     
*   Xiaochen Peng   Email: xpeng15 at asu dot edu
********************************************************************************/

#include <cstdio>
#include <vector>
#include <random>
#include "Param.h"
#include "Array.h"
#include "NeuroSim.h"

extern Param *param;

extern std::vector< std::vector<double> > weight1;
extern std::vector< std::vector<double> > weight2;

extern Array *arrayIH;
extern Array *arrayHO;

/* Weights initialization */
void WeightInitialize() {
    srand(2);
    /* Initialize weights for the input layer */
    for (int i = 0; i < param->nHide; i++) {
        for (int j = 0; j < param->nInput; j++) {
            weight1[i][j] = (double)(rand() % 7 +(-3) ) / 3;   // random number: 0, 0.33, 0.66 or 1
            //printf("weight 1 is %.4f\n", weight1[i][j]);
        }
    }
    /* Initialize weights for the hidden layer */
    for (int i = 0; i < param->nOutput; i++) {
        for (int j = 0; j < param->nHide; j++) {
            weight2[i][j] = (double)(rand() % 7 +(-3) ) / 3;   // random number: 0, 0.33, 0.66 or 1
        }
    }
}

/* Conductance initialization (map weight to RRAM conductance or SRAM data) */
void WeightToConductance() {
    /* Erase the weight of arrayIH */
    for (int col=0; col<param->nHide; col++) {
        for (int row=0; row<param->nInput; row++) {
            arrayIH->WriteCell(col, row, -(param->maxWeight-param->minWeight), 0 /* delta_W=-(param->maxWeight-param->minWeight) will completely erase */, param->maxWeight, param->minWeight, false);
        }
    }
    /* Erase the weight of arrayHO */
    for (int col=0; col<param->nOutput; col++) {
        for (int row=0; row<param->nHide; row++) {
            arrayHO->WriteCell(col, row, -(param->maxWeight-param->minWeight), 0 /* delta_W=-(param->maxWeight-param->minWeight) will completely erase */, param->maxWeight, param->minWeight, false);
        }
    }
    /* Write weight to arrayIH */
    for (int col=0; col<param->nHide; col++) {
        for (int row=0; row<param->nInput; row++) {
            arrayIH->WriteCell(col, row, weight1[col][row], weight1[col][row], param->maxWeight, param->minWeight, false);
        }
    }
    /* Write weight to arrayHO */
    for (int col=0; col<param->nOutput; col++) {
        for (int row=0; row<param->nHide; row++) {
            arrayHO->WriteCell(col, row, weight2[col][row], weight2[col][row], param->maxWeight, param->minWeight, false);
        }
    }
}

/* Mapping from analog current to digital output*/
int CurrentToDigits(double I /* current */, double Imax /* max current */) {
    return (int)(I / (Imax/param->pSumMaxHardware));
}

/* Mapping from hardware digital output to algorithm value*/
double DigitsToAlgorithm(int outputDigits /* output digits from ADC */, double pSumMaxAlgorithm /* max value of partial weighted sum in algorithm */) {
    return ((double)outputDigits / param->pSumMaxHardware) * pSumMaxAlgorithm;
}

