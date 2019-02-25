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
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include "formula.h"
#include "Param.h"
#include "Cell.h"
#include "Array.h"

extern Param *param;
extern Array *arrayIH;
extern Array *arrayHO;
extern std::vector< std::vector<double> > Input;
extern std::vector< std::vector<int> > dInput;
extern std::vector< std::vector<double> > testInput;
extern std::vector< std::vector<int> > dTestInput;
extern std::vector< std::vector<double> > Output;
extern std::vector< std::vector<double> > testOutput;

extern std::vector< std::vector<double> > weight1;
extern std::vector< std::vector<double> > weight2;
extern std::vector< std::vector<double> > deltaWeight1;
extern std::vector< std::vector<double> > deltaWeight2;
extern std::vector<std::vector<double> >  totalDeltaWeight1;
extern std::vector<std::vector<double> >  totalDeltaWeight1_abs;
extern std::vector<std::vector<double> >  totalDeltaWeight2;
extern std::vector<std::vector<double> >  totalDeltaWeight2_abs;

/* Read trainging data from file */
void ReadTrainingDataFromFile(const char *trainPatchFileName, const char *trainLabelFileName) {
	FILE *fp_patch = fopen(trainPatchFileName, "r");
	FILE *fp_label = fopen(trainLabelFileName, "r");

	if (!fp_patch) {
		std::cout << trainPatchFileName << " cannot be found!\n";
		exit(-1);
	}
	if (!fp_label) {
		std::cout << trainLabelFileName << " cannot be found!\n";
		exit(-1);
	}

	int i = 0;
	int j = 0;
	while (fscanf(fp_patch, "%lf", &Input[i][j]) != EOF){
		Input[i][j] = truncate(Input[i][j], param->numInputLevel - 1, param->BWthreshold);
		dInput[i][j] = round(Input[i][j] * (param->numInputLevel - 1));
		i += 1;
		if (i%param->numMnistTrainImages == 0){
			j += 1;
			i = 0;
		}
	}

	i = 0;
	j = 0;
	int k = 0;
	while (fscanf(fp_label, "%d", &k) != EOF){
		Output[i][k] = 1;
		i += 1;
	}
	fclose(fp_patch);
	fclose(fp_label);
}

/* Read testing data from file */
void ReadTestingDataFromFile(const char *testPatchFileName, const char *testLabelFileName) {
	FILE *fp_patch = fopen(testPatchFileName, "r");
	FILE *fp_label = fopen(testLabelFileName, "r");

	if (!fp_patch) {
		std::cout << testPatchFileName << " cannot be found!\n";
		exit(-1);
	}
	if (!fp_label) {
		std::cout << testLabelFileName << " cannot be found!\n";
		exit(-1);
	}

	int i = 0;
	int j = 0;
	while (fscanf(fp_patch, "%lf", &testInput[i][j]) != EOF){
		testInput[i][j] = truncate(testInput[i][j], param->numInputLevel - 1, param->BWthreshold);
		dTestInput[i][j] = round(testInput[i][j] * (param->numInputLevel - 1));
		i += 1;
		if (i%param->numMnistTestImages == 0){
			j += 1;
			i = 0;
		}
	}
	i = 0;
	j = 0;
	int k = 0;
	while (fscanf(fp_label, "%d", &k) != EOF){
		testOutput[i][k] = 1;
		i += 1;
	}

	fclose(fp_patch);
	fclose(fp_label);
}

/* Print weight to file */
void PrintWeightToFile(const char *str) {
	/* Print weight1 */
	char printWeight1FileName[50];
	sprintf(printWeight1FileName, "%s1.csv", str);
	FILE *fp_dw1 = fopen(printWeight1FileName, "w");
	fprintf(fp_dw1, "minWeight=%f, maxWeight=%f\n", param->minWeight, param->maxWeight);
	for (int j = 0; j < param->nHide; j++){
		for (int k = 0; k < param->nInput; k++){
			    fprintf(fp_dw1, "%f,", weight1[j][k]);
    }
		fprintf(fp_dw1, "\n");
	}
	fclose(fp_dw1);
	/* Print weight2 */
	char printWeight2FileName[50];
	sprintf(printWeight2FileName, "%s2.csv", str);
	FILE *fp_dw2 = fopen(printWeight2FileName, "w");
	fprintf(fp_dw2, "minWeight=%f, maxWeight=%f\n", param->minWeight, param->maxWeight);
	for (int j = 0; j < param->nOutput; j++){
		for (int k = 0; k < param->nHide; k++){
			fprintf(fp_dw2, "%f,", weight2[j][k]);
    }
		fprintf(fp_dw2, "\n");
	}
	fclose(fp_dw2);
}

