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

// This file cannot be compiled alone. Only include this file in main.cpp.

/* Global variables */
Param *param = new Param(); // Parameter set

/* Inputs of training set */
std::vector< std::vector<double> >
Input(param->numMnistTrainImages, std::vector<double>(param->nInput));
/* Outputs of training set */
std::vector< std::vector<double> >
Output(param->numMnistTrainImages, std::vector<double>(param->nOutput));

/* Weights from input to hidden layer */
std::vector< std::vector<double> >
weight1(param->nHide, std::vector<double>(param->nInput));
/* Weights from hidden layer to output layer */
std::vector< std::vector<double> >
weight2(param->nOutput, std::vector<double>(param->nHide));

/* Weight change of weight1 */
std::vector< std::vector<double> >
deltaWeight1(param->nHide, std::vector<double>(param->nInput));

/* Weight change of weight2 */
std::vector< std::vector<double> >
deltaWeight2(param->nOutput, std::vector<double>(param->nHide));

/*the variables to track the ΔW*/
std::vector< std::vector<double> >
totalDeltaWeight1(param->nHide, std::vector<double>(param->nInput));
std::vector< std::vector<double> >
totalDeltaWeight1_abs(param->nHide, std::vector<double>(param->nInput));
/*the variables to track the ΔW*/
std::vector< std::vector<double> >
totalDeltaWeight2(param->nOutput, std::vector<double>(param->nHide));
std::vector< std::vector<double> >
totalDeltaWeight2_abs(param->nOutput, std::vector<double>(param->nHide));

/* Inputs of testing set */
std::vector< std::vector<double> >
testInput(param->numMnistTestImages, std::vector<double>(param->nInput));
/* Outputs of testing set */
std::vector< std::vector<double> >
testOutput(param->numMnistTestImages, std::vector<double>(param->nOutput));

/* Digitized inputs of training set (an integer between 0 to 2^numBitInput-1) */
std::vector< std::vector<int> >
dInput(param->numMnistTrainImages, std::vector<int>(param->nInput));
/* Digitized inputs of testing set (an integer between 0 to 2^numBitInput-1) */
std::vector< std::vector<int> >
dTestInput(param->numMnistTestImages, std::vector<int>(param->nInput));

// the arrays for optimization
std::vector< std::vector<double> > 
gradSquarePrev1(param->nHide, std::vector<double>(param->nInput));
std::vector< std::vector<double> >
gradSquarePrev2(param->nOutput, std::vector<double>(param->nHide));
std::vector< std::vector<double> > 
gradSum1(param->nHide, std::vector<double>(param->nInput));
std::vector< std::vector<double> >
gradSum2(param->nOutput, std::vector<double>(param->nHide));
std::vector< std::vector<double> >
momentumPrev1(param->nHide, std::vector<double>(param->nInput));
std::vector< std::vector<double> >
momentumPrev2(param->nOutput, std::vector<double>(param->nHide));


/* # of correct prediction */
int correct = 0;

/* Synaptic array between input and hidden layer */
Array *arrayIH = new Array(param->nHide, param->nInput, param->arrayWireWidth);
/* Synaptic array between hidden and output layer */
Array *arrayHO = new Array(param->nOutput, param->nHide, param->arrayWireWidth);

/* Random number generator engine */
std::mt19937 gen;

/* NeuroSim */
SubArray *subArrayIH;   // NeuroSim synaptic core for arrayIH
SubArray *subArrayHO;   // NeuroSim synaptic core for arrayHO
/* Global properties of subArrayIH */
InputParameter inputParameterIH;
Technology techIH;
MemCell cellIH;
/* Global properties of subArrayHO */
InputParameter inputParameterHO;
Technology techHO;
MemCell cellHO;
/* Neuron peripheries below subArrayIH */
Adder adderIH(inputParameterIH, techIH, cellIH);
Mux muxIH(inputParameterIH, techIH, cellIH);
RowDecoder muxDecoderIH(inputParameterIH, techIH, cellIH);
DFF dffIH(inputParameterIH, techIH, cellIH);
Subtractor subtractorIH(inputParameterIH, techIH, cellIH);
/* Neuron peripheries below subArrayHO */
Adder adderHO(inputParameterHO, techHO, cellHO);
Mux muxHO(inputParameterHO, techHO, cellHO);
RowDecoder muxDecoderHO(inputParameterHO, techHO, cellHO);
DFF dffHO(inputParameterHO, techHO, cellHO);
Subtractor subtractorHO(inputParameterHO, techHO, cellHO);
