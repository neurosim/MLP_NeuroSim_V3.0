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

#include <string>

#ifndef PARAM_H_
#define PARAM_H_

class Param {
public:
	Param();

	/* MNIST dataset */
	int numMnistTrainImages;// # of training images in MNIST
	int numMnistTestImages;	// # of testing images in MNIST
	
	/* Algorithm parameters */
	int numTrainImagesPerEpoch;	// # of training images per epoch
	int totalNumEpochs;	// Total number of epochs
	int interNumEpochs;	// Internal number of epochs (print out the results every interNumEpochs)
	int nInput;     // # of neurons in input layer
	int nHide;      // # of neurons in hidden layer
	int nOutput;	// # of neurons in output layer
	double alpha1;		// Learning rate for the synapses from input to hidden layer
	double alpha2;		// Learning rate for the synapses from hidden to output layer
	double maxWeight;	// Upper bound of weight value
	double minWeight;	// Lower bound of weight value
    char* optimization_type;

	/* Hardware parameters */
	bool useHardwareInTrainingFF;   // Use hardware in the feed forward part of training or not (true: realistic hardware, false: ideal software)
	bool useHardwareInTrainingWU;   // Use hardware in the weight update part of training or not (true: realistic hardware, false: ideal software)
	bool useHardwareInTraining;		// Use hardware in the training or not
	bool useHardwareInTestingFF;    // Use hardware in the feed forward part of testing or not (true: realistic hardware, false: ideal software)
	int numBitInput;		// # of bits of the input data (=1 for black and white data)
	int numBitPartialSum;	// # of bits of the digital output (partial weighted sum output)
	int pSumMaxHardware;	// Max digital output value of partial weighted sum
	int numInputLevel;	// # of levels of the input data
	int numWeightBit;	// # of weight bits (only for pure algorithm, SRAM and digital RRAM hardware)
	double BWthreshold; // The black and white threshold for numBitInput=1
	double Hthreshold;	// The spiking threshold for the hidden layer (da1 in Train.cpp and Test.cpp)
	int numColMuxed;	// How many columns share 1 read circuit (for analog RRAM) or 1 S/A (for digital RRAM)
	int numWriteColMuxed;	// How many columns share 1 write column decoder driver (for digital RRAM)
	bool writeEnergyReport;	// Report write energy calculation or not
	bool NeuroSimDynamicPerformance; // Report the dynamic performance (latency and energy) in NeuroSim or not
	bool relaxArrayCellHeight;	// True: relax the array cell height to standard logic cell height in the synaptic array
	bool relaxArrayCellWidth;	// True: relax the array cell width to standard logic cell width in the synaptic array
	double arrayWireWidth;	// Array wire width (nm)
	int processNode;	// Technology node (nm)
	double clkFreq;		// Clock frequency (Hz)
};

#endif
