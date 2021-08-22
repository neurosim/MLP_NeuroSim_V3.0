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
*   Pai-Yu Chen	    Email: pchen72 at asu dot edu 
*                    
*   Xiaochen Peng   Email: xpeng15 at asu dot edu
********************************************************************************/

#include <cmath>
#include <iostream>
#include "constant.h"
#include "formula.h"
#include "MultilevelSenseAmp.h"

using namespace std;

MultilevelSenseAmp::MultilevelSenseAmp(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), currentSenseAmp(_inputParameter, _tech, _cell), FunctionUnit() {
	initialized = false;
}

void MultilevelSenseAmp::Initialize(int _numCol, int _levelOutput, double _clkFreq, int _numReadCellPerOperationNeuro) {
	if (initialized) {
		cout << "[MultilevelSenseAmp] Warning: Already initialized!" << endl;
    } else {

	numCol = _numCol;
	levelOutput = _levelOutput;                // # of bits for A/D output ... 
	clkFreq = _clkFreq;
	numReadCellPerOperationNeuro = _numReadCellPerOperationNeuro;
	
	// Initialize SenseAmp
	currentSenseAmp.Initialize(levelOutput*numCol, false, false, clkFreq, numReadCellPerOperationNeuro);        // use real-traced mode ... 

	initialized = true;
	}
}

void MultilevelSenseAmp::CalculateArea(double heightArray, double widthArray, AreaModify _option) {
	if (!initialized) {
		cout << "[MultilevelSenseAmp] Error: Require initialization first!" << endl;
	} else {
		
		if (widthArray && _option==NONE) {
			currentSenseAmp.CalculateUnitArea();
			currentSenseAmp.CalculateArea(widthArray);
		} else {
			cout << "[MultilevelSenseAmp] Error: No width assigned for the multiSenseAmp circuit" << endl;
			exit(-1);
		}
		
		
		
		// Assume the Current Mirrors are on the same row and the total width of them is smaller than the adder or DFF
		area = currentSenseAmp.area;
		width = widthArray;
		height = area / width;
		
		// Modify layout
		newHeight = heightArray;
		newWidth = widthArray;
		switch (_option) {
			case MAGIC:
				MagicLayout();
				break;
			case OVERRIDE:
				OverrideLayout();
				break;  
			default:    // NONE
				break;
		}

	}
}

void MultilevelSenseAmp::CalculateLatency(double numRead) {
	if (!initialized) {
		cout << "[MultilevelSenseAmp] Error: Require initialization first!" << endl;
	} else {
		readLatency = 0;
		
		readLatency += 1e-9;
		
		readLatency *= numRead;
	}
}

void MultilevelSenseAmp::CalculatePower(double numof1, double numof2, double numof3, double numof4, double numof5, double numof6, double numof7, double numof8, double numof9, double numof10, double numRead) {
	if (!initialized) {
		cout << "[MultilevelSenseAmp] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;
		
		if (numof1+numof2+numof3+numof4+numof5+numof6+numof7+numof8+numof9+numof10 != 0) {     // trace mode
			if (levelOutput == 2) {
			    readDynamicEnergy += numof1*7.08e-15 + numof2*6.87e-15 + numof3*6.83-15 + numof4*6.82e-15 + numof5*6.81e-15 + numof6*6.81e-15 + numof7*6.80e-15 + numof8*6.80e-15 + numof9*6.80e-15 + numof10*6.80e-15;
		    } else if (levelOutput == 4) {
			    readDynamicEnergy += numof1*2.12e-14 + numof2*2.06e-14 + numof3*2.05e-14 + numof4*2.05e-14 + numof5*2.04e-14 + numof6*2.04e-14 + numof7*2.04e-14 + numof8*2.04e-14 + numof9*2.04e-14 + numof10*2.04e-14;
		    } else if (levelOutput == 8) {
			    readDynamicEnergy += numof1*4.96e-14 + numof2*4.82e-14 + numof3*4.79e-14 + numof4*4.78e-14 + numof5*4.78e-14 + numof6*4.77e-14 + numof7*4.77e-14 + numof8*4.77e-14 + numof9*4.77e-14 + numof10*4.77e-14;
		    } else if (levelOutput == 16) {
			    readDynamicEnergy += numof1*1.07e-13 + numof2*1.03e-13 + numof3*1.03e-13 + numof4*1.03e-13 + numof5*1.02e-13 + numof6*1.02e-13 + numof7*1.02e-13 + numof8*1.02e-13 + numof9*1.02e-13 + numof10*1.02e-13;
		    } else if (levelOutput == 32) {
			    readDynamicEnergy += numof1*2.20e-13 + numof2*2.14e-13 + numof3*2.13e-13 + numof4*2.12e-13 + numof5*2.12e-13 + numof6*2.12e-13 + numof7*2.12e-13 + numof8*2.12e-13 + numof9*2.12e-13 + numof10*2.12e-13;
		    } else if (levelOutput == 64) {
			    readDynamicEnergy += numof1*4.50e-13 + numof2*4.37e-13 + numof3*4.34e-13 + numof4*4.33e-13 + numof5*4.33e-13 + numof6*4.33e-13 + numof7*4.32e-13 + numof8*4.32e-13 + numof9*4.32e-13 + numof10*4.32e-13;
		    }
		} else {    // average case
		    
			// 64*64_8levels: 6.019e-13,     64*64_4levels: 2.255e-13,    128*128_8levels: 6.902e-13
			readDynamicEnergy = 6.902e-13; 
            readDynamicEnergy *= numCol;			
		    readDynamicEnergy *= numRead;
		}

		if (!readLatency) {
			//cout << "[MultilevelSenseAmp] Error: Need to calculate read latency first" << endl;
		} else {
			readPower = readDynamicEnergy/readLatency;
		}
	}
}

void MultilevelSenseAmp::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}


