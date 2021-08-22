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
#include "CurrentSenseAmp.h"

using namespace std;

CurrentSenseAmp::CurrentSenseAmp(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit() {
	// TODO Auto-generated constructor stub
	initialized = false;
	invalid = false;
}

CurrentSenseAmp::~CurrentSenseAmp() {
	// TODO Auto-generated destructor stub
}

void CurrentSenseAmp::Initialize(int _numCol, bool _parallel, bool _rowbyrow, double _clkFreq, int _numReadCellPerOperationNeuro) {
	if (initialized)
		cout << "[Current Sense Amp] Warning: Already initialized!" << endl;

	numCol = _numCol;
	parallel = _parallel;
	rowbyrow = _rowbyrow;
	clkFreq = _clkFreq;
	numReadCellPerOperationNeuro = _numReadCellPerOperationNeuro;
	
	widthNmos = MIN_NMOS_SIZE * tech.featureSize;
	widthPmos = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;

	initialized = true;
}

void CurrentSenseAmp::CalculateUnitArea() {
	if (!initialized) {
		cout << "[CurrentSenseAmp] Error: Require initialization first!" << endl;
	} else {
		double hNmosL, wNmosL, hNmosS, wNmosS, hNmosM, wNmosM;
		
		CalculateGateArea(INV, 1, widthNmos, 0, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hNmosS, &wNmosS);
		CalculateGateArea(INV, 1, widthNmos*2, 0, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hNmosM, &wNmosM);
		CalculateGateArea(INV, 1, widthNmos*6, 0, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hNmosL, &wNmosL);
		
		areaUnit = (hNmosL * wNmosL) * 4 + (hNmosS * wNmosS) * 8 + (hNmosM * wNmosM) * 4;
	}
}

void CurrentSenseAmp::CalculateArea(double widthArray) {	// adjust CurrentSenseAmp area by fixing S/A width
	if (!initialized) {
		cout << "[CurrentSenseAmp] Error: Require initialization first!" << endl;
	} else {
		double x = sqrt(areaUnit/HEIGHT_WIDTH_RATIO_LIMIT); // area = HEIGHT_WIDTH_RATIO_LIMIT * x^2
		if (widthArray > x)   // Limit W/H <= HEIGHT_WIDTH_RATIO_LIMIT
			widthArray = x;
		
		area = areaUnit * numCol;
		width = widthArray;
		height = area/widthArray;
		
	}
}


void CurrentSenseAmp::CalculateLatency(double numRead) {
	if (!initialized) {
		cout << "[Current Sense Amp] Error: Require initialization first!" << endl;
	} else {
		if (parallel) {
			readLatency = 1e-9;
		}else {
		    readLatency = 5e-9;
		}
		readLatency *= numRead;
	}
}

void CurrentSenseAmp::CalculatePower(double numRead) {
	if (!initialized) {
		cout << "[Current Sense Amp] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		if (parallel) {
			readDynamicEnergy = 52.01e-15;
		}else {
		    readDynamicEnergy = 181.3e-15;
		}
		readDynamicEnergy *= MIN(numReadCellPerOperationNeuro, numCol) * numRead;
		if (!readLatency) {
			//cout << "[CurrentSenseAmp] Error: Need to calculate read latency first" << endl;
		} else {
			readPower = readDynamicEnergy/readLatency;
		}
	}
}

void CurrentSenseAmp::PrintProperty(const char* str) {
	//cout << "Current Sense Amplifier Properties:" << endl;
	FunctionUnit::PrintProperty(str);
}

