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
#include "Precharger.h"

using namespace std;

Precharger::Precharger(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit() {
	initialized = false;
}

void Precharger::Initialize(int _numCol, double _resLoad, double _activityColWrite, int _numReadCellPerOperationNeuro, int _numWriteCellPerOperationNeuro) {
	if (initialized)
		cout << "[Precharger] Warning: Already initialized!" << endl;

	numCol = _numCol;
	resLoad = _resLoad;
	activityColWrite = _activityColWrite;
	numReadCellPerOperationNeuro = _numReadCellPerOperationNeuro;
	numWriteCellPerOperationNeuro = _numWriteCellPerOperationNeuro;
	
	widthPMOSBitlineEqual = MIN_NMOS_SIZE * tech.featureSize;
	widthPMOSBitlinePrecharger = 6 * tech.featureSize;
	
	initialized = true;
}

void Precharger::CalculateArea(double _newHeight, double _newWidth, AreaModify _option) {
	if (!initialized) {
		cout << "[Precharger] Error: Require initialization first!" << endl;
	} else {
		double hBitlinePrecharger, wBitlinePrecharger;
		double hBitlineEqual, wBitlineEqual;
		CalculateGateArea(INV, 1, 0, widthPMOSBitlinePrecharger, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hBitlinePrecharger, &wBitlinePrecharger);
		CalculateGateArea(INV, 1, 0, widthPMOSBitlineEqual, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hBitlineEqual, &wBitlineEqual);
		
		double hUnit = hBitlinePrecharger + hBitlineEqual * 2;
		double wUnit = MAX(wBitlinePrecharger, wBitlineEqual);

		if (_newWidth && _option==NONE) {
			int numRowUnit;  // Number of rows of unit
			int numUnitPerRow;
			numUnitPerRow = (int)(_newWidth/wUnit);
			if (numUnitPerRow > numCol) {
				numUnitPerRow = numCol;
			}
			numRowUnit = (int)ceil((double)numCol/numUnitPerRow);
			width = _newWidth;
			height = numRowUnit * hUnit;
		} else {
			width = numCol * wUnit;
			height = hUnit;
		}
		area = height * width;
		
		// Modify layout
		newHeight = _newHeight;
		newWidth = _newWidth;
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
		
		// Capacitance
		capOutputBitlinePrecharger = CalculateDrainCap(widthPMOSBitlinePrecharger, PMOS, hBitlinePrecharger, tech) + CalculateDrainCap(widthPMOSBitlineEqual, PMOS, hBitlineEqual, tech);
	}
}

void Precharger::CalculateLatency(double _rampInput, double _capLoad, double numRead, double numWrite){
	if (!initialized) {
		cout << "[Precharger] Error: Require initialization first!" << endl;
	} else {
		readLatency = 0;
		writeLatency = 0;

		rampInput = _rampInput;
		capLoad = _capLoad;
		double tr;		/* time constant */
		double gm;		/* transconductance */
		double beta;	/* for horowitz calculation */
		double resPullUp;
		double tau;

		resPullUp = CalculateOnResistance(widthPMOSBitlinePrecharger, PMOS, inputParameter.temperature, tech);
		tau = resPullUp * (capLoad + capOutputBitlinePrecharger) + resLoad * capLoad / 2;
		gm = CalculateTransconductance(widthPMOSBitlinePrecharger, PMOS, tech);
		beta = 1 / (resPullUp * gm);
		readLatency += horowitz(tau, beta, 1e20, &rampOutput);
		writeLatency = readLatency;

		readLatency *= numRead;
		writeLatency *= numWrite;
	}
}

void Precharger::CalculatePower(double numRead, double numWrite) {
	if (!initialized) {
		cout << "[Precharger] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;
		writeDynamicEnergy = 0;
		
		/* Leakage power */
		leakage = CalculateGateLeakage(INV, 1, 0, widthPMOSBitlinePrecharger, inputParameter.temperature, tech) * tech.vdd * numCol;

		/* Dynamic energy */
		// Read
		readDynamicEnergy = capLoad * tech.vdd * tech.vdd * MIN(numReadCellPerOperationNeuro, numCol) * 2;   // BL and BL_bar
		readDynamicEnergy *= numRead;
		
		// Write
		writeDynamicEnergy = capLoad * tech.vdd * tech.vdd * MIN(numWriteCellPerOperationNeuro, numCol*activityColWrite);
		writeDynamicEnergy *= numWrite;
	}
}

void Precharger::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

