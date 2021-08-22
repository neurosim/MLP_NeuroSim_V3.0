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
#include "SRAMWriteDriver.h"

using namespace std;

SRAMWriteDriver::SRAMWriteDriver(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit() {
	initialized = false;
}

void SRAMWriteDriver::Initialize(int _numCol, double _activityColWrite, int _numWriteCellPerOperationNeuro){
	if (initialized)
		cout << "[SRAMWriteDriver] Warning: Already initialized!" << endl;
	
	numCol = _numCol;
	activityColWrite = _activityColWrite;
	numWriteCellPerOperationNeuro = _numWriteCellPerOperationNeuro;
	
	widthInvN = MIN_NMOS_SIZE * tech.featureSize;
	widthInvP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	
	initialized = true;
}

void SRAMWriteDriver::CalculateArea(double _newHeight, double _newWidth, AreaModify _option) {
	if (!initialized) {
		cout << "[SRAMWriteDriver] Error: Require initialization first!" << endl;
	} else {
		double hInv, wInv;
		// INV
		CalculateGateArea(INV, 1, widthInvN, widthInvP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hInv, &wInv);

		double hUnit = hInv * 3;
		double wUnit = wInv;

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
		// INV
		CalculateGateCapacitance(INV, 1, widthInvN, widthInvP, hInv, tech, &capInvInput, &capInvOutput);
	}
}

void SRAMWriteDriver::CalculateLatency(double _rampInput, double _capLoad, double _resLoad, double numWrite){
	if (!initialized) {
		cout << "[SRAMWriteDriver] Error: Require initialization first!" << endl;
	} else {
		rampInput = _rampInput;
		capLoad = _capLoad;
		resLoad = _resLoad;
		writeLatency = 0;
		double tr;		/* time constant */
		double gm;		/* transconductance */
		double beta;	/* for horowitz calculation */
		double resPullUp, resPullDown;
		double rampInvOutput;
		
		// 1st stage INV (Pullup)
		resPullUp = CalculateOnResistance(widthInvP, PMOS, inputParameter.temperature, tech);
		tr = resPullUp * (capInvOutput + capInvInput);
		gm = CalculateTransconductance(widthInvP, PMOS, tech);
		beta = 1 / (resPullUp * gm);
		writeLatency += horowitz(tr, beta, rampInput, &rampInvOutput);
		
		// 2nd stage INV (Pulldown)
		resPullDown = CalculateOnResistance(widthInvN, NMOS, inputParameter.temperature, tech);
		tr = resPullDown * (capLoad + capInvOutput) + resLoad * capLoad / 2;
		gm = CalculateTransconductance(widthInvN, NMOS, tech);
		beta = 1 / (resPullDown * gm);
		writeLatency += horowitz(tr, beta, rampInvOutput, &rampOutput);

		writeLatency *= numWrite;
	}
}

void SRAMWriteDriver::CalculatePower(double numWrite) {
	if (!initialized) {
		cout << "[SRAMWriteDriver] Error: Require initialization first!" << endl;
	} else {
		/* Leakage power */
		leakage = CalculateGateLeakage(INV, 1, widthInvN, widthInvP, inputParameter.temperature, tech) * tech.vdd * 3 * numCol;

		/* Write Dynamic energy */
		// After the precharger precharges the BL and BL_bar to Vdd, the write driver only discharges one of them to zero, so there is no energy consumption on the BL and BL_bar
		// Assuming the write data is 0, the energy consumption is on the output of 1st INV and the input of 2nd INV
		writeDynamicEnergy = (capInvOutput + capInvInput) * tech.vdd * tech.vdd * MIN(numWriteCellPerOperationNeuro, numCol*activityColWrite);
		writeDynamicEnergy *= numWrite;
	}
}

void SRAMWriteDriver::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

