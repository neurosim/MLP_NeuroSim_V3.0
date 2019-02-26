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
#include "WLDecoderOutput.h"

using namespace std;

WLDecoderOutput::WLDecoderOutput(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit(){
	initialized = false;
}

void WLDecoderOutput::Initialize(int _numWLRow) {
	if (initialized)
		cout << "[WLDecoderOutput] Warning: Already initialized!" << endl;

	numWLRow = _numWLRow;

	// NOR2
	widthNorN = MIN_NMOS_SIZE * tech.featureSize;
    widthNorP = 2 * tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	
	// INV
	widthInvN = MIN_NMOS_SIZE * tech.featureSize;
	widthInvP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	
	// TG
	widthTgN = MIN_NMOS_SIZE * tech.featureSize;
	widthTgP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;

	// NMOS
	widthNmos = MIN_NMOS_SIZE * tech.featureSize;

	initialized = true;
}

void WLDecoderOutput::CalculateArea(double _newHeight, double _newWidth, AreaModify _option) {
	if (!initialized) {
		cout << "[WLDecoderOutput] Error: Require initialization first!" << endl;
	} else {
		double hNor, wNor, hInv, wInv, hTg, wTg, hNmos, wNmos;
		double minCellHeight = MAX_TRANSISTOR_HEIGHT * tech.featureSize;

		// NOR2
		CalculateGateArea(NOR, 2, widthNorN, widthNorP, minCellHeight, tech, &hNor, &wNor);
		// INV
		CalculateGateArea(INV, 1, widthInvN, widthInvP, minCellHeight, tech, &hInv, &wInv);
		// TG
		CalculateGateArea(INV, 1, widthTgN, widthTgP, minCellHeight, tech, &hTg, &wTg);
		// NMOS
		CalculateGateArea(INV, 1, widthNmos, 0, minCellHeight, tech, &hNmos, &wNmos);

		if (_newHeight && _option==NONE) {
			int numColTg = 1;
			
			// Adjust TG height based on _newHeight if necessary
			double TgHeight = _newHeight / numWLRow;
			if (TgHeight < minCellHeight) {
				numColTg = (int)ceil(minCellHeight/TgHeight); // Calculate the number of columns
				if (numColTg > numWLRow) {
					cout << "[WLDecoderOutput] Error: pass gate height is even larger than the array height" << endl;
				}
				TgHeight = _newHeight / (int)ceil((double)numWLRow/numColTg);
			}
			CalculateGateArea(INV, 1, widthTgN, widthTgP, TgHeight, tech, &hTg, &wTg);
			
			double hUnit = MAX(MAX(hInv,hNor), MAX(hNmos,hTg));
			double wUnit = wNor + wInv + wNmos + wTg;
			int numUnitPerCol = (int)(_newHeight / hUnit);
			int numColUnit = (int)ceil((double)numWLRow / numUnitPerCol);
			if (numUnitPerCol > numWLRow) {
				numUnitPerCol = numWLRow;
			}
			height = _newHeight;
			width = (wNor + wInv + wNmos + wTg) * numColUnit;

		} else {
			height = MAX(MAX(hInv,hNor), MAX(hTg,hNmos)) * numWLRow;
			width = wNor + wInv + wTg + wNmos;
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

		// Resistance
		// TG
		double resTgN, resTgP;
		resTgN = CalculateOnResistance(widthTgN, NMOS, inputParameter.temperature, tech);
		resTgP = CalculateOnResistance(widthTgP, PMOS, inputParameter.temperature, tech);
		resTg = 1/(1/resTgN + 1/resTgP);

		// Capacitance
		// NOR2
		CalculateGateCapacitance(NOR, 2, widthNorN, widthNorP, hNor, tech, &capNorInput, &capNorOutput);
		// INV
		CalculateGateCapacitance(INV, 1, widthInvN, widthInvP, hInv, tech, &capInvInput, &capInvOutput);
		// TG
		capTgGateN = CalculateGateCap(widthTgN, tech);
		capTgGateP = CalculateGateCap(widthTgP, tech);
		CalculateGateCapacitance(INV, 1, widthTgN, widthTgP, hTg, tech, NULL, &capTgDrain);
		// NMOS
		CalculateGateCapacitance(INV, 1, widthNmos, 0, hNmos, tech, &capNmosGate, &capNmosDrain);
	}
}

void WLDecoderOutput::CalculateLatency(double _rampInput, double _capLoad, double _resLoad, double numRead, double numWrite) {
	if (!initialized) {
		cout << "[WLDecoderOutput] Error: Require initialization first!" << endl;
	} else {
		readLatency = 0;
		writeLatency = 0;
		
		rampInput = _rampInput;
		capLoad = _capLoad;
		resLoad = _resLoad;
		double resPullUp;
		double capOutput;
		double tr;		/* time constant */
		double gm;		/* transconductance */
		double beta;	/* for horowitz calculation */
		
		// 1st stage: NOR2
		resPullUp = CalculateOnResistance(widthNorP, PMOS, inputParameter.temperature, tech) * 2;
		tr = resPullUp * (capNorOutput + capInvInput + capTgGateP + capNmosGate);
		gm = CalculateTransconductance(widthNorP, PMOS, tech);
		beta = 1 / (resPullUp * gm);
		readLatency += horowitz(tr, beta, rampInput, NULL);
		
		// TG delay
		capOutput = capTgDrain + capNmosDrain;
		tr = resTg * (capOutput + capLoad) + resLoad * capLoad / 2;
		readLatency += horowitz(tr, 0, 1e20, &rampOutput);	// get from chargeLatency in the original SubArray.cpp
		
		readLatency *= numRead;
		writeLatency = readLatency / numRead * numWrite;
	}
}

void WLDecoderOutput::CalculatePower(double numRead, double numWrite) {
	if (!initialized) {
		cout << "[WLDecoderOutput] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;
		writeDynamicEnergy = 0;

		// Leakage power
		// NOR2
		leakage += CalculateGateLeakage(NOR, 2, widthNorN, widthNorP, inputParameter.temperature, tech) * tech.vdd * numWLRow;
		// INV
		leakage += CalculateGateLeakage(INV, 1, widthInvN, widthInvP, inputParameter.temperature, tech) * tech.vdd * numWLRow;
		// NMOS
		leakage += CalculateGateLeakage(INV, 1, widthNmos, 0, inputParameter.temperature, tech) * tech.vdd * numWLRow;
		
		// Read dynamic energy
		readDynamicEnergy += capNorInput * tech.vdd * tech.vdd * numWLRow;
		readDynamicEnergy += (capInvOutput + capTgGateN) * tech.vdd * tech.vdd * (numWLRow-1);
		readDynamicEnergy += capTgDrain * cell.accessVoltage * cell.accessVoltage * (numWLRow-1);
		readDynamicEnergy *= numRead;
		if (!readLatency) {
			//cout << "[WLDecoderOutput] Error: Need to calculate read latency first" << endl;
		} else {
			readPower = readDynamicEnergy/readLatency;
		}
		
		// Write dynamic energy (only one row activated)
		writeDynamicEnergy += capNorInput * tech.vdd * tech.vdd;
		writeDynamicEnergy += (capInvOutput + capTgGateN) * tech.vdd * tech.vdd;
		writeDynamicEnergy += capTgDrain * cell.accessVoltage * cell.accessVoltage;
		writeDynamicEnergy *= numWrite;
		if (!writeLatency) {
			//cout << "[WLDecoderOutput] Error: Need to calculate write latency first" << endl;
		} else {
			writePower = writeDynamicEnergy/writeLatency;
		}

	}
}

void WLDecoderOutput::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

