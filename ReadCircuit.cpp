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
#include "ReadCircuit.h"

using namespace std;

ReadCircuit::ReadCircuit(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit() {
	initialized = false;
}

void ReadCircuit::Initialize(ReadCircuitMode _mode, int _numReadCol, int _maxNumIntBit, SpikingMode _spikingMode, double _clkFreq) {
	if (initialized)
		cout << "[ReadCircuit] Warning: Already initialized!" << endl;
	
	mode = _mode;
	maxNumIntPerCycle = (int)pow(2, maxNumIntBit);
	voltageIntThreshold = cell.readVoltage * RATIO_READ_THRESHOLD_VS_VOLTAGE;
	numReadCol = _numReadCol;
	maxNumIntBit = _maxNumIntBit;
	spikingMode = _spikingMode;
	clkFreq = _clkFreq;

	// Oscillation neron parameters
	Vhold = 0.5;	// Hold voltage
	Vth = 0.7;		// Threshold voltage
	R_OSC_OFF = 1e6;	// NbO2 Roff
	
	if (mode == CMOS) {
		Vrow = 0.65;					// Read voltage applied at rows
		Vcol = Vrow - cell.readVoltage;	// Column offset voltage
		if (spikingMode == NONSPIKING)	// Actually the original design of read circuit should be spiking (count spikes)
			numDff = maxNumIntBit;
		else
			numDff = (int)pow(2, maxNumIntBit);
	} else {	// OSCILLATION
		Vrow = 1.2;
		Vcol = 1.2 - cell.readVoltage;
		if (spikingMode == NONSPIKING)
			numDff = maxNumIntBit;
		else
			numDff = (int)pow(2, maxNumIntBit);
	}

	// Special DFF with RESET function
	widthDffTgN = MIN_NMOS_SIZE * tech.featureSize;
	widthDffTgP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	widthDffInvN = MIN_NMOS_SIZE * tech.featureSize;
	widthDffInvP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	widthDffNorN = MIN_NMOS_SIZE * tech.featureSize;
	widthDffNorP = 2 * tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;

	// Read circuit body
	widthTgN = 33.33 * tech.featureSize;
	widthTgP = 33.33 * tech.featureSize;
	widthNmos1 = 5 * tech.featureSize;
	widthPmos1 = 9 * tech.featureSize;
	widthNmos2 = 5 * tech.featureSize;
	widthNmos3 = 3.33 * tech.featureSize;
	widthPmos3 = 6.66 * tech.featureSize;
	widthNmos4 = 3.33 * tech.featureSize;
	widthPmos4 = 3.33 * tech.featureSize;
	widthNmos5 = 3.33 * tech.featureSize;
	widthPmos5 = 3.33 * tech.featureSize;
	widthNmos6 = 3.33 * tech.featureSize;
	widthNmos7 = 6.11 * tech.featureSize;
	widthNmos8 = 8.33 * tech.featureSize;
	widthPmos8 = 9 * tech.featureSize;
	// Buffer
	widthInvN = MIN_NMOS_SIZE * tech.featureSize;
	widthInvP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;

	initialized = true;
}

void ReadCircuit::CalculateUnitArea() {
	if (!initialized) {
		cout << "[ReadCircuit] Error: Require initialization first!" << endl;
	} else {
		
		double hTg, wTg, h1, w1, h2, w2, h3, w3, h4, w4, h5, w5, h6, w6, h7, w7, h8, w8, hBufInv, wBufInv;
		// Read circuit body
		if (mode == CMOS) {
			// Analytical estimation
			//CalculateGateArea(INV, 1, widthTgN, widthTgP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hTg, &wTg);
			//CalculateGateArea(INV, 1, widthNmos1, widthPmos1, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &h1, &w1);
			//CalculateGateArea(INV, 1, widthNmos2, 0, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &h2, &w2);
			//CalculateGateArea(INV, 1, widthNmos3, widthPmos3, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &h3, &w3);
			//CalculateGateArea(INV, 1, widthNmos4, widthPmos4, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &h4, &w4);
			//CalculateGateArea(INV, 1, widthNmos5, widthPmos5, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &h5, &w5);
			//CalculateGateArea(INV, 1, widthNmos6, 0, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &h6, &w6);
			//CalculateGateArea(INV, 1, widthNmos7, 0, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &h7, &w7);
			//CalculateGateArea(INV, 1, widthNmos8, widthPmos8, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &h8, &w8);
			//CalculateGateArea(INV, 1, widthInvN, widthInvP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hBufInv, &wBufInv);
			//wReadBody = wTg + w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + wBufInv * 2;
			//hReadBody = hTg;
			
			// Just use the dimension from layout
			hReadBody = MAX_TRANSISTOR_HEIGHT * tech.featureSize;
			wReadBody = 20 * (POLY_WIDTH + MIN_GAP_BET_GATE_POLY) * tech.featureSize;   // Need 21 polys

		} else {	// mode==OSCILLATION, only one INV
			CalculateGateArea(INV, 1, widthInvN, widthInvP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hBufInv, &wBufInv);
			wReadBody = wBufInv;
			hReadBody = hBufInv;
		}
		areaReadBody = hReadBody * wReadBody;

		// Assume DFF size is 12 minimum-size standard cells put together
		double hMinCell, wMinCell;
		CalculateGateArea(INV, 1, MIN_NMOS_SIZE * tech.featureSize, tech.pnSizeRatio*MIN_NMOS_SIZE * tech.featureSize, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hMinCell, &wMinCell);
		hDff = hMinCell;
		wDff = wMinCell * 12;
		areaDff = hDff * wDff;
		
		hUnit = hReadBody + hDff * numDff;
		wUnit = MAX(wReadBody, wDff);
		areaUnit = hUnit * wUnit;

		// Capacitance
		// DFF
		capDffTgGateN = CalculateGateCap(widthDffTgN, tech);
		capDffTgGateP = CalculateGateCap(widthDffTgP, tech);
		CalculateGateCapacitance(INV, 1, widthDffTgN, widthDffTgP, hMinCell, tech, NULL, &capDffTgDrain);
		CalculateGateCapacitance(INV, 1, widthDffInvN, widthDffInvP, hMinCell, tech, &capDffInvInput, &capDffInvOutput);
		CalculateGateCapacitance(NOR, 2, widthDffNorN, widthDffNorP, hMinCell, tech, &capNorInput, &capNorOutput);
		// Read circuit body
		capTgGateN = CalculateGateCap(widthTgN, tech);
		capTgGateP = CalculateGateCap(widthTgP, tech);
		CalculateGateCapacitance(INV, 1, widthTgN, widthTgP, hTg, tech, NULL, &capTgDrain);
		CalculateGateCapacitance(INV, 1, widthInvN, widthInvP, hBufInv, tech, &capInvInput, &capInvOutput);
		CalculateGateCapacitance(INV, 1, widthPmos4, 0, h4, tech, &capPmosGate, &capPmosDrain);
	}
}

void ReadCircuit::CalculateArea(double _newWidth) {	// Just add up the area of all the components
	if (!initialized) {
		cout << "[ReadCircuit] Error: Require initialization first!" << endl;
	} else {
		if (_newWidth) {	// Need to put read circuit units in multiple rows if there is a limitation on total width
			if (wUnit > _newWidth) {
				cout << "[ReadCircuit] Error: width too small even for 1 read circuit" << endl;
			}
			numUnitPerRow = (int)(_newWidth/wUnit);
			if (numUnitPerRow > numReadCol) {
				numUnitPerRow = numReadCol;
			}
			numRowUnit = (int)ceil((double)numReadCol / numUnitPerRow);
			width = _newWidth;
		} else {	// Just one row if there is no limitation on total width
			numUnitPerRow = numReadCol;
			numRowUnit = 1;
			width = wUnit * numUnitPerRow;
		}

		height = hUnit * numRowUnit;
		area = height * width;
	}
}

void ReadCircuit::CalculateLatency(double numRead) {
	if (!initialized) {
		cout << "[ReadCircuit] Error: Require initialization first!" << endl;
	} else {
		readLatency = 0;
		
		readLatency += cell.readPulseWidth;
		readLatency += 1/clkFreq;
		
		readLatency *= numRead;
	}
}

void ReadCircuit::CalculatePower(double numRead) {
	if (!initialized) {
		cout << "[ReadCircuit] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;
		
		// Leakage (DFF) (rough calculation)
		leakage += CalculateGateLeakage(INV, 1, widthDffInvN, widthDffInvP, inputParameter.temperature, tech) * tech.vdd * 8 * numDff;
		// Leakage (Read circuit body)
		if (mode == CMOS) {
			// Analytical result
			leakage += CalculateGateLeakage(INV, 1, widthNmos1, widthPmos1, inputParameter.temperature, tech) * tech.vdd;
			leakage += CalculateGateLeakage(INV, 1, widthNmos2, 0, inputParameter.temperature, tech) * tech.vdd;
			leakage += CalculateGateLeakage(INV, 1, widthNmos3, widthPmos3, inputParameter.temperature, tech) * tech.vdd;
			leakage += CalculateGateLeakage(INV, 1, widthNmos4, widthPmos4, inputParameter.temperature, tech) * tech.vdd;
			leakage += CalculateGateLeakage(INV, 1, widthNmos5, widthPmos5, inputParameter.temperature, tech) * tech.vdd;
			leakage += CalculateGateLeakage(INV, 1, widthNmos6, 0, inputParameter.temperature, tech) * tech.vdd;
			leakage += CalculateGateLeakage(INV, 1, widthNmos7, 0, inputParameter.temperature, tech) * tech.vdd;
			leakage += CalculateGateLeakage(INV, 1, widthNmos8, widthPmos8, inputParameter.temperature, tech) * tech.vdd;
			// Buffer
			leakage += CalculateGateLeakage(INV, 1, widthInvN, widthInvP, inputParameter.temperature, tech) * tech.vdd * 2;

			// SPICE result	(65nm tech node)
			//leakage += 104.9e-6;

		} else {	// mode==OSCILLATION, only one INV
			// Analytical result
			//leakage += CalculateGateLeakage(INV, 1, widthInvN, widthInvP, inputParameter.temperature, tech) * tech.vdd;

			// SPICE result (65nm tech node)
			leakage += 35.84e-9;
		}

		leakage *= numReadCol;

		// Dynamic energy (currently just import values)
		if (mode == CMOS) {
			readDynamicEnergy = 1.346e-12;	// 65nm tech node
		} else {
			readDynamicEnergy = 0.265e-12;	// 65nm tech node
		}
		readDynamicEnergy *= numReadCol;
		readDynamicEnergy *= numRead;

		if (!readLatency) {
			//cout << "[ReadCircuit] Error: Need to calculate read latency first" << endl;
		} else {
			readPower = readDynamicEnergy/readLatency;
		}
	}
}

void ReadCircuit::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

