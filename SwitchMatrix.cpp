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

#include <cmath>
#include <iostream>
#include "constant.h"
#include "formula.h"
#include "SwitchMatrix.h"

using namespace std;

SwitchMatrix::SwitchMatrix(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), dff(_inputParameter, _tech, _cell), FunctionUnit() {
	initialized = false;
}

void SwitchMatrix::Initialize(int _mode, int _numOutput, double _resTg, double _activityRowRead, double _activityColWrite, int _numWriteCellPerOperationNeuro, double _numWritePulse, double _clkFreq){
	if (initialized)
		cout << "[SwitchMatrix] Warning: Already initialized!" << endl;
	
	mode = _mode;
	numOutput = _numOutput;
	activityRowRead = _activityRowRead;
	activityColWrite = _activityColWrite;
	numWriteCellPerOperationNeuro = _numWriteCellPerOperationNeuro;
	numWritePulse = _numWritePulse;
	clkFreq = _clkFreq;
    
	// DFF
	dff.Initialize(numOutput, clkFreq);

	// TG  resTg = cell.resMemCellOn / numLoad * IR_DROP_TOLERANCE;
	resTg = _resTg;      // given actual TG resistance
	
	// Why use pre-defined resTg? Becasue we want to define TG resistance according to loading and performance ...
	
	widthTgN = CalculateOnResistance(tech.featureSize, NMOS, inputParameter.temperature, tech) * tech.featureSize / (resTg*2);
	// R~(1/W), calculate actual TG width based on feature-sized TG resistance and given actual TG resistance 
	
	widthTgP = CalculateOnResistance(tech.featureSize, PMOS, inputParameter.temperature, tech) * tech.featureSize / (resTg*2);
	// assuming resTgN = resTgP, so resTgN = resTgP = 2*resTg (connected in parallel)

	initialized = true;
}

void SwitchMatrix::CalculateArea(double _newHeight, double _newWidth, AreaModify _option) {
	if (!initialized) {
		cout << "[SwitchMatrix] Error: Require initialization first!" << endl;
	} else {
		double hTg, wTg;
		
		if (mode == ROW_MODE) {	// Connect to rows
			double minCellHeight = MAX_TRANSISTOR_HEIGHT * tech.featureSize;

			if (_newHeight && _option==NONE) {
				if (_newHeight < minCellHeight) {
					cout << "[SwitchMatrix] Error: pass gate height is even larger than the array height" << endl;
				}
				int numTgPairPerCol = (int)(_newHeight / minCellHeight);// Get max # Tg pair per column (this is not the final # Tg pair per column because the last column may have less # Tg)
				numColTgPair = (int)ceil((double)numOutput / numTgPairPerCol);	// Get min # columns based on this max # Tg pair per column
				numTgPairPerCol = (int)ceil((double)numOutput / numColTgPair);	// Get # Tg pair per column based on this min # columns
				TgHeight = _newHeight / numTgPairPerCol;
				CalculateGateArea(INV, 1, widthTgN, widthTgP, TgHeight, tech, &hTg, &wTg);
				
				// DFF
				dff.CalculateArea(_newHeight, NULL, NONE);
				
				height = _newHeight;
				width = (wTg * 2) * numColTgPair + dff.width;

			} else {
				CalculateGateArea(INV, 1, widthTgN, widthTgP, minCellHeight, tech, &hTg, &wTg); // Pass gate with folding
				height = hTg * numOutput;
				dff.CalculateArea(height, NULL, NONE);	// Need to give the height information, otherwise by default the area calculation of DFF is in column mode
				width = (wTg * 2) + dff.width;
			}
		} else {	// Connect to columns
			if (_newWidth && _option==NONE) {
				numRowTgPair = 1;
				double minCellWidth = 2 * (POLY_WIDTH + MIN_GAP_BET_GATE_POLY) * tech.featureSize; // min standard cell width for 1 Tg
				if (minCellWidth > _newWidth) {
					cout << "[SwitchMatrix] Error: pass gate width is even larger than the array width" << endl;
				}

				int numTgPairPerRow = (int)(_newWidth / (minCellWidth*2));    // Get max # Tg pair per row (this is not the final # Tg pair per row because the last row may have less # Tg)
				numRowTgPair = (int)ceil((double)numOutput / numTgPairPerRow); // Get min # rows based on this max # Tg pair per row
				numTgPairPerRow = (int)ceil((double)numOutput / numRowTgPair);     // Get # Tg pair per row based on this min # rows
				TgWidth = _newWidth / numTgPairPerRow / 2;	// division of 2 because there are 2 Tg in one pair
				int numFold = (int)(TgWidth / (0.5*minCellWidth)) - 1;  // get the max number of folding

				// widthTgN, widthTgP and numFold can determine the height and width of each pass gate
				CalculatePassGateArea(widthTgN, widthTgP, tech, numFold, &hTg, &wTg);
				
				// DFF
				dff.CalculateArea(NULL, _newWidth, NONE);

				width = _newWidth;
				height = hTg * numRowTgPair + dff.height;

			} else {
				// Default (pass gate with folding=1)
				CalculatePassGateArea(widthTgN, widthTgP, tech, 1, &hTg, &wTg);
				width = wTg * 2 * numOutput;
				dff.CalculateArea(NULL, width, NONE);
				height = hTg + dff.height;
			}
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
		// TG
		capTgGateN = CalculateGateCap(widthTgN, tech);
		capTgGateP = CalculateGateCap(widthTgP, tech);
		CalculateGateCapacitance(INV, 1, widthTgN, widthTgP, hTg, tech, NULL, &capTgDrain);
		
	}
}

void SwitchMatrix::CalculateLatency(double _rampInput, double _capLoad, double _resLoad, double numRead, double numWrite) {	// For simplicity, assume shift register is ideal
	if (!initialized) {
		cout << "[SwitchMatrix] Error: Require initialization first!" << endl;
	} else {
		rampInput = _rampInput;
		capLoad = _capLoad;
		resLoad = _resLoad;
		double capOutput;
		double tr;  /* time constant */
		readLatency = 0;

		// DFF
		dff.CalculateLatency(1e20, numRead);

		// TG
		capOutput = capTgDrain * 3;
		tr = resTg * (capOutput + capLoad) + resLoad * capLoad / 2;
		readLatency += horowitz(tr, 0, rampInput, &rampOutput);	// get from chargeLatency in the original SubArray.cpp
		
		readLatency *= numRead;
		readLatency += dff.readLatency;

		writeLatency = cell.writePulseWidth;
		writeLatency *= numWrite;
		writeLatency += dff.readLatency;	// Use DFF read latency here because no write in the DFF module
	}
}

void SwitchMatrix::CalculatePower(double numRead, double numWrite) {
	if (!initialized) {
		cout << "[SwitchMatrix] Error: Require initialization first!" << endl;
	} else {
		
		leakage = 0;
		readDynamicEnergy = 0;
		writeDynamicEnergy = 0;
		
		// DFF
		dff.CalculatePower(numRead, numOutput);	// Use numOutput since every DFF will pass signal (either 0 or 1)

		// Leakage power
		leakage += dff.leakage;	// Only DFF has leakage

		// Read dynamic energy
		if (mode == ROW_MODE) {
			readDynamicEnergy += (capTgDrain * 3) * cell.readVoltage * cell.readVoltage * numOutput * activityRowRead;
			readDynamicEnergy += (capTgGateN + capTgGateP) * tech.vdd * tech.vdd * numOutput * activityRowRead;
		} // No read energy in COL_MODE
		readDynamicEnergy *= numRead;
		readDynamicEnergy += dff.readDynamicEnergy;
		
		if (!readLatency) {
			//cout << "[Switch Matrix] Error: Need to calculate read latency first" << endl;
		} else {
			readPower = readDynamicEnergy/readLatency;
		}
		
		// Write dynamic energy (2-step write and average case half SET and half RESET)
		if (cell.accessType == CMOS_access) {	// 1T1R

			if (mode == ROW_MODE) {	// Connects to rows
				writeDynamicEnergy += (capTgGateN + capTgGateP) * tech.vdd * tech.vdd * 2;	// Selected row in LTP, *2 means switching from one selected row to another
			} else {	// Connects to columns
				// LTP
				writeDynamicEnergy += (capTgDrain * 3) * cell.writeVoltage * cell.writeVoltage * numWritePulse * MIN(numWriteCellPerOperationNeuro, numOutput*activityColWrite) / 2;   // Selected columns
				writeDynamicEnergy += (capTgDrain * 3) * cell.writeVoltage * cell.writeVoltage * (numOutput - MIN(numWriteCellPerOperationNeuro, numOutput*activityColWrite)/2);   // Unselected columns 
				// LTD
				writeDynamicEnergy += (capTgDrain * 3) * cell.writeVoltage * cell.writeVoltage * numWritePulse * MIN(numWriteCellPerOperationNeuro, numOutput*activityColWrite) / 2;   // Selected columns
				
				writeDynamicEnergy += (capTgGateN + capTgGateP) * tech.vdd * tech.vdd * numOutput;
			}

		} else {	// Cross-point
			
			if (mode == ROW_MODE) { // Connects to rows
				writeDynamicEnergy += (capTgDrain * 3) * cell.writeVoltage * cell.writeVoltage;   // Selected row in LTP
				writeDynamicEnergy += (capTgDrain * 3) * cell.writeVoltage/2 * cell.writeVoltage/2 * (numOutput-1);   // Unselected rows in LTP and LTD
				writeDynamicEnergy += (capTgGateN + capTgGateP) * tech.vdd * tech.vdd * numOutput;
			} else {    // Connects to columns
				writeDynamicEnergy += (capTgDrain * 3) * cell.writeVoltage * cell.writeVoltage * numWritePulse * MIN(numWriteCellPerOperationNeuro, numOutput*activityColWrite) / 2;   // Selected columns in LTP
				writeDynamicEnergy += (capTgDrain * 3) * cell.writeVoltage * cell.writeVoltage * numWritePulse * MIN(numWriteCellPerOperationNeuro, numOutput*activityColWrite) / 2;   // Selected columns in LTD
				writeDynamicEnergy += (capTgDrain * 3) * cell.writeVoltage/2 * cell.writeVoltage/2 * numOutput;   // Total unselected columns in LTP and LTD within the 2-step write
				writeDynamicEnergy += (capTgGateN + capTgGateP) * tech.vdd * tech.vdd * numOutput;
			}

		}
		writeDynamicEnergy *= numWrite;
		writeDynamicEnergy += dff.readDynamicEnergy;	// Use DFF read energy here because no write in the DFF module
		if (!writeLatency) {
			//cout << "[Switch Matrix] Error: Need to calculate write latency first" << endl;
		} else {
			writePower = writeDynamicEnergy/writeLatency;
		}
	}
}

void SwitchMatrix::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

