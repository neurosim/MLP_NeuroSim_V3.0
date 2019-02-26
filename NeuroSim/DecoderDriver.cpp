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
#include "DecoderDriver.h"

using namespace std;

DecoderDriver::DecoderDriver(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit() {
	initialized = false;
}

void DecoderDriver::Initialize(int _mode, int _numOutput /* # of array rows/columns */, int numLoad) {
	if (initialized)
		cout << "[Decoder Driver] Warning: Already initialized!" << endl;

	mode = _mode;
	numOutput = _numOutput;
	
	// INV
	widthInvN = MIN_NMOS_SIZE * tech.featureSize;
	widthInvP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;

	// TG
	resTg = cell.resMemCellOn / numLoad * IR_DROP_TOLERANCE;
	widthTgN = CalculateOnResistance(tech.featureSize, NMOS, inputParameter.temperature, tech)
				* tech.featureSize / (resTg*2);
	widthTgP = CalculateOnResistance(tech.featureSize, PMOS, inputParameter.temperature, tech)
				* tech.featureSize / (resTg*2);;
	
	initialized = true;
}

void DecoderDriver::CalculateArea(double _newHeight, double _newWidth, AreaModify _option) {
	if (!initialized) {
		cout << "[Decoder Driver] Error: Require initialization first!" << endl;
	} else {
		double hInv, wInv, hTg, wTg;
		double minCellHeight = MAX_TRANSISTOR_HEIGHT * tech.featureSize;
		double minCellWidth = 2 * (POLY_WIDTH + MIN_GAP_BET_GATE_POLY) * tech.featureSize;
		
		// TG
		if (_newHeight && _option==NONE) {
			if (_newHeight < minCellHeight) {
				cout << "[DecoderDriver] Error: pass gate height is even larger than the array height" << endl;
			}

			int numTgPerCol = (int)(_newHeight / minCellHeight);    // Get max # Tg per column (this is not the final # Tg per column because the last column may have less # Tg)
			numColTg = (int)ceil((double)numOutput / numTgPerCol); // Get min # columns based on this max # Tg per column
			numTgPerCol = (int)ceil((double)numOutput / numColTg);     // Get # Tg per column based on this min # columns
			TgHeight = _newHeight / numTgPerCol;
			CalculateGateArea(INV, 1, widthTgN, widthTgP, TgHeight, tech, &hTg, &wTg);

		} else if (_newWidth && _option==NONE) {
			if (minCellWidth > _newWidth) {
				cout << "[DecoderDriver] Error: pass gate width is even larger than the array width" << endl;
			}

			int numTgPerRow = (int)(_newWidth / minCellWidth);    // Get max # Tg per row (this is not the final # Tg per row because the last row may have less # Tg)
			numRowTg = (int)ceil((double)numOutput / numTgPerRow); // Get min # rows based on this max # Tg per row
			numTgPerRow = (int)ceil((double)numOutput / numRowTg);     // Get # Tg per row based on this min # rows
			TgWidth = _newWidth / numTgPerRow;
			int numFold = (int)(TgWidth / (0.5*minCellWidth)) - 1;  // get the max number of folding
			CalculatePassGateArea(widthTgN, widthTgP, tech, numFold, &hTg, &wTg);

		} else {
			CalculateGateArea(INV, 1, widthTgN, widthTgP, minCellHeight, tech, &hTg, &wTg);
		}
		
		// INV
		CalculateGateArea(INV, 1, widthInvN, widthInvP, minCellHeight, tech, &hInv, &wInv);

		double hUnit, wUnit;
		if (cell.accessType == CMOS_access) {	// 1T1R
			if (mode == ROW_MODE) {
				hUnit = MAX(hInv, hTg);
				wUnit = wInv + wTg * 3;
			} else {
				hUnit = hInv + hTg * 3;
				wUnit = MAX(wInv, wTg);
			}
		} else {	// Cross-point
			if (mode == ROW_MODE) {
				hUnit = MAX(hInv, hTg);
				wUnit = wInv + wTg * 2;
			} else {
				hUnit = hInv + hTg * 2;
				wUnit = MAX(wInv, wTg);
			}
		}

		if (mode == ROW_MODE) {	// Connect to rows
			if (_newHeight && _option==NONE) {
				int numColUnit, numUnitPerCol;
				numUnitPerCol = (int)(_newHeight/hUnit);
				numColUnit = (int)ceil((double)numOutput/numUnitPerCol);
				if (numColUnit > numOutput) {
					numColUnit = numOutput;
				}
				height = _newHeight;
				width = wUnit * numColUnit;
			} else {
				height = hUnit * numOutput;
				width = wUnit;
			}
		} else {	// Connect to columns
			if (_newWidth && _option==NONE) {
				int numRowUnit, numUnitPerRow;
				numUnitPerRow = (int)(_newWidth/wUnit);
				numRowUnit = (int)ceil((double)numOutput/numUnitPerRow);
				if (numRowUnit > numOutput) {
					numRowUnit = numOutput;
				}
				height = hUnit * numRowUnit;
				width = _newWidth;
			} else {
				height = hUnit;
				width = wUnit * numOutput;
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
		// INV
		CalculateGateCapacitance(INV, 1, widthInvN, widthInvP, hInv, tech, &capInvInput, &capInvOutput);
		// TG
		capTgGateN = CalculateGateCap(widthTgN, tech);
		capTgGateP = CalculateGateCap(widthTgP, tech);
		CalculateGateCapacitance(INV, 1, widthTgN, widthTgP, hTg, tech, NULL, &capTgDrain);
	}
}

void DecoderDriver::CalculateLatency(double _rampInput, double _capLoad1, double _capLoad2, double _resLoad, double numRead, double numWrite) {
	if (!initialized) {
		cout << "[Decoder Driver] Error: Require initialization first!" << endl;
	} else {
		readLatency = 0;
		capLoad1 = _capLoad1;   // worst-case load (1T1R SL, crosspoint WL/BL)
		capLoad2 = _capLoad2;   // 1T1R BL (which does not include transistor drain cap), doesn't matter for crosspoint
		resLoad = _resLoad;
		
		rampInput = _rampInput;
		double capOutput;
		double tr;	/* time constant */
		double gm;	/* transconductance */
		double beta;	/* for horowitz calculation */
		
		// TG
		capOutput = capTgDrain * 2;
		tr = resTg * (capOutput + capLoad1) + resLoad * capLoad1 / 2;
		readLatency += horowitz(tr, 0, rampInput, &rampOutput); // get from chargeLatency in the original SubArray.cpp
		readLatency *= numRead;

		writeLatency = cell.writePulseWidth;
		writeLatency *= numWrite;
	}
}

void DecoderDriver::CalculatePower(double numReadCellPerOp, double numWriteCellPerOp, double numRead, double numWrite) {
	if (!initialized) {
		cout << "[Decoder Driver] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;
		writeDynamicEnergy = 0;

		// Leakage power
		leakage += CalculateGateLeakage(INV, 1, widthInvN, widthInvP, inputParameter.temperature, tech) * tech.vdd * numOutput;
		
		// Read dynamic energy
		if (cell.accessType == CMOS_access) {  // 1T1R
			// Selected SLs and BLs are floating
			// Unselected SLs and BLs are GND
			readDynamicEnergy += (capInvInput + capTgGateN * 2 + capTgGateP) * tech.vdd * tech.vdd * numReadCellPerOp;
			readDynamicEnergy += (capInvOutput + capTgGateP * 2 + capTgGateN) * tech.vdd * tech.vdd * numReadCellPerOp;
		} else {	// Crosspoint
			// For WL decoder driver, the selected WLs are GND
			// For BL decoder driver, the selected BLs are floating
			// No matter which one, the unselected WLs/BLs are read voltage
			readDynamicEnergy += (capInvInput + capTgGateN + capTgGateP) * tech.vdd * tech.vdd * numReadCellPerOp;
			readDynamicEnergy += (capInvOutput + capTgGateP + capTgGateN) * tech.vdd * tech.vdd * numReadCellPerOp;
			readDynamicEnergy += (capTgDrain * 2) * cell.readVoltage * cell.readVoltage * (numOutput-numReadCellPerOp);
		}
		readDynamicEnergy *= numRead;
		if (!readLatency) {
			//cout << "[Decoder Driver] Error: Need to calculate read latency first" << endl;
		} else {
			readPower = readDynamicEnergy/readLatency;
		}
		
		// Write dynamic energy
		if (cell.accessType == CMOS_access) {  // 1T1R
			// Worst case: RESET operation (because SL cap is larger than BL cap)
			writeDynamicEnergy += (capInvInput + capTgGateN * 2 + capTgGateP) * tech.vdd * tech.vdd * numWriteCellPerOp;
			writeDynamicEnergy += (capInvOutput + capTgGateP * 2 + capTgGateN) * tech.vdd * tech.vdd * numWriteCellPerOp;
			writeDynamicEnergy += (capTgDrain * 2) * cell.writeVoltage * cell.writeVoltage * numWriteCellPerOp;
		} else {    // Crosspoint
			writeDynamicEnergy += (capInvInput + capTgGateN + capTgGateP) * tech.vdd * tech.vdd * numWriteCellPerOp;
			writeDynamicEnergy += (capInvOutput + capTgGateP + capTgGateN) * tech.vdd * tech.vdd * numWriteCellPerOp;
			if (mode == ROW_MODE) {	// Connects to rows
				writeDynamicEnergy += (capTgDrain * 2) * cell.writeVoltage * cell.writeVoltage * numWriteCellPerOp;
				writeDynamicEnergy += (capTgDrain * 2) * cell.writeVoltage/2 * cell.writeVoltage/2 * (numOutput-numWriteCellPerOp);
			} else {	// Connects to columns
				writeDynamicEnergy += (capTgDrain * 2) * cell.writeVoltage/2 * cell.writeVoltage/2 * (numOutput-numWriteCellPerOp);
			}
		}
		writeDynamicEnergy *= numWrite;
		if (!writeLatency) {
			//cout << "[Decoder Driver] Error: Need to calculate write latency first" << endl;
		} else {
			writePower = writeDynamicEnergy/writeLatency;
		}

	}
}

void DecoderDriver::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

