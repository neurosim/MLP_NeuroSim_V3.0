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
#include "RowDecoder.h"

using namespace std;

RowDecoder::RowDecoder(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit(){
	initialized = false;
}

void RowDecoder::Initialize(DecoderMode _mode, int _numAddrRow, bool _MUX) {
	if (initialized)
		cout << "[Row Decoder] Warning: Already initialized!" << endl;
	
	mode = _mode;
	numAddrRow = _numAddrRow;
	MUX = _MUX;

	// Use 2-bit predecoding
	// INV
	widthInvN = MIN_NMOS_SIZE * tech.featureSize;
	widthInvP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	numInv = numAddrRow;	// The INV at outpur driver stage does not count here

	// NAND2
	widthNandN = 2 * MIN_NMOS_SIZE * tech.featureSize;
	widthNandP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	numNand = 4 * (int)(floor(numAddrRow/2));

	// NOR (ceil(N/2) inputs)
	widthNorN = MIN_NMOS_SIZE * tech.featureSize;
	widthNorP = (int)ceil((double)numAddrRow/2) * tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	if (numAddrRow > 2)
		numNor = pow(2, numAddrRow);
	else
		numNor = 0;

	// Number of M3 for connection between NAND2 and NOR stages (if numAddrRow > 2)
	if (numAddrRow > 2)
		numMetalConnection = numNand + (numAddrRow%2) * 2;
	else
		numMetalConnection = 0;
	
	// Output driver INV
	widthDriverInvN = 3 * MIN_NMOS_SIZE * tech.featureSize;
	widthDriverInvP = 3 * tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;

	initialized = true;
}

void RowDecoder::CalculateArea(double _newHeight, double _newWidth, AreaModify _option) {
	if (!initialized) {
		cout << "[Row Decoder Area] Error: Require initialization first!" << endl;
	} else {
		double hInv, wInv, hNand, wNand, hNor, wNor, hDriverInv, wDriverInv;
		// INV
		CalculateGateArea(INV, 1, widthInvN, widthInvP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hInv, &wInv);
		// NAND2
		CalculateGateArea(NAND, 2, widthNandN, widthNandP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hNand, &wNand);
		// NOR (ceil(N/2) inputs)
		CalculateGateArea(NOR, (int)ceil((double)numAddrRow/2), widthNorN, widthNorP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hNor, &wNor);
		// Output Driver INV
		CalculateGateArea(INV, 1, widthDriverInvN, widthDriverInvP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hDriverInv, &wDriverInv);

		if (mode == REGULAR_ROW) {	// Connect to rows
			if (_newHeight && _option==NONE) {
				// NOR
				numColNor = 0;	// Number of columns of NOR
				numNorPerCol = (int)(_newHeight/hNor);
				if (numNorPerCol > numNor) {
					numNorPerCol = numNor;
				}
				if (numNorPerCol > 0) {	// Prevent division by 0
					numColNor = (int)ceil((double)numNor/numNorPerCol);
				}
				// NAND2
				numColNand = 0;  // Number of columns of NAND2
				numNandPerCol = (int)(_newHeight/hNand);
				if (numNandPerCol > numNand) {
					numNandPerCol = numNand;
				}
				if (numNandPerCol > 0) {	// Prevent division by 0
					numColNand = (int)ceil((double)numNand/numNandPerCol);
				}
				// INV
				numColInv = 0;  // Number of columns of INV
				numInvPerCol = (int)(_newHeight/hInv);
				if (numInvPerCol > numInv) {
					numInvPerCol = numInv;
				}
				numColInv = (int)ceil((double)numInv/numInvPerCol);
				if (numNorPerCol * numNandPerCol * numInvPerCol == 0 && numAddrRow > 2) { // If either of them is zero and not because of no NOR and no NAND2
					cout << "[Row Decoder] Error: logic gate height is even larger than the assigned height" << endl;
				}
				height = _newHeight;
				width = wInv * numColInv + wNand * numColNand + M3_PITCH * numMetalConnection * tech.featureSize + wNor * numColNor;
				if (MUX) {    // Mux enable circuit (NAND + INV) + INV
					width += (wNand + wInv * 2) * numColNor;
				} else {    // REGULAR: 2 INV as output driver
					width += (wDriverInv * 2) * numColNor;
				}

			} else {
				height = MAX(hNor*numNor, hNand*numNand);
				width = wInv + wNand + M3_PITCH * numMetalConnection * tech.featureSize + wNor;
				if (MUX) {	// Mux enable circuit (NAND + INV) + INV
					width += wNand + wInv * 2;
				} else {	// REGULAR: 2 INV as output driver
					width += wDriverInv * 2;
				}
			}
		} else {	// mode==REGULAR_COL
			if (_newWidth && _option==NONE) {
				// NOR
				numRowNor = 0;  // Number of rows of NOR
				numNorPerRow = (int)(_newWidth/wNor);
				if (numNorPerRow > numNor) {
					numNorPerRow = numNor;
				}
				if (numNorPerRow > 0) {	// Prevent division by 0
					numRowNor = (int)ceil((double)numNor/numNorPerRow);
				}
				// NAND2
				numRowNand = 0;  // Number of rows of NAND2
				numNandPerRow = (int)(_newWidth/wNand);
				if (numNandPerRow > numNand) {
					numNandPerRow = numNand;
				}
				if (numNandPerRow > 0) {	// Prevent division by 0
					numRowNand = (int)ceil((double)numNand/numNandPerRow);
				}
				// INV
				numRowInv = 0;  // Number of rows of INV
				numInvPerRow = (int)(_newWidth/wInv);
				if (numInvPerRow > numInv) {
					numInvPerRow = numInv;
				}
				numRowInv = (int)ceil((double)numInv/numInvPerRow);

				if (numNorPerRow * numNandPerRow * numInvPerRow == 0 && numAddrRow > 2) { // If either of them is zero and not because of no NOR and no NAND2
					cout << "[Row Decoder] Error: logic gate width is even larger than the assigned width" << endl;
				}
				width = _newWidth;
				height = hInv * numRowInv + hNand * numRowNand + M2_PITCH * numMetalConnection * tech.featureSize + hNor * numRowNor;
				if (MUX) {    // Mux enable circuit (NAND + INV) + INV
					height += (hNand + hInv * 2) * numRowNor;
				} else {    // REGULAR: 2 INV as output driver
					height += (hDriverInv * 2) * numRowNor;
				}
			} else {
				height = hInv + hNand + M2_PITCH * numMetalConnection * tech.featureSize + hNor;
				width = MAX(wNor*numNor, wNand*numNand);
				if (MUX) {    // Mux enable circuit (NAND + INV) + INV
					height += hNand + hInv * 2;
				} else {    // REGULAR: 2 INV as output driver
					height += hDriverInv * 2;
				}
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
			default:	// NONE
				break;
		}
		
		// Capacitance
		// INV
		CalculateGateCapacitance(INV, 1, widthInvN, widthInvP, hInv, tech, &capInvInput, &capInvOutput);
		// NAND2
		if (numNand) {
			CalculateGateCapacitance(NAND, 2, widthNandN, widthNandP, hNand, tech, &capNandInput, &capNandOutput);
		} else {
			capNandInput = capNandOutput = 0;
		}
		// NOR (ceil(N/2) inputs)
		if (numNor) {
			CalculateGateCapacitance(NOR, (int)ceil((double)numAddrRow/2), widthNorN, widthNorP, hNor, tech, &capNorInput, &capNorOutput);
		} else {
			capNorInput = capNorOutput = 0;
		}
		// Output Driver INV
		CalculateGateCapacitance(INV, 1, widthDriverInvN, widthDriverInvP, hDriverInv, tech, &capDriverInvInput, &capDriverInvOutput);
	}
}

void RowDecoder::CalculateLatency(double _rampInput, double _capLoad1, double _capLoad2, double numRead, double numWrite) {
	if (!initialized) {
		cout << "[Row Decoder Latency] Error: Require initialization first!" << endl;
	} else {
		rampInput = _rampInput;
		capLoad1 = _capLoad1;   // REGULAR: general capLoad, MUX: the NMOS Tg gates
		capLoad2 = _capLoad2;   // MUX: the PMOS Tg gates
		readLatency = 0;
		writeLatency = 0;

		double resPullDown, resPullUp;
		double tr;		/* time constant */
		double gm;		/* transconductance */
		double beta;	/* for horowitz calculation */
		double rampInvOutput = 1e20;
		double rampNandOutput = 1e20;
		double rampNorOutput = 1e20;
		
		// INV
		resPullDown = CalculateOnResistance(widthInvN, NMOS, inputParameter.temperature, tech);	// doesn't matter pullup/pulldown?
		if (numNand)
			tr = resPullDown * (capInvOutput + capNandInput * 2);	// one address line connects to 2 NAND inputs
		else
			tr = resPullDown * (capInvOutput + capLoad1);
		gm = CalculateTransconductance(widthInvN, NMOS, tech);
		beta = 1 / (resPullDown * gm);
		readLatency += horowitz(tr, beta, rampInput, &rampInvOutput);
		
		if (!numNand)
			rampOutput = rampInvOutput;

		// NAND2
		if (numNand) {
			resPullDown = CalculateOnResistance(widthNandN, NMOS, inputParameter.temperature, tech) * 2;
			if (numNor)
				tr = resPullDown * (capNandOutput + capNorInput * numNor/4);
			else
				tr = resPullDown * (capNandOutput + capLoad1);
			gm = CalculateTransconductance(widthNandN, NMOS, tech);
			beta = 1 / (resPullDown * gm);
			readLatency += horowitz(tr, beta, rampInvOutput, &rampNandOutput);
			if (!numNor)
				rampOutput = rampNandOutput;
		}
		
		// NOR (ceil(N/2) inputs)
		if (numNor) {
			resPullUp = CalculateOnResistance(widthNorP, PMOS, inputParameter.temperature, tech) * 2;
			if (MUX)
				tr = resPullUp * (capNorOutput + capNandInput);
			else
				tr = resPullUp * (capNorOutput + capInvInput);
			gm = CalculateTransconductance(widthNorP, PMOS, tech);
			beta = 1 / (resPullUp * gm);
			readLatency += horowitz(tr, beta, rampNandOutput, &rampNorOutput);
		}

		// Output driver or Mux enable circuit
		if (MUX) {	// Mux enable circuit (NAND + INV) + INV
			// 1st NAND
			resPullDown = CalculateOnResistance(widthNandN, NMOS, inputParameter.temperature, tech);
			tr = resPullDown * (capNandOutput + capInvInput);
			gm = CalculateTransconductance(widthNandN, NMOS, tech);
			beta = 1 / (resPullDown * gm);
			readLatency += horowitz(tr, beta, rampNorOutput, &rampNandOutput);
			// 2nd INV
			resPullUp = CalculateOnResistance(widthInvP, PMOS, inputParameter.temperature, tech);
			tr = resPullUp * (capInvOutput + capInvInput + capLoad1);
			gm = CalculateTransconductance(widthInvP, PMOS, tech);
			beta = 1 / (resPullUp * gm);
			readLatency += horowitz(tr, beta, rampNandOutput, &rampInvOutput);
			// 3rd INV
			resPullDown = CalculateOnResistance(widthInvN, NMOS, inputParameter.temperature, tech);
			tr = resPullDown * (capInvOutput + capLoad2);
			gm = CalculateTransconductance(widthInvN, NMOS, tech);
			beta = 1 / (resPullDown * gm);
			readLatency += horowitz(tr, beta, rampInvOutput, &rampOutput);
		} else {	// REGULAR: 2 INV as output driver
			// 1st INV
			resPullDown = CalculateOnResistance(widthDriverInvN, NMOS, inputParameter.temperature, tech);
			tr = resPullDown * (capDriverInvOutput + capDriverInvInput);
			gm = CalculateTransconductance(widthDriverInvN, NMOS, tech);
			beta = 1 / (resPullDown * gm);
			readLatency += horowitz(tr, beta, rampNorOutput, &rampInvOutput);
			// 2nd INV
			resPullUp = CalculateOnResistance(widthDriverInvP, PMOS, inputParameter.temperature, tech);
			tr = resPullUp * (capDriverInvOutput + capLoad1);
			gm = CalculateTransconductance(widthDriverInvP, PMOS, tech);
			beta = 1 / (resPullUp * gm);
			readLatency += horowitz(tr, beta, rampInvOutput, &rampOutput);
		}

		readLatency *= numRead;

		writeLatency = readLatency / numRead * numWrite;
	}
}

void RowDecoder::CalculatePower(double numRead, double numWrite) {
	if (!initialized) {
		cout << "[Row Decoder] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;
		writeDynamicEnergy = 0;
		// Leakage power
		// INV
		leakage += CalculateGateLeakage(INV, 1, widthInvN, widthInvP, inputParameter.temperature, tech) * tech.vdd * numInv;
		// NAND2
		leakage += CalculateGateLeakage(NAND, 2, widthNandN, widthNandP, inputParameter.temperature, tech) * tech.vdd * numNand;
		// NOR (ceil(N/2) inputs)
		leakage += CalculateGateLeakage(NOR, (int)ceil((double)numAddrRow/2), widthNorN, widthNorP, inputParameter.temperature, tech) * tech.vdd * numNor;
		// Output driver or Mux enable circuit
		if (MUX) {
			leakage += CalculateGateLeakage(NAND, 2, widthNandN, widthNandP, inputParameter.temperature, tech) * tech.vdd * numNor;
			leakage += CalculateGateLeakage(INV, 1, widthInvN, widthInvP, inputParameter.temperature, tech) * tech.vdd * 2 * numNor;
		} else {
			leakage += CalculateGateLeakage(INV, 1, widthDriverInvN, widthDriverInvP, inputParameter.temperature, tech) * tech.vdd * 2 * numNor;
		}

		// Read dynamic energy for both memory and neuro modes (rough calculation assuming all addr from 0 to 1)
		// INV
		readDynamicEnergy += (capInvInput + capNandInput * 2) * tech.vdd * tech.vdd * (int)floor(numAddrRow/2)*2;
		readDynamicEnergy += (capInvInput + capNorInput * numNor/2) * tech.vdd * tech.vdd * (numAddrRow - (int)floor(numAddrRow/2)*2);	// If numAddrRow is odd number
		// NAND2
		readDynamicEnergy += (capNandOutput + capNorInput * numNor/4) * tech.vdd * tech.vdd * numNand/4;	// every (NAND * 4) group has one NAND output activated
		// NOR (ceil(N/2) inputs)
		if (MUX)
			readDynamicEnergy += (capNorOutput + capNandInput) * tech.vdd * tech.vdd;	// one NOR output activated
		else
			readDynamicEnergy += (capNorOutput + capInvInput) * tech.vdd * tech.vdd;	// one NOR output activated
		// Output driver or Mux enable circuit
		if (MUX) {
			readDynamicEnergy += (capNandOutput + capInvInput) * tech.vdd * tech.vdd;
			readDynamicEnergy += (capInvOutput + capInvInput) * tech.vdd * tech.vdd;
			readDynamicEnergy += capInvOutput * tech.vdd * tech.vdd;
		} else {
			readDynamicEnergy += (capDriverInvInput + capDriverInvOutput) * tech.vdd * tech.vdd * 2;
		}

		readDynamicEnergy *= numRead;
		if (!readLatency) {
			//cout << "[Row Decoder] Error: Need to calculate read latency first" << endl;
		} else {
			readPower = readDynamicEnergy/readLatency;
		}
		
		// Write dynamic energy
		writeDynamicEnergy = readDynamicEnergy / numRead * numWrite;
		if (!writeLatency) {
			//cout << "[Row Decoder] Error: Need to calculate write latency first" << endl;
		} else {
			writePower = writeDynamicEnergy/writeLatency;
		}

	}
}

void RowDecoder::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

