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
#include "ShiftAdd.h"

using namespace std;

ShiftAdd::ShiftAdd(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), adder(_inputParameter, _tech, _cell), dff(_inputParameter, _tech, _cell), FunctionUnit() {
	initialized = false;
}

void ShiftAdd::Initialize(int _numUnit, int _numAdderBit, double _clkFreq, SpikingMode _spikingMode, int _numReadPulse) {
	if (initialized)
		cout << "[ShiftAdd] Warning: Already initialized!" << endl;
	
	numUnit = _numUnit;
	numAdderBit = _numAdderBit;
	numAdder = numUnit;
	clkFreq = _clkFreq;
	spikingMode = _spikingMode;
	numReadPulse = _numReadPulse;
	
	if (spikingMode == NONSPIKING) {	// NONSPIKING: binary format
		numDff = (numAdderBit+1 + numReadPulse-1) * numUnit;	// numAdderBit+1 because the adder output is 1 bit more than the input, and numReadPulse-1 is for shift-and-add extension (shift register)
		dff.Initialize(numDff, clkFreq);
		adder.Initialize(numAdderBit, numAdder);
	} else {	// SPIKING: count spikes
		numBitPerDff = pow(2, numAdderBit);
		numDff = numBitPerDff * numUnit;	// numUnit shift registers in total
		dff.Initialize(numDff, clkFreq);
	}

	/* Currently ignore INV and NAND in shift-add circuit */
	// PISO shift register (https://en.wikipedia.org/wiki/Shift_register)
	// INV
	widthInvN = MIN_NMOS_SIZE * tech.featureSize;
	widthInvP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	numInv = numUnit;
	// NAND2
	widthNandN = 2 * MIN_NMOS_SIZE * tech.featureSize;
	widthNandP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	numNand = 3 * (numDff/numUnit-1) * numUnit;	// numDff/numUnit means the number of DFF for each shift register

	initialized = true;
}

void ShiftAdd::CalculateArea(double _newHeight, double _newWidth, AreaModify _option) {
	if (!initialized) {
		cout << "[ShiftAdd] Error: Require initialization first!" << endl;
	} else {
		double hInv, wInv, hNand, wNand;
		
		// Adder
		if (_newWidth && _option==NONE) {
			if (spikingMode == NONSPIKING) {   // NONSPIKING: binary format
				adder.CalculateArea(NULL, _newWidth, NONE);
				dff.CalculateArea(NULL, _newWidth, NONE);
			} else {    // SPIKING: count spikes
				dff.CalculateArea(NULL, _newWidth, NONE);
			}
		} else {
			cout << "[ShiftAdd] Error: No width assigned for the shift-and-add circuit" << endl;
			exit(-1);
		}
		
		// Assume the INV and NAND2 are on the same row and the total width of them is smaller than the adder or DFF
		if (spikingMode == NONSPIKING) {	// NONSPIKING: binary format
			height = adder.height + tech.featureSize*MAX_TRANSISTOR_HEIGHT /* INV and NAND2 */ + dff.height;
			width = _newWidth;
		} else {	// SPIKING: count spikes
			height = tech.featureSize*MAX_TRANSISTOR_HEIGHT /* INV and NAND2 */ + dff.height;
			width = _newWidth;
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

	}
}

void ShiftAdd::CalculateLatency(double numRead) {
	if (!initialized) {
		cout << "[ShiftAdd] Error: Require initialization first!" << endl;
	} else {
		readLatency = 0;
		
		// Assume the delay of INV and NAND2 are negligible
		if (spikingMode == NONSPIKING) {   // NONSPIKING: binary format
			// We can shift and add the weighted sum data in the next vector pulse integration cycle
			// Thus the shift-and-add time can be partially hidden by the vector pulse integration time at the next cycle
			// But there is at least one time of shift-and-add, which is at the last vector pulse cycle
			adder.CalculateLatency(1e20, dff.capTgDrain, 1);
			dff.CalculateLatency(1e20, 1);
			double shiftAddLatency = adder.readLatency + dff.readLatency;
			if (shiftAddLatency > cell.readPulseWidth)    // Completely hidden in the vector pulse cycle if smaller
				readLatency += (shiftAddLatency - cell.readPulseWidth) * (numRead - 1);
			readLatency += shiftAddLatency;    // At least need one time of shift-and-add
		} else {	// SPIKING: count spikes
			// We can shift out the weighted sum data in the next vector pulse integration cycle
			// Thus the shiftout time can be partially hidden by the vector pulse integration time at the next cycle
			// But there is at least one time of shiftout, which is at the last vector pulse cycle
			dff.CalculateLatency(1e20, numBitPerDff);	// Need numBitPerDff cycles to shift out the weighted sum data
			double shiftLatency = dff.readLatency;
			if (shiftLatency > cell.readPulseWidth)	// Completely hidden in the vector pulse cycle if smaller
				readLatency += (shiftLatency - cell.readPulseWidth) * (numRead - 1);
			readLatency += shiftLatency;	// At least need one time of shiftout
		}
	}
}

void ShiftAdd::CalculatePower(double numRead) {
	if (!initialized) {
		cout << "[ShiftAdd] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;
		if (spikingMode == NONSPIKING) {	// NONSPIKING: binary format
			adder.CalculatePower(numRead, numAdder);
			dff.CalculatePower(numRead, numDff);
			readDynamicEnergy += adder.readDynamicEnergy;
			readDynamicEnergy += dff.readDynamicEnergy;
			leakage += adder.leakage;
			leakage += dff.leakage;
		} else {	// SPIKING: count spikes
			dff.CalculatePower(numRead, numDff);
			readDynamicEnergy += dff.readDynamicEnergy;
			leakage += dff.leakage;
		}

		if (!readLatency) {
			//cout << "[ShiftAdd] Error: Need to calculate read latency first" << endl;
		} else {
			readPower = readDynamicEnergy/readLatency;
		}
	}
}

void ShiftAdd::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

