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
#include "typedef.h"
#include "formula.h"
#include "Adder.h"

using namespace std;

Adder::Adder(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit() {
	initialized = false;
}

void Adder::Initialize(int _numBit, int _numAdder){
	if (initialized)
		cout << "[Adder] Warning: Already initialized!" << endl;
	
	numBit = _numBit;
	numAdder = _numAdder;
	
	widthNandN = 2 * MIN_NMOS_SIZE * tech.featureSize;
	widthNandP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;

	initialized = true;
}

void Adder::CalculateArea(double _newHeight, double _newWidth, AreaModify _option) {
	if (!initialized) {
		cout << "[Adder] Error: Require initialization first!" << endl;
	} else {
		double hNand, wNand;
		// NAND2
		CalculateGateArea(NAND, 2, widthNandN, widthNandP, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hNand, &wNand);
		
		if (_newHeight && _option==NONE) {   // Adder in multiple columns given the total height
			hAdder = hNand;
			wAdder = wNand * 9 * numBit;

			// Calculate the number of adder per column
			int numAdderPerCol = (int)(_newHeight/hAdder);
			if (numAdderPerCol > numAdder) {
				numAdderPerCol = numAdder;
			}
			int numColAdder = (int)ceil((double)numAdder / numAdderPerCol);
			height = _newHeight;
			width = wAdder * numColAdder;
			
		} else if (_newWidth && _option==NONE) { // Adder in multiple rows given the total width
			hAdder = hNand * numBit;
			wAdder = wNand * 9;

			// Calculate the number of adder per row
			int numAdderPerRow = (int)(_newWidth/wAdder);
			if (numAdderPerRow > numAdder) {
				numAdderPerRow = numAdder;
			}
			int numRowAdder = (int)ceil((double)numAdder / numAdderPerRow);
			width = _newWidth;
			height = hAdder * numRowAdder;

		} else {    // Assume one row of adder by default
			hAdder = hNand;
			wAdder = wNand * 9 * numBit;
			width = wAdder * numAdder;
			height = hAdder;
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
		
		// NAND2 capacitance
		CalculateGateCapacitance(NAND, 2, widthNandN, widthNandP, hNand, tech, &capNandInput, &capNandOutput);
	}
}

void Adder::CalculateLatency(double _rampInput, double _capLoad, double numRead){
	if (!initialized) {
		cout << "[Adder] Error: Require initialization first!" << endl;
	} else {
		readLatency = 0;
		rampInput = _rampInput;
		capLoad = _capLoad;
		double tr;		/* time constant */
		double gm;		/* transconductance */
		double beta;	/* for horowitz calculation */
		double resPullUp, resPullDown;
		double readLatencyIntermediate = 0;
		double ramp[10];
		
		ramp[0] = rampInput;

		// Calibration data pattern is A=1111111..., B=1000000... and Cin=1
		// 1st
		resPullDown = CalculateOnResistance(widthNandN, NMOS, inputParameter.temperature, tech) * 2;
		tr = resPullDown * (capNandOutput + capNandInput * 3);
		gm = CalculateTransconductance(widthNandN, NMOS, tech);
		beta = 1 / (resPullDown * gm);
		readLatency += horowitz(tr, beta, ramp[0], &ramp[1]);
		
		// 2nd
		resPullUp = CalculateOnResistance(widthNandP, PMOS, inputParameter.temperature, tech);
		tr = resPullUp * (capNandOutput + capNandInput * 2);
		gm = CalculateTransconductance(widthNandP, PMOS, tech);
		beta = 1 / (resPullUp * gm);
		readLatency += horowitz(tr, beta, ramp[1], &ramp[2]);
		
		// 3rd
		resPullDown = CalculateOnResistance(widthNandN, NMOS, inputParameter.temperature, tech) * 2;
		tr = resPullDown * (capNandOutput + capNandInput * 3);
		gm = CalculateTransconductance(widthNandN, NMOS, tech);
		beta = 1 / (resPullDown * gm);
		readLatencyIntermediate += horowitz(tr, beta, ramp[2], &ramp[3]);

		// 4th
		resPullUp = CalculateOnResistance(widthNandP, PMOS, inputParameter.temperature, tech);
		tr = resPullUp * (capNandOutput + capNandInput * 2);
		gm = CalculateTransconductance(widthNandP, PMOS, tech);
		beta = 1 / (resPullUp * gm);
		readLatencyIntermediate += horowitz(tr, beta, ramp[3], &ramp[4]);
		
		if (numBit > 2) {
			readLatency += readLatencyIntermediate * (numBit - 2);
		}
		
		// 5th
		resPullDown = CalculateOnResistance(widthNandN, NMOS, inputParameter.temperature, tech) * 2;
		tr = resPullDown * (capNandOutput + capNandInput * 3);
		gm = CalculateTransconductance(widthNandN, NMOS, tech);
		beta = 1 / (resPullDown * gm);
		readLatency += horowitz(tr, beta, ramp[4], &ramp[5]);

		// 6th
		resPullUp = CalculateOnResistance(widthNandP, PMOS, inputParameter.temperature, tech);
		tr = resPullUp * (capNandOutput + capNandInput);
		gm = CalculateTransconductance(widthNandP, PMOS, tech);
		beta = 1 / (resPullUp * gm);
		readLatency += horowitz(tr, beta, ramp[5], &ramp[6]);
		
		// 7th
		resPullDown = CalculateOnResistance(widthNandN, NMOS, inputParameter.temperature, tech) * 2;
		tr = resPullDown * (capNandOutput + capLoad);
		gm = CalculateTransconductance(widthNandN, NMOS, tech);
		beta = 1 / (resPullDown * gm);
		readLatency += horowitz(tr, beta, ramp[6], &ramp[7]);

		readLatency *= numRead;
		rampOutput = ramp[7];
	}
}

void Adder::CalculatePower(double numRead, int numAdderPerOperation) {
	if (!initialized) {
		cout << "[Adder] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;
		
		/* Leakage power */
		leakage += CalculateGateLeakage(NAND, 2, widthNandN, widthNandP, inputParameter.temperature, tech) * tech.vdd * 9 * numBit * numAdder;

		/* Read Dynamic energy */
		// Calibration data pattern of critical path is A=1111111..., B=1000000... and Cin=1
		// Only count 0 to 1 transition for energy
		// First stage
		readDynamicEnergy += (capNandInput * 6) * tech.vdd * tech.vdd;    // Input of 1 and 2 and Cin
        readDynamicEnergy += (capNandOutput * 2) * tech.vdd * tech.vdd;  // Output of S[0] and 5
		// Second and later stages
		readDynamicEnergy += (capNandInput * 7) * tech.vdd * tech.vdd * (numBit-1);
		readDynamicEnergy += (capNandOutput * 3) * tech.vdd * tech.vdd * (numBit-1);
		
		// Hidden transition
		// First stage
		readDynamicEnergy += (capNandOutput + capNandInput) * tech.vdd * tech.vdd * 2;	// #2 and #3
		readDynamicEnergy += (capNandOutput + capNandInput * 2) * tech.vdd * tech.vdd;	// #4
		readDynamicEnergy += (capNandOutput + capNandInput * 3) * tech.vdd * tech.vdd;	// #5
		readDynamicEnergy += (capNandOutput + capNandInput) * tech.vdd * tech.vdd;		// #6
		// Second and later stages
		readDynamicEnergy += (capNandOutput + capNandInput * 3) * tech.vdd * tech.vdd * (numBit-1);	// # 1
		readDynamicEnergy += (capNandOutput + capNandInput) * tech.vdd * tech.vdd * (numBit-1);		// # 3
		readDynamicEnergy += (capNandOutput + capNandInput) * tech.vdd * tech.vdd * 2 * (numBit-1);		// #6 and #7
		
		readDynamicEnergy *= MIN(numAdderPerOperation, numAdder) * numRead;

	}
}

void Adder::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}


