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
#include "Subtractor.h"

using namespace std;

Subtractor::Subtractor(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit() {
	initialized = false;
}

void Subtractor::Initialize(int _numBit, int _numSubtractor){
	if (initialized)
		cout << "[Subtractor] Warning: Already initialized!" << endl;
	
	numBit = _numBit;
	numSubtractor = _numSubtractor;
	
	// NAND
	widthNandN = 2 * MIN_NMOS_SIZE * tech.featureSize;
	widthNandP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	
	// INV
	widthInvN = MIN_NMOS_SIZE * tech.featureSize;
	widthInvP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	
	// NOR
	widthNorN = MIN_NMOS_SIZE * tech.featureSize;
	widthNorP = 2 * tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;

	initialized = true;
}

void Subtractor::CalculateArea(double _newHeight, double _newWidth, AreaModify _option) {
	if (!initialized) {
		cout << "[Subtractor] Error: Require initialization first!" << endl;
	} else {
		double hNand, wNand, hInv, wInv, hNor, wNor;
		// NAND2
		CalculateGateArea(NAND, 2, widthNandN, widthNandP, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hNand, &wNand);
		// INV
		CalculateGateArea(INV, 1, widthInvN, widthInvP, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hInv, &wInv);
		// NOR
		CalculateGateArea(NOR, 2, widthNorN, widthNorP, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hNor, &wNor);
		
		if (_newHeight && _option==NONE) {   // Subtractor in multiple columns given the total height
			hSubtractor = MAX(hNand, hNor);
			wSubtractor = (wNand*10 + wNor + wInv) * numBit;   // a full subtactor contains 10* NAND, 1* NOR and 1* INV

			// Calculate the number of Subtractor per column
			int numSubtractorPerCol = (int)(_newHeight/hSubtractor);
			if (numSubtractorPerCol > numSubtractor) {
				numSubtractorPerCol = numSubtractor;
			}
			int numColSubtractor = (int)ceil((double)numSubtractor / numSubtractorPerCol);
			height = _newHeight;
			width = wSubtractor * numColSubtractor;
			
		} else if (_newWidth && _option==NONE) { // Subtractor in multiple rows given the total width
			hSubtractor = MAX(hNand, hNor) * numBit;
			wSubtractor = wNand*10 + wNor + wInv;

			// Calculate the number of Subtractor per row
			int numSubtractorPerRow = (int)(_newWidth/wSubtractor);
			if (numSubtractorPerRow > numSubtractor) {
				numSubtractorPerRow = numSubtractor;
			}
			int numRowSubtractor = (int)ceil((double)numSubtractor / numSubtractorPerRow);
			width = _newWidth;
			height = hSubtractor * numRowSubtractor;

		} else {    // Assume one row of Subtractor by default
			hSubtractor = MAX(hNand, hNor);
			wSubtractor = (wNand*10 + wNor + wInv) * numBit;
			width = wSubtractor * numSubtractor;
			height = hSubtractor;
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
		// NOR2 capacitance
		CalculateGateCapacitance(NOR, 2, widthNorN, widthNorP, hNor, tech, &capNorInput, &capNorOutput);
		// INV capacitance
		CalculateGateCapacitance(INV, 1, widthInvN, widthInvP, hInv, tech, &capInvInput, &capInvOutput);
	}
}

void Subtractor::CalculateLatency(double _rampInput, double _capLoad, double numRead){
	if (!initialized) {
		cout << "[Subtractor] Error: Require initialization first!" << endl;
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

		// Calibration data pattern is X=000000..., Y=111111... and Bin=1
		// 1st
		resPullUp = CalculateOnResistance(widthNandP, PMOS, inputParameter.temperature, tech);
		tr = resPullUp * (capNandOutput + capNandInput*2);
		gm = CalculateTransconductance(widthNandP, PMOS, tech);
		beta = 1 / (resPullUp * gm);
		readLatency += horowitz(tr, beta, ramp[0], &ramp[1]);
        //printf("%f, %f, %f, %f,%.4e, %f\n",resPullUp,tr,gm,beta,readLatency,&ramp[1]);
		
		// 2nd
		resPullDown = CalculateOnResistance(widthNandN, NMOS, inputParameter.temperature, tech)*2;
		tr = resPullDown * (capNandOutput + capNandInput*2);
		gm = CalculateTransconductance(widthNandN, NMOS, tech);
		beta = 1 / (resPullDown * gm);
		readLatency += horowitz(tr, beta, ramp[1], &ramp[2]);
		
		// 3rd
		resPullUp = CalculateOnResistance(widthNandP, PMOS, inputParameter.temperature, tech);
		tr = resPullUp * (capNandOutput + capNandInput*2);
		gm = CalculateTransconductance(widthNandP, PMOS, tech);
		beta = 1 / (resPullUp * gm);
		readLatency += horowitz(tr, beta, ramp[2], &ramp[3]);

		// 4th
		resPullDown = CalculateOnResistance(widthNandN, NMOS, inputParameter.temperature, tech)*2;
		tr = resPullDown * (capNandOutput + capNandInput*2);
		gm = CalculateTransconductance(widthNandN, NMOS, tech);
		beta = 1 / (resPullDown * gm);
		readLatency += horowitz(tr, beta, ramp[3], &ramp[4]);
		
		// 5th
		resPullUp = CalculateOnResistance(widthNandP, PMOS, inputParameter.temperature, tech);
		tr = resPullUp * (capNandOutput + capNandInput*2);
		gm = CalculateTransconductance(widthNandP, PMOS, tech);
		beta = 1 / (resPullUp * gm);
		readLatency += horowitz(tr, beta, ramp[4], &ramp[5]);

		// 6th
		resPullDown = CalculateOnResistance(widthNandN, NMOS, inputParameter.temperature, tech)*2;
		tr = resPullDown * (capNandOutput + capNorInput);
		gm = CalculateTransconductance(widthNandN, NMOS, tech);
		beta = 1 / (resPullDown * gm);
		readLatency += horowitz(tr, beta, ramp[5], &ramp[6]);
		
		// 7th
		resPullDown = CalculateOnResistance(widthNorN, NMOS, inputParameter.temperature, tech);
		tr = resPullDown * (capNorOutput + capInvInput);
		gm = CalculateTransconductance(widthNorN, NMOS, tech);
		beta = 1 / (resPullDown * gm);
		readLatency += horowitz(tr, beta, ramp[6], &ramp[7]);
		
		// 8th
		resPullUp = CalculateOnResistance(widthInvP, PMOS, inputParameter.temperature, tech);
		tr = resPullUp * (capInvOutput + capLoad);
		gm = CalculateTransconductance(widthInvP, PMOS, tech);
		beta = 1 / (resPullUp * gm);
		readLatency += horowitz(tr, beta, ramp[7], &ramp[8]);

        //printf("readLatency is %.4e\n", readLatency);		
		readLatency *= numBit;
		readLatency *= numRead;
		rampOutput = ramp[8];
        //printf("number of bit is %.4e\n", numBit);
        //printf("number of read is %.4e\n", numRead);       
	}
}

void Subtractor::CalculatePower(double numRead, int numSubtractorPerOperation) {
	if (!initialized) {
		cout << "[Subtractor] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;
		
		/* Leakage power */
		leakage += CalculateGateLeakage(NAND, 2, widthNandN, widthNandP, inputParameter.temperature, tech) * tech.vdd * 10 * numBit * numSubtractor;
		leakage += CalculateGateLeakage(NOR, 2, widthNorN, widthNorP, inputParameter.temperature, tech) * tech.vdd * numBit * numSubtractor;
		leakage += CalculateGateLeakage(INV, 1, widthInvN, widthInvP, inputParameter.temperature, tech) * tech.vdd * numBit * numSubtractor;

		/* Read Dynamic energy */
		// Calibration data pattern of critical path is X=000000..., Y=111111... and Bin=1
		// Only count 0 to 1 transition for energy

		// I/O transition
		readDynamicEnergy += (capNandInput * 4) * tech.vdd * tech.vdd * numBit;     // Input of Y[n] and Bin
		readDynamicEnergy += (capInvOutput * 1) * tech.vdd * tech.vdd * numBit;     // Bout
		
		// Hidden transition
		// First stage
		readDynamicEnergy += (capNandOutput + capNandInput*2) * tech.vdd * tech.vdd;	
		readDynamicEnergy += (capNandOutput + capNandInput) * tech.vdd * tech.vdd;	    
		readDynamicEnergy += (capNandOutput + capNandInput*2) * tech.vdd * tech.vdd;	
		readDynamicEnergy += (capNandOutput + capNandInput) * tech.vdd * tech.vdd;		
		readDynamicEnergy += (capNandOutput + capNandInput*3) * tech.vdd * tech.vdd;	
		readDynamicEnergy += (capNandOutput + capNorInput) * tech.vdd * tech.vdd;	

		// Rest stages	
		readDynamicEnergy *= (numBit-1);	
		
		readDynamicEnergy *= MIN(numSubtractorPerOperation, numSubtractor) * numRead;

	}
}

void Subtractor::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}


