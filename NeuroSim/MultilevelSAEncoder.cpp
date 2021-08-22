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
#include "MultilevelSAEncoder.h"

using namespace std;

MultilevelSAEncoder::MultilevelSAEncoder(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit() {
	initialized = false;
}

void MultilevelSAEncoder::Initialize(int _numLevel, int _numEncoder){
	if (initialized)
		cout << "[MultilevelSAEncoder] Warning: Already initialized!" << endl;
	
	numEncoder = _numEncoder;      // number of encoder needed
	numLevel= _numLevel;           // number of levels from MultilevelSA
	numInput = numLevel / 2;       // number of NAND gate in encoder
	numGate = log2(numLevel);      // number of NAND gate in encoder 
	
	widthInvN = MIN_NMOS_SIZE * tech.featureSize;
	widthInvP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	widthNandN = 2 * MIN_NMOS_SIZE * tech.featureSize;
	widthNandP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;

	initialized = true;
}

void MultilevelSAEncoder::CalculateArea(double _newHeight, double _newWidth, AreaModify _option) {
	if (!initialized) {
		cout << "[MultilevelSAEncoder] Error: Require initialization first!" << endl;
	} else {
        double wEncoder, hEncoder, wNand, hNand, wNandLg, hNandLg, wInv, hInv;
		
		// NAND2
		CalculateGateArea(NAND, 2, widthNandN, widthNandP, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hNand, &wNand);
		// INV
		CalculateGateArea(INV, 1, widthInvN, widthInvP, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hInv, &wInv);
		// Large NAND in Encoder
		CalculateGateArea(NAND, numInput, widthNandN, widthNandP, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hNandLg, &wNandLg);
		
		wEncoder = 2*wInv + wNand + wNandLg;
		hEncoder = max( (numLevel-1)*hInv, (numLevel-1)*hNand );
	    
		int numEncoderPerRow = (int)(_newWidth/wEncoder);
		if (numEncoderPerRow > numEncoder) {
			numEncoderPerRow = numEncoder;
		}
		int numRowEncoder = (int)ceil((double)numEncoder / numEncoderPerRow);
		width = _newWidth;
		height = hEncoder * numRowEncoder;
		
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
		// NAND2
		CalculateGateCapacitance(NAND, 2, widthNandN, widthNandP, hNand, tech, &capNandInput, &capNandOutput);
		// Large NAND in Encoder
		CalculateGateCapacitance(NAND, numInput, widthNandN, widthNandP, hNandLg, tech, &capNandLgInput, &capNandLgOutput);
	}
}

void MultilevelSAEncoder::CalculateLatency(double _rampInput, double numRead){
	if (!initialized) {
		cout << "[MultilevelSAEncoder] Error: Require initialization first!" << endl;
	} else {
		readLatency = 0;
		rampInput = _rampInput;
		double tr;		/* time constant */
		double gm;		/* transconductance */
		double beta;	/* for horowitz calculation */
		double resPullUp, resPullDown;
		double readLatencyIntermediate = 0;
		double ramp[10];
		
		ramp[0] = rampInput;

		// 1st INV to NAND2
		resPullDown = CalculateOnResistance(widthInvN, NMOS, inputParameter.temperature, tech) * 2;
		tr = resPullDown * (capInvOutput + capNandInput * 2);
		gm = CalculateTransconductance(widthNandN, NMOS, tech);
		beta = 1 / (resPullDown * gm);
		readLatency += horowitz(tr, beta, ramp[0], &ramp[1]);
		
		// 2nd NAND2 to Large NAND
		resPullUp = CalculateOnResistance(widthNandP, PMOS, inputParameter.temperature, tech);
		tr = resPullUp * (capNandOutput + capNandLgInput * numInput);
		gm = CalculateTransconductance(widthNandP, PMOS, tech);
		beta = 1 / (resPullUp * gm);
		readLatency += horowitz(tr, beta, ramp[1], &ramp[2]);
		
		// 3rd large NAND to INV
		resPullDown = CalculateOnResistance(widthNandN, NMOS, inputParameter.temperature, tech) * 2;
		tr = resPullDown * (capNandLgOutput + capInvInput);
		gm = CalculateTransconductance(widthNandN, NMOS, tech);
		beta = 1 / (resPullDown * gm);
		readLatencyIntermediate += horowitz(tr, beta, ramp[2], &ramp[3]);

		// 4th INV
		resPullUp = CalculateOnResistance(widthInvP, PMOS, inputParameter.temperature, tech);
		tr = resPullUp * capInvOutput;
		gm = CalculateTransconductance(widthNandP, PMOS, tech);
		beta = 1 / (resPullUp * gm);
		readLatencyIntermediate += horowitz(tr, beta, ramp[3], &ramp[4]);
		
		readLatency *= numRead;
		rampOutput = ramp[4];
	}
}

void MultilevelSAEncoder::CalculatePower(double numRead) {
	if (!initialized) {
		cout << "[MultilevelSAEncoder] Error: Require initialization first!" << endl;
	} else {
		readDynamicEnergy = 0;
		leakage = 0;

		leakage = CalculateGateLeakage(INV, 1, widthInvN, widthInvP, inputParameter.temperature, tech) * tech.vdd * (numLevel+numGate) * numEncoder
		          + CalculateGateLeakage(NAND, 2, widthNandN, widthNandP, inputParameter.temperature, tech) * tech.vdd * (numLevel+numGate) * numEncoder
				  + CalculateGateLeakage(NAND, numInput, widthNandN, widthNandP, inputParameter.temperature, tech) * tech.vdd * numGate * numEncoder;
		
		readDynamicEnergy += (capInvInput + capInvOutput) * tech.vdd * tech.vdd * (numLevel+numGate) * numEncoder;
		readDynamicEnergy += (capNandInput + capNandOutput) * tech.vdd * tech.vdd * (numLevel+numGate) * numEncoder;
		readDynamicEnergy += (capNandLgInput + capNandLgOutput) * tech.vdd * tech.vdd * numGate * numEncoder;
		readDynamicEnergy *= numRead;
		
		if (!readLatency) {
			//cout << "[MultilevelSenseAmp] Error: Need to calculate read latency first" << endl;
		} else {
			readPower = readDynamicEnergy/readLatency;
		}
	}
}

void MultilevelSAEncoder::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}


