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
#include "SenseAmp.h"

using namespace std;

SenseAmp::SenseAmp(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit() {
	initialized = false;
}

void SenseAmp::Initialize(int _numCol, bool _currentSense, double _senseVoltage, double _pitchSenseAmp, double _clkFreq, int _numReadCellPerOperationNeuro) {
	if (initialized)
		cout << "[SenseAmp] Warning: Already initialized!" << endl;

	numCol = _numCol;
	currentSense = _currentSense;
	senseVoltage = _senseVoltage;
	pitchSenseAmp = _pitchSenseAmp;
	clkFreq = _clkFreq;
	numReadCellPerOperationNeuro = _numReadCellPerOperationNeuro;

	if (pitchSenseAmp <= tech.featureSize * 6) {
		puts("[SenseAmp] Error: pitch too small, cannot do the layout");
	}

	initialized = true;
}

void SenseAmp::CalculateArea(double _newHeight, double _newWidth, AreaModify _option) {
	if (!initialized) {
		cout << "[SenseAmp] Error: Require initialization first!" << endl;
	} else {
		height = width = area = 0;
		if (currentSense) {
			// TODO
			puts("[SenseAmp] Do not support current sensing yet");
			exit(-1);
		} 
		double hSenseP, wSenseP, hSenseN, wSenseN, hSenseIso, wSenseIso, hSenseEn, wSenseEn;
		
		// Exchange width and height as in the original code
		CalculateGateArea(INV, 1, 0, W_SENSE_P * tech.featureSize, pitchSenseAmp, tech, &wSenseP, &hSenseP);
		CalculateGateArea(INV, 1, 0, W_SENSE_ISO * tech.featureSize, pitchSenseAmp, tech, &wSenseIso, &hSenseIso);
		CalculateGateArea(INV, 1, W_SENSE_N * tech.featureSize, 0, pitchSenseAmp, tech, &wSenseN, &hSenseN);
		CalculateGateArea(INV, 1, W_SENSE_EN * tech.featureSize, 0, pitchSenseAmp, tech, &wSenseEn, &hSenseEn);
		// Just sum up the area of all components
		area += (wSenseP * hSenseP) * 2 + (wSenseN * hSenseN) * 2 + wSenseIso * hSenseIso + wSenseEn * hSenseEn;
		area *= numCol;
		
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
		capLoad = CalculateGateCap((W_SENSE_P + W_SENSE_N) * tech.featureSize, tech)
				+ CalculateDrainCap(W_SENSE_N * tech.featureSize, NMOS, pitchSenseAmp, tech)
				+ CalculateDrainCap(W_SENSE_P * tech.featureSize, PMOS, pitchSenseAmp, tech)
				+ CalculateDrainCap(W_SENSE_ISO * tech.featureSize, PMOS, pitchSenseAmp, tech)
				+ CalculateDrainCap(W_SENSE_MUX * tech.featureSize, NMOS, pitchSenseAmp, tech);
	}
}

void SenseAmp::CalculateLatency(double numRead) {
	if (!initialized) {
		cout << "[SenseAmp] Error: Require initialization first!" << endl;
	} else {
		readLatency = 0;

		/* Voltage sense amplifier */
		double gm = CalculateTransconductance(W_SENSE_N * tech.featureSize, NMOS, tech)
				+ CalculateTransconductance(W_SENSE_P * tech.featureSize, PMOS, tech);
		double tau = capLoad / gm;
		readLatency += tau * log(tech.vdd / senseVoltage);
		readLatency += 1/clkFreq;   // Clock time for S/A enable

		readLatency *= numRead;
	}
}

void SenseAmp::CalculatePower(double numRead) {
	if (!initialized) {
		cout << "[SenseAmp] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;

		/* Voltage sense amplifier */
		// Leakage
		double idleCurrent =  CalculateGateLeakage(INV, 1, W_SENSE_EN * tech.featureSize, 0, inputParameter.temperature, tech) * tech.vdd;
		leakage += idleCurrent * tech.vdd * numCol;
		
		// Dynamic energy
		readDynamicEnergy += capLoad * tech.vdd * tech.vdd;
		readDynamicEnergy *= MIN(numReadCellPerOperationNeuro, numCol) * numRead;
	}
}

void SenseAmp::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

