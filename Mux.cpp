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
#include "Mux.h"

using namespace std;

Mux::Mux(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit() {
	initialized = false;
}

void Mux::Initialize(int _numInput, int _numSelection, double _resTg, bool _digital){
	if (initialized)
		cout << "[Mux] Warning: Already initialized!" << endl;

	numInput = _numInput;
	numSelection = _numSelection;	/* Typically numColMuxed */
	digital = _digital;

	// Mux
	if (digital) {	// Assume digital Mux has standard NMOS and PMOS width
		widthTgN = MIN_NMOS_SIZE * tech.featureSize;
		widthTgP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
		resTg = 1 / (1/CalculateOnResistance(widthTgN, NMOS, inputParameter.temperature, tech) 
					+ 1/CalculateOnResistance(widthTgP, PMOS, inputParameter.temperature, tech));
	} else {
		resTg = _resTg;
		widthTgN = CalculateOnResistance(tech.featureSize, NMOS, inputParameter.temperature, tech)
								* tech.featureSize / (resTg*2);
		widthTgP = CalculateOnResistance(tech.featureSize, PMOS, inputParameter.temperature, tech)
								* tech.featureSize / (resTg*2);
	}
	initialized = true;
}

void Mux::CalculateArea(double _newHeight, double _newWidth, AreaModify _option){
	if (!initialized) {
		cout << "[Mux] Error: Require initialization first!" << endl;
	} else {
		double hTg, wTg;
		int numTg = numInput * numSelection;

		// TG
		if (digital) {	// Digital Mux
			CalculateGateArea(INV, 1, widthTgN, widthTgP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hTg, &wTg);
			if (_newWidth && _option==NONE) { // Tg in multiple rows given the total width
				// Calculate the number of Tg per row
				int numTgPerRow = (int)(_newWidth/wTg);
				if (numTgPerRow > numTg) {
					numTgPerRow = numTg;
				}
				numRowTg = (int)ceil((double)numTg / numTgPerRow);
				width = _newWidth;
				height = hTg * numRowTg;
			} else {    // Assume one row of Tg by default
				width = wTg * numTg;
				height = hTg;
			}
		} else {	// Analog Mux
			if (_newWidth && _option==NONE) {
				numRowTg = 1;
				double minCellWidth = 2 * (POLY_WIDTH + MIN_GAP_BET_GATE_POLY) * tech.featureSize; // min standard cell width
				if (minCellWidth > _newWidth) {
					cout << "[Mux] Error: pass gate width is even larger than the array width" << endl;
				}

				int numTgPerRow = (int)(_newWidth / minCellWidth);	// Get max # Tg per row (this is not the final # Tg per row because the last row may have less # Tg)
				numRowTg = (int)ceil((double)numTg / numTgPerRow);	// Get min # rows based on this max # Tg per row
				numTgPerRow = (int)ceil((double)numTg / numRowTg);	// Get # Tg per row based on this min # rows
				TgWidth = _newWidth / numTgPerRow;
				int numFold = (int)(TgWidth / (0.5*minCellWidth)) - 1;  // Get the max number of folding

				// widthTgN, widthTgP and numFold can determine the height and width of each pass gate
				CalculatePassGateArea(widthTgN, widthTgP, tech, numFold, &hTg, &wTg);

				width = _newWidth;
				height = hTg * numRowTg;

			} else {
				// Default (just use pass gate without folding)
				CalculatePassGateArea(widthTgN, widthTgP, tech, 1, &hTg, &wTg);
				height = hTg;
				width = wTg * numTg;
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

		widthTgShared = width/numInput;

		capTgGateN = CalculateGateCap(widthTgN, tech);
		capTgGateP = CalculateGateCap(widthTgP, tech);
		CalculateGateCapacitance(INV, 1, widthTgN, widthTgP, hTg, tech, NULL, &capTgDrain);
	}
}

void Mux::CalculateLatency(double _rampInput, double _capLoad, double numRead) {  // rampInput is from SL/BL, not fron EN signal
	if (!initialized) {
		cout << "[Mux] Error: Require initialization first!" << endl;
	} else {
		rampInput = _rampInput;
		capLoad = _capLoad;
		double resPullDown;
		double tr;  	/* time constant */
		double gm;  	/* transconductance */
		double beta;    /* for horowitz calculation */
		double rampNandOutput;
		readLatency = 0;

		// TG
		tr = resTg*2 * (capTgDrain + 0.5*capTgGateN + 0.5*capTgGateP + capLoad);	// Calibration: use resTg*2 (only one transistor is transmitting signal in the pass gate) may be more accurate, and include gate cap because the voltage at the source of NMOS and drain of PMOS is changing (assuming Cg = 0.5Cgs + 0.5Cgd)
		readLatency += 2.3 * tr;	// 2.3 means charging from 0% to 90%

		readLatency *= numRead;
	}
}

void Mux::CalculatePower(double numRead) {
	if (!initialized) {
		cout << "[Mux] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;

		// TG gates only
		readDynamicEnergy += capTgGateN * numInput * tech.vdd * tech.vdd;	// Selected pass gates (OFF to ON)
		//readDynamicEnergy += (capTgDrain * 2) * numInput * cell.readVoltage * cell.readVoltage;	// Selected pass gates (OFF to ON)
		readDynamicEnergy += capTgGateP * numInput * tech.vdd * tech.vdd;	// Deselected pass gates (ON to OFF)
		
		readDynamicEnergy *= numRead;
		if (!readLatency) {
			//cout << "[Mux] Error: Need to calculate read latency first" << endl;
		} else {
			readPower = readDynamicEnergy/readLatency;
		}

	}
}

void Mux::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

