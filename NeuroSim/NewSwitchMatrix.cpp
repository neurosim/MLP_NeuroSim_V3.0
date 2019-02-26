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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// This NewSwitchMatrix is used for new BNN Parallel RRAM mode, only for WL and BL... Row connected //////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


#include <cmath>
#include <iostream>
#include "constant.h"
#include "formula.h"
#include "NewSwitchMatrix.h"

using namespace std;

NewSwitchMatrix::NewSwitchMatrix(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), dff(_inputParameter, _tech, _cell), FunctionUnit() {
	// TODO Auto-generated constructor stub
	initialized = false;
}

NewSwitchMatrix::~NewSwitchMatrix() {
	// TODO Auto-generated destructor stub
}

void NewSwitchMatrix::Initialize(int _numOutput, double _activityRowRead, double _clkFreq, bool _XNOR){
	if (initialized)
		cout << "[NewSwitchMatrix] Warning: Already initialized!" << endl;
	
	numOutput = _numOutput;
	activityRowRead = _activityRowRead;
	clkFreq = _clkFreq;
	XNOR = _XNOR;
    
	// DFF
	dff.Initialize(numOutput, clkFreq);     // used for scan in ....
    
	widthTgN = MIN_NMOS_SIZE * tech.featureSize;
	widthTgP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	
	
	initialized = true;
}

void NewSwitchMatrix::CalculateArea(double _newHeight, double _newWidth, AreaModify _option) {
	if (!initialized) {
		cout << "[NewSwitchMatrix] Error: Require initialization first!" << endl;
	} else {
		double hTg, wTg;
		
		double minCellHeight = MAX_TRANSISTOR_HEIGHT * tech.featureSize;   // min standard cell height for 1 Tg 
		if (_newHeight && _option==NONE) {
			if (_newHeight < minCellHeight) {
				cout << "[NewSwitchMatrix] Error: pass gate height is even larger than the array height" << endl;
			}

		    int numTgPerCol = (int)(_newHeight / minCellHeight);	// Get max # Tg per column (this is not the final # Tg per column because the last column may have less # Tg)
		    numColTg = (int)ceil((double)numOutput / numTgPerCol);	// Get min # columns based on this max # Tg per column
		    numTgPerCol = (int)ceil((double)numOutput / numColTg);		// Get # Tg per column based on this min # columns
		    TgHeight = _newHeight / numTgPerCol;        // release TG height
		    CalculateGateArea(INV, 1, widthTgN, widthTgP, TgHeight, tech, &hTg, &wTg);         // calculate released TG layout height and width
		
		    // DFF
		    dff.CalculateArea(_newHeight, NULL, NONE);
		    //dff.CalculateArea(_newHeight, NULL, MAGIC);

				
		    height = _newHeight;
		    width = (wTg * 4) * numColTg + dff.width;      // add switch matrix and dff together

		} else {       // MAGIC or OVERRIDE ...
			CalculateGateArea(INV, 1, widthTgN, widthTgP, minCellHeight, tech, &hTg, &wTg); // Pass gate with folding
			height = hTg * numOutput;
			dff.CalculateArea(height, NULL, NONE);	// Need to give the height information, otherwise by default the area calculation of DFF is in column mode
			width = (wTg * 4) + dff.width;
		}

		
		if (XNOR) {
			area = height*width/2;
		} else {
			area = height * width;
		}

	    // Modify layout
	    newHeight = _newHeight;
	    newWidth = _newWidth;
	    switch (_option) {
		    case MAGIC:
			    MagicLayout();       // if MAGIC, call Magiclayout() in FunctionUnit.cpp
			    break;
		    case OVERRIDE:
			    OverrideLayout();    // if OVERRIDE, call Overridelayout() in FunctionUnit.cpp
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


void NewSwitchMatrix::CalculateLatency(double _rampInput, double _capLoad, double _resLoad, double numRead, double numWrite) {	// For simplicity, assume shift register is ideal
	if (!initialized) {
		cout << "[NewSwitchMatrix] Error: Require initialization first!" << endl;
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
		capOutput = capTgDrain * 5;         // pass 2 TG, 5 loading drain capacitance
		tr = resTg * (capOutput + capLoad) + resLoad * capLoad / 2;     // elmore delay model
		readLatency += horowitz(tr, 0, rampInput, &rampOutput);	// get from chargeLatency in the original SubArray.cpp
		
		readLatency *= numRead;
		readLatency += dff.readLatency;

		writeLatency = cell.writePulseWidth;     // write latency determined by write pulse width
		writeLatency *= numWrite;
		writeLatency += dff.readLatency;	// Use DFF read latency here because no write in the DFF module
	}
}

void NewSwitchMatrix::CalculatePower(double numRead, double numWrite) {      
	if (!initialized) {
		cout << "[NewSwitchMatrix] Error: Require initialization first!" << endl;
	} else {
		
		leakage = 0;
		readDynamicEnergy = 0;
		writeDynamicEnergy = 0;
		
		// DFF
		dff.CalculatePower(numRead, numOutput);	// Use numOutput since every DFF will pass signal (either 0 or 1)

		// Leakage power
		leakage += dff.leakage;	// Only DFF has leakage, assuming TG do not have leakage

		// Read dynamic energy
		readDynamicEnergy += (capTgDrain * 2) * cell.accessVoltage * cell.accessVoltage * numOutput * activityRowRead;   // 1 TG pass Vaccess to CMOS gate to select the row
		readDynamicEnergy += (capTgDrain * 5) * cell.readVoltage * cell.readVoltage * numOutput * activityRowRead;    // 2 TG pass Vread to BL, total loading is 5 Tg Drain capacitance
		readDynamicEnergy += (capTgGateN + capTgGateP) * 3 * tech.vdd * tech.vdd * numOutput * activityRowRead;    // open 3 TG when selected

		readDynamicEnergy *= numRead;
		readDynamicEnergy += dff.readDynamicEnergy;
		
		if (!readLatency) {
			//cout << "[Switch Matrix] Error: Need to calculate read latency first" << endl;
		} else {
			readPower = readDynamicEnergy/readLatency;
		}
		
		// Write dynamic energy (2-step write and average case half SET and half RESET)
		// 1T1R
		// connect to rows, when writing, pass GND to BL, no transmission energy acrossing BL
		writeDynamicEnergy += (capTgDrain * 2) * cell.accessVoltage * cell.accessVoltage;    // 1 TG pass Vaccess to CMOS gate to select the row
		writeDynamicEnergy += (capTgGateN + capTgGateP) * 2 * tech.vdd * tech.vdd * 2;    // open 2 TG when Q selected, and *2 means switching from one selected row to another
        writeDynamicEnergy += (capTgGateN + capTgGateP) * tech.vdd * tech.vdd;    // always open one TG when writing		
	}
	writeDynamicEnergy *= numWrite;
	writeDynamicEnergy += dff.readDynamicEnergy;	// Use DFF read energy here because no write in the DFF module
	if (!writeLatency) {
		//cout << "[Switch Matrix] Error: Need to calculate write latency first" << endl;
	} else {
		writePower = writeDynamicEnergy/writeLatency;
	}
}


void NewSwitchMatrix::PrintProperty(const char* str) {
	//cout << "NewSwitchMatrix Properties:" << endl;
	FunctionUnit::PrintProperty(str);
}

