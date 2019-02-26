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

#ifndef SUBARRAY_H_
#define SUBARRAY_H_

#include "typedef.h"
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"
#include "formula.h"
#include "FunctionUnit.h"
#include "Adder.h"
#include "RowDecoder.h"
#include "Mux.h"
#include "WLDecoderOutput.h"
#include "DFF.h"
#include "VoltageSenseAmp.h"
#include "Precharger.h"
#include "SenseAmp.h"
#include "DecoderDriver.h"
#include "SRAMWriteDriver.h"
#include "ReadCircuit.h"
#include "SwitchMatrix.h"
#include "ShiftAdd.h"
#include "Subtractor.h"

// include the .h file for the new circuit component
#include "WLNewDecoderDriver.h"
#include "NewSwitchMatrix.h"

class SubArray: public FunctionUnit {
public:
	SubArray(InputParameter& _inputParameter, Technology& _tech, MemCell& _cell);
	virtual ~SubArray() {}
	InputParameter& inputParameter;
	Technology& tech;
	MemCell& cell;

	/* Functions */
	void PrintProperty();
	void Initialize(int _numRow, int _numCol, double _unitWireRes);
	void CalculateArea();
	void CalculateLatency(double _rampInput);
	void CalculatePower();

	/* Properties */
	bool initialized;	// Initialization flag
	int numRow;			// Number of rows
	int numCol;			// Number of columns
	int numColMuxed;	// How many columns share 1 read circuit (for neuro mode with analog RRAM) or 1 S/A (for memory mode or neuro mode with digital RRAM)
	int numWriteColMuxed;	// // How many columns share 1 write column decoder driver (for memory or neuro mode with digital RRAM)
	double heightArray;
	double widthArray;
	double areaArray;
	double readDynamicEnergyArray, writeDynamicEnergyArray;

	double lengthRow;	// Length of rows, Unit: m
	double lengthCol;	// Length of columns, Unit: m
	double unitWireRes;	// Unit wire resistance, Unit ohm/m
	double capRow1;		// Capacitance of row1 (BL for 1T1R, WL for Cross-point), Unit: F
	double capRow2;		// Capacitance of row2 (WL for 1T1R), Unit: F
	double capCol;		// Capacitance of column, Unit: F
	double resRow;		// Row resistance, Unit: ohm
	double resCol;		// Column resistance, Unit: ohm
	double resCellAccess;	// Resistance of access device, Unit: ohm
	double capCellAccess;	// Capacitance of access device, Unit: ohm
	double colDelay;	// Column delay, Unit: s
    double widthAccessTransistor; // for the 3T1C cell only

	double activityRowWrite;	// Activity for # of rows in the write
	double activityColWrite;	// Activity for # of columns in the write
	double activityRowRead;		// Activity for # of rows in the read
	int numReadPulse;		// # of read pulses for the input vector
	double numWritePulse;	// Average number of write pulse
	int maxNumWritePulse;	// Max # of write pulses for the device
	int maxNumIntBit;		// Max # bits for the integrate-and-fire neuron
	bool digitalModeNeuro;	// Digital RRAM in neuro mode
	ReadCircuitMode readCircuitMode;
	int numWriteCellPerOperationNeuro;  // For multifunctional mode or neuro only mode
	double clkFreq;			// Clock frequency (Hz)
	int numCellPerSynapse;  // Number of cells per synapse
	int avgWeightBit;		// Average weight for each synapse (value can range from 0 to numCellPerSynapse)
	int numReadCellPerOperationNeuro;   // Parameter for SRAM (use numBitPerCell cells to represent one D)
	SpikingMode spikingMode;	// NONSPIKING: input data using pulses in binary representation
								// SPIKING: input data using # of pulses
	bool shiftAddEnable;    // 0 for partition because the shift-and-add circuit will be after the last A&R stage

	bool relaxArrayCellHeight;	// true: relax the memory cell height to match the height of periperal circuit unit that connects to the row (ex: standard cell height in the last stage of row decoder) if the latter is larger
	bool relaxArrayCellWidth;	// true: relax the memory cell width to match the width of periperal circuit unit that connects to the column (ex: pass gate width in the column mux) if the latter is larger
    
    // parallel readout for STT-MRAM
    // for digital neuromode only
    bool parallelRead;
    
	/* Circuit modules */
	RowDecoder  wlDecoder;
	WLDecoderOutput wlDecoderOutput;
	Mux     mux, mux_2to1;
	RowDecoder  muxDecoder, muxDecoder_2to1;
	SwitchMatrix    slSwitchMatrix; 
	SwitchMatrix    blSwitchMatrix;     // In hybridcell case, it is for MSB BL
	SwitchMatrix    wlSwitchMatrix;     // In hybridcell case, it is for MSB access transistor gate
    SwitchMatrix    plSwitchMatrix;     //  controls the two power gate for 3T1C+2PCM and 2T1F cells;
    SwitchMatrix    wlSwitchMatrix_LSB; // for the LSB access transistor gate 
    SwitchMatrix    blSwitchMatrix_LSB; // for the LSB access transistor source
    
    // add the new switch matrix for BL and WL
    NewSwitchMatrix wlBlSwitchMatrix;
    // add the new WL decoder driver
    WLNewDecoderDriver wlDecoderDriverNew;
    
	ReadCircuit readCircuit;
	VoltageSenseAmp voltageSenseAmp;
	Precharger precharger;
	SenseAmp senseAmp;
	RowDecoder colDecoder;
	DecoderDriver wlDecoderDriver;
	DecoderDriver colDecoderDriver;
	SRAMWriteDriver sramWriteDriver;
	DFF dff;
	Adder adder;
	ShiftAdd shiftAdd;
	Subtractor subtractor;
};

#endif /* SUBARRAY_H_ */
