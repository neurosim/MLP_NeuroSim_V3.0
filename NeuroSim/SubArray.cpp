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
#include "SubArray.h"

using namespace std;

SubArray::SubArray(InputParameter& _inputParameter, Technology& _tech, MemCell& _cell):
						inputParameter(_inputParameter), tech(_tech), cell(_cell), //initialize the class parameters in Subarray class
						wlDecoder(_inputParameter, _tech, _cell),
						wlDecoderOutput(_inputParameter, _tech, _cell),
						mux(_inputParameter, _tech, _cell),
                        mux_2to1(_inputParameter,_tech,_cell),
						muxDecoder(_inputParameter, _tech, _cell),
                        muxDecoder_2to1(_inputParameter, _tech, _cell),
						slSwitchMatrix(_inputParameter, _tech, _cell),
						blSwitchMatrix(_inputParameter, _tech, _cell),
						wlSwitchMatrix(_inputParameter, _tech, _cell),
						plSwitchMatrix(_inputParameter, _tech, _cell),
						blSwitchMatrix_LSB(_inputParameter, _tech, _cell),
						wlSwitchMatrix_LSB(_inputParameter, _tech, _cell),
                        // initialized the new circuit components
                        wlBlSwitchMatrix(_inputParameter, _tech, _cell),
                        wlDecoderDriverNew(_inputParameter, _tech, _cell),
						readCircuit(_inputParameter, _tech, _cell),
						voltageSenseAmp(_inputParameter, _tech, _cell),
						precharger(_inputParameter, _tech, _cell),
						senseAmp(_inputParameter, _tech, _cell),
						colDecoder(_inputParameter, _tech, _cell),
						wlDecoderDriver(_inputParameter, _tech, _cell),
						colDecoderDriver(_inputParameter, _tech, _cell),
						sramWriteDriver(_inputParameter, _tech, _cell),
						adder(_inputParameter, _tech, _cell),
						dff(_inputParameter, _tech, _cell),
						shiftAdd(_inputParameter, _tech, _cell),
						subtractor(_inputParameter, _tech, _cell){
	initialized = false;
	readDynamicEnergyArray = writeDynamicEnergyArray = 0;
}

void SubArray::Initialize(int _numRow, int _numCol, double _unitWireRes){  //initialization module
	if (initialized)
		cout << "[Subarray] Warning: Already initialized!" << endl;  //avioding initialize twice
	
	numRow = _numRow;    //import parameters
	numCol = _numCol;
	unitWireRes = _unitWireRes;
	
	double MIN_CELL_HEIGHT = MAX_TRANSISTOR_HEIGHT;  //set real layout cell height
	double MIN_CELL_WIDTH = (MIN_GAP_BET_GATE_POLY + POLY_WIDTH) * 2;  //set real layout cell width
	if (cell.memCellType == Type::SRAM) {  //if array is SRAM
		if (relaxArrayCellWidth) {  //if want to relax the cell width
			lengthRow = (double)numCol * MAX(cell.widthInFeatureSize, MIN_CELL_WIDTH) * tech.featureSize;
		} else { //if not relax the cell width
			lengthRow = (double)numCol * cell.widthInFeatureSize * tech.featureSize;
		}
		if (relaxArrayCellHeight) {  //if want to relax the cell height
			lengthCol = (double)numRow * MAX(cell.heightInFeatureSize, MIN_CELL_HEIGHT) * tech.featureSize;
		} else {  //if not relax the cell height
			lengthCol = (double)numRow * cell.heightInFeatureSize * tech.featureSize;
		}
	
	} else if (cell.memCellType == Type::RRAM) {  //if array is RRAM
		double cellHeight = cell.heightInFeatureSize;  //set RRAM cell height. The parameters are set in Cell.cpp
		double cellWidth = cell.widthInFeatureSize;  //set RRAM cell width
		if (cell.accessType == CMOS_access) {  // 1T1R
			if (relaxArrayCellWidth) {
				lengthRow = (double)numCol * MAX(cellWidth, MIN_CELL_WIDTH*2) * tech.featureSize;	// Width*2 because generally switch matrix has 2 pass gates per column, even the SL/BL driver has 2 pass gates per column in traditional 1T1R memory
			} else {
				lengthRow = (double)numCol * cellWidth * tech.featureSize;
			}
			if (relaxArrayCellHeight) {
				lengthCol = (double)numRow * MAX(cellHeight, MIN_CELL_HEIGHT) * tech.featureSize;
			} else {
				lengthCol = (double)numRow * cellHeight * tech.featureSize;
			}
		} else {	// Cross-point, if enter anything else except 'CMOS_access'
			if (relaxArrayCellWidth) {
				lengthRow = (double)numCol * MAX(cellWidth*cell.featureSize, MIN_CELL_WIDTH*2*tech.featureSize);	// Width*2 because generally switch matrix has 2 pass gates per column, even the SL/BL driver has 2 pass gates per column in traditional 1T1R memory
			} else {
				lengthRow = (double)numCol * cellWidth * cell.featureSize;
			}
			if (relaxArrayCellHeight) {
				lengthCol = (double)numRow * MAX(cellHeight*cell.featureSize, MIN_CELL_HEIGHT*tech.featureSize);
			} else {  
				lengthCol = (double)numRow * cellHeight * cell.featureSize;
			}
		}
	}
      //finish setting array size
	
	capRow1 = lengthRow * 0.2e-15/1e-6;	// BL for 1T1R, WL for Cross-point and SRAM
	capRow2 = lengthRow * 0.2e-15/1e-6;	// WL for 1T1R
	capCol = lengthCol * 0.2e-15/1e-6;
	
	resRow = lengthRow * unitWireRes; 
	resCol = lengthCol * unitWireRes;
	

	//start to initializing the subarray periphery cicuit modules
	if (cell.memCellType == Type::SRAM) {  //if array is SRAM
		
		//firstly calculate the CMOS resistance and capacitance
		resCellAccess = CalculateOnResistance(cell.widthAccessCMOS * tech.featureSize, NMOS, inputParameter.temperature, tech);
		capCellAccess = CalculateDrainCap(cell.widthAccessCMOS * tech.featureSize, NMOS, cell.widthInFeatureSize * tech.featureSize, tech);
		cell.capSRAMCell = capCellAccess + CalculateDrainCap(cell.widthSRAMCellNMOS * tech.featureSize, NMOS, cell.widthInFeatureSize * tech.featureSize, tech) + CalculateDrainCap(cell.widthSRAMCellPMOS * tech.featureSize, PMOS, cell.widthInFeatureSize * tech.featureSize, tech) + CalculateGateCap(cell.widthSRAMCellNMOS * tech.featureSize, tech) + CalculateGateCap(cell.widthSRAMCellPMOS * tech.featureSize, tech);

		senseAmp.Initialize(numCol, false, cell.minSenseVoltage, lengthRow/numCol, clkFreq, numReadCellPerOperationNeuro);
		int adderBit = (int)ceil(log2(numRow)) + avgWeightBit;	// Assume that if numRow=64 and avgWeightBit=3 then it will be a 9-bit adder
		int numAdder = numCol/numCellPerSynapse; // each cell has its own adder for SRAM

		dff.Initialize((adderBit+1)*numAdder, clkFreq);	// +1 because the adder output is 1 bit more than the input
		adder.Initialize(adderBit, numAdder);
		subtractor.Initialize(adderBit, numAdder);   // for subtracting Dummy Column 
		shiftAdd.Initialize(numAdder, adderBit+1, clkFreq, spikingMode, numReadPulse);
		
		// Note that SRAM write driver and precharger connect to BL and BL_bar in parallel (ignore the Mux that selects between these two)
		capCol += capCellAccess * numRow + senseAmp.capLoad + mux.capTgDrain * (2 + numColMuxed - 1);	// assuming Mux is already enabled
		
		wlDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numRow)), false);
		precharger.Initialize(numCol, resCol, activityColWrite,	numReadCellPerOperationNeuro, numWriteCellPerOperationNeuro);
		sramWriteDriver.Initialize(numCol, activityColWrite, numWriteCellPerOperationNeuro);
    } 
    else if (cell.memCellType == Type::RRAM) 
    {
		if (cell.accessType == CMOS_access) 
        {	// 1T1R
			if (!cell.resCellAccess)    // If not defined
				cell.resCellAccess = cell.resistanceOn * IR_DROP_TOLERANCE;    //calculate access CMOS resistance
			    cell.widthAccessCMOS = CalculateOnResistance(tech.featureSize, NMOS, inputParameter.temperature, tech) / cell.resCellAccess;   //get access CMOS width
			if (cell.widthAccessCMOS > cell.widthInFeatureSize) {	// Place transistor vertically
				printf("Transistor width of 1T1R=%.2fF is larger than the assigned cell width=%.2fF in layout\n", cell.widthAccessCMOS, cell.widthInFeatureSize);
				exit(-1);
            }

			cell.resMemCellOn = cell.resCellAccess + cell.resistanceOn;       //calculate single memory cell resistance_ON
			cell.resMemCellOff = cell.resCellAccess + cell.resistanceOff;      //calculate single memory cell resistance_OFF
			cell.resMemCellAvg = cell.resCellAccess + cell.resistanceAvg;      //calculate single memory cell resistance_AVG

			capRow2 += CalculateGateCap(cell.widthAccessCMOS * tech.featureSize, tech) * numCol;          //sum up all the gate cap of access CMOS, as the row cap
			capCol += CalculateDrainCap(cell.widthAccessCMOS * tech.featureSize, NMOS, cell.widthInFeatureSize * tech.featureSize, tech) * numRow;	// If capCol is found to be too large, increase cell.widthInFeatureSize to relax the limit

			if (digitalModeNeuro) {  //if digital mode pseudo-1T1R
				double capBL = lengthCol * 0.2e-15/1e-6;
				//for traditional parallel pseudo-1T1R, need wlDecoder, colDecoder 
                // and colDecoderDriver, MUX and muxDecoder, and VSA, Adder, ShiftAdd and DFF
				//TODO for new parallel pseudo-1T1R, sturcture will be quite different, need to change...

				wlDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numRow)), false);           //initialize WL Decoder
				colDecoderDriver.Initialize(COL_MODE, numCol, 1);          //initialize column decoder dirver, only need to drive one row
				if (numWriteColMuxed > 1) {      //if more than one column share 1 writer driver
					colDecoder.Initialize(REGULAR_COL, (int)ceil(log2(numWriteColMuxed)), false);     //initialize column decoder
				}
				// Here numColMuxed can mean how many synapses share 1 adder 
                // or how many columns share 1 S/A
				int numAdder = (int)ceil(((double)numCol/numCellPerSynapse)/numColMuxed);   // numCol is divisible by numCellPerSynapse
				int numInput = numAdder * numCellPerSynapse;        //XXX input number of MUX, 
				double resTg = cell.resMemCellOn / numRow * IR_DROP_TOLERANCE;     //transmission gate resistance
				mux.Initialize(numInput, numColMuxed, resTg, false);

				if (numColMuxed > 1) {       //if more than one column share 1 S/A
					muxDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numColMuxed)), true);
				}
				voltageSenseAmp.Initialize(numInput, clkFreq);
				//XXX adder should have X more bit to avoid over flow
				int adderBit = (int)ceil(log2(numRow)) + avgWeightBit;  // Assume that if numRow=64 and avgWeightBit=3 then it will be a 9-bit adder
				dff.Initialize((adderBit+1)*numAdder, clkFreq); // +1 because the adder output is 1 bit more than the input
				adder.Initialize(adderBit, numAdder);
				subtractor.Initialize(adderBit, numAdder);   // for subtracting Dummy Column
                
                if (parallelRead==true)
                {
                double resTg = cell.resMemCellOn * IR_DROP_TOLERANCE;
				slSwitchMatrix.Initialize(COL_MODE, numCol, resTg, activityRowRead, activityColWrite, numWriteCellPerOperationNeuro, numWritePulse, clkFreq);     //SL use switch matrix
				
				resTg = cell.resMemCellOn / numCol * IR_DROP_TOLERANCE;
				/* The previous design with seperate WL decoder and BL decoder
                blSwitchMatrix.Initialize(ROW_MODE, numRow, resTg, activityRowRead, activityColWrite, numWriteCellPerOperationNeuro, numWritePulse, clkFreq);    //BL use switch matrix
                wlDecoderOutput.Initialize(numRow); 
                 */
                wlBlSwitchMatrix.Initialize(2*numRow, activityRowRead, clkFreq, false);
                readCircuit.Initialize(readCircuitMode, (int)ceil((double)numCol/numColMuxed), maxNumIntBit, spikingMode, clkFreq);     
                }
                // need to consider input vector bit???
				shiftAdd.Initialize(numAdder, adderBit+1, clkFreq, spikingMode, numReadPulse);
			} 
            else {  //analog mode pesudo-1T1R
				wlDecoderOutput.Initialize(numRow);     //WL decoder follower    
				wlDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numRow)), false);
				
				double resTg = cell.resMemCellOn * IR_DROP_TOLERANCE;
				slSwitchMatrix.Initialize(COL_MODE, numCol, resTg, activityRowRead, activityColWrite, numWriteCellPerOperationNeuro, numWritePulse, clkFreq);     //SL use switch matrix
				
				resTg = cell.resMemCellOn / numCol * IR_DROP_TOLERANCE;
				blSwitchMatrix.Initialize(ROW_MODE, numRow, resTg, activityRowRead, activityColWrite, numWriteCellPerOperationNeuro, numWritePulse, clkFreq);    //BL use switch matrix
				
				int numInput = (int)ceil((double)numCol/numColMuxed);     //input number of mux (num of column/ num of column that share one SA)
				resTg = cell.resMemCellOn / numRow * IR_DROP_TOLERANCE;
				mux.Initialize(numInput, numColMuxed, resTg, false);

				if (numColMuxed > 1) {    //if more than one column share one SA
					muxDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numColMuxed)), true);
				}
				readCircuit.Initialize(readCircuitMode, (int)ceil((double)numCol/numColMuxed), maxNumIntBit, spikingMode, clkFreq);
				subtractor.Initialize(readCircuit.maxNumIntBit, numInput);   // for subtracting Dummy Column 
				shiftAdd.Initialize((int)ceil((double)numCol/numColMuxed), readCircuit.maxNumIntBit, clkFreq, spikingMode, numReadPulse);
			}

		} else {	// Cross-point			
			// The nonlinearity is from the selector, assuming RRAM itself is linear
			if (cell.nonlinearIV) {   //introduce nonlinearity to the RRAM resistance
				cell.resMemCellOn = cell.resistanceOn;
				cell.resMemCellOff = cell.resistanceOff;
				cell.resMemCellOnAtHalfVw = NonlinearResistance(cell.resistanceOn, cell.nonlinearity, cell.writeVoltage, cell.readVoltage, cell.writeVoltage/2);
				cell.resMemCellOffAtHalfVw = NonlinearResistance(cell.resistanceOff, cell.nonlinearity, cell.writeVoltage, cell.readVoltage, cell.writeVoltage/2);
				cell.resMemCellOnAtVw = NonlinearResistance(cell.resistanceOn, cell.nonlinearity, cell.writeVoltage, cell.readVoltage, cell.writeVoltage);
				cell.resMemCellOffAtVw = NonlinearResistance(cell.resistanceOff, cell.nonlinearity, cell.writeVoltage, cell.readVoltage, cell.writeVoltage);
				cell.resMemCellAvg = cell.resistanceAvg;
				cell.resMemCellAvgAtHalfVw = (cell.resMemCellOnAtHalfVw + cell.resMemCellOffAtHalfVw) / 2;
				cell.resMemCellAvgAtVw = (cell.resMemCellOnAtVw + cell.resMemCellOffAtVw) / 2;
			} else {  //simply assume RRAM resistance is linear
				cell.resMemCellOn = cell.resistanceOn;
				cell.resMemCellOff = cell.resistanceOff;
				cell.resMemCellOnAtHalfVw = cell.resistanceOn;
				cell.resMemCellOffAtHalfVw = cell.resistanceOff;
				cell.resMemCellOnAtVw = cell.resistanceOn;
				cell.resMemCellOffAtVw = cell.resistanceOff;
				cell.resMemCellAvg = cell.resistanceAvg;
				cell.resMemCellAvgAtHalfVw = cell.resistanceAvg;
				cell.resMemCellAvgAtVw = cell.resistanceAvg;
			}

			if (digitalModeNeuro) {  //digital mode cross-point, need wlDecoder, colDecoder and colDecoderDriver, MUX and VSA, Adder, ShiftAdd and DFF
				wlDecoderDriver.Initialize(ROW_MODE, numRow, numCol);
				wlDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numRow)), false);
				colDecoderDriver.Initialize(COL_MODE, numCol, numRow);
				if (numWriteColMuxed > 1) {    //if more than one column share 1 write driver
					colDecoder.Initialize(REGULAR_COL, (int)ceil(log2(numWriteColMuxed)), false);
				}
				// Here numColMuxed can mean how many synapses share 1 adder or how many columns share 1 S/A
				int numAdder = (int)ceil(((double)numCol/numCellPerSynapse)/numColMuxed);   // numCol is divisible by numCellPerSynapse
				int numInput = numAdder * numCellPerSynapse;    ///XXX input number of mux
				double resTg = cell.resMemCellOnAtVw / numRow * IR_DROP_TOLERANCE;
				mux.Initialize(numInput, numColMuxed, resTg, false);
				if (numColMuxed > 1) {
					muxDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numColMuxed)), true);
				}
				voltageSenseAmp.Initialize(numInput, clkFreq);
				int adderBit = (int)ceil(log2(numRow)) + avgWeightBit;  // Assume that if numRow=64 and avgWeightBit=3 then it will be a 9-bit adder
				dff.Initialize((adderBit+1)*numAdder, clkFreq); // +1 because the adder output is 1 bit more than the input
				adder.Initialize(adderBit, numAdder);
				subtractor.Initialize(adderBit, numAdder);   // for subtracting Dummy Column 
				shiftAdd.Initialize(numAdder, adderBit+1, clkFreq, spikingMode, numReadPulse);
			} else {  //analog mode cross-point
				double resTg = cell.resMemCellOnAtVw / numCol * IR_DROP_TOLERANCE;
				wlSwitchMatrix.Initialize(ROW_MODE, numRow, resTg, activityRowRead, activityColWrite, numWriteCellPerOperationNeuro, numWritePulse, clkFreq);
				
				resTg = cell.resMemCellOnAtVw / numRow * IR_DROP_TOLERANCE;
				blSwitchMatrix.Initialize(COL_MODE, numCol, resTg, activityRowRead, activityColWrite, numWriteCellPerOperationNeuro, numWritePulse, clkFreq);

				int numInput = (int)ceil((double)numCol/numColMuxed);
				resTg = cell.resMemCellOnAtVw / numRow * IR_DROP_TOLERANCE;
				mux.Initialize(numInput, numColMuxed, resTg, false);

				if (numColMuxed > 1) {
					muxDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numColMuxed)), true);
				}
				
				readCircuit.Initialize(readCircuitMode, (int)ceil((double)numCol/numColMuxed), maxNumIntBit, spikingMode, clkFreq);
				subtractor.Initialize(readCircuit.maxNumIntBit, numInput);   // for subtracting Dummy Column
				shiftAdd.Initialize((int)ceil((double)numCol/numColMuxed), readCircuit.maxNumIntBit, clkFreq, spikingMode, numReadPulse);
			}
		}
	}
   initialized = true;  //finish initialization
}

void SubArray::CalculateArea() {  //calculate layout area for total design
	if (!initialized) 
    {
		cout << "[Subarray] Error: Require initialization first!" << endl;  //ensure initialization first
	} 
    else {  //if initialized, start to do calculation
		if (cell.memCellType == Type::SRAM) {       
			
			// Array only
			heightArray = lengthCol;
			widthArray = lengthRow;
			areaArray = heightArray * widthArray;
			
			//precharger and writeDriver are always needed for all different designs
			precharger.CalculateArea(NULL, widthArray, NONE);
			sramWriteDriver.CalculateArea(NULL, widthArray, NONE);

			wlDecoder.CalculateArea(heightArray, NULL, NONE);
			senseAmp.CalculateArea(NULL, widthArray, MAGIC);
			adder.CalculateArea(NULL, widthArray, NONE);
			dff.CalculateArea(NULL, widthArray, NONE);
			subtractor.CalculateArea(NULL, widthArray, NONE);
			if (shiftAddEnable) {
				shiftAdd.CalculateArea(NULL, widthArray, NONE);
			}
			
			height = precharger.height + sramWriteDriver.height + heightArray + senseAmp.height + mux.height + adder.height + dff.height + shiftAdd.height + subtractor.height;
			width = wlDecoder.width + widthArray;
			area = height * width;
			usedArea = areaArray + wlDecoder.area + precharger.area + sramWriteDriver.area + senseAmp.area + mux.totalArea + muxDecoder.totalArea + adder.area + dff.area + shiftAdd.area + subtractor.area;
			emptyArea = area - usedArea;

	    }
        else if (cell.memCellType == Type::RRAM) {
			if (cell.accessType == CMOS_access) {	// 1T1R
				
				// Array only
				heightArray = lengthCol;
				widthArray = lengthRow;
				areaArray = heightArray * widthArray;
				
				if (digitalModeNeuro) {
          if(this->parallelRead==true)
          {
          // the component for parallel readout
          // the new WL-BL switchmatrix
          wlBlSwitchMatrix.CalculateArea(heightArray,NULL,NONE);
          slSwitchMatrix.CalculateArea(NULL, widthArray, NONE);
				  mux.CalculateArea(NULL, widthArray, NONE);
 					muxDecoder.CalculateArea(NULL, NULL, NONE);
 					double minMuxHeight = MAX(muxDecoder.height, mux.height);
 					mux.CalculateArea(minMuxHeight, widthArray, OVERRIDE);
 					subtractor.CalculateArea(NULL, widthArray, NONE);

          readCircuit.CalculateUnitArea();
					readCircuit.CalculateArea(mux.width);
					if (shiftAddEnable) {
						shiftAdd.CalculateArea(NULL, widthArray, NONE);
					}
                                               
            height = slSwitchMatrix.height + heightArray + mux.height  + readCircuit.height + shiftAdd.height + subtractor.height;
  					width = MAX(wlBlSwitchMatrix.width, muxDecoder.width) + widthArray;
  					area = height * width;
  					usedArea = areaArray + slSwitchMatrix.area + wlBlSwitchMatrix.area + mux.area + muxDecoder.area + readCircuit.area + shiftAdd.area + subtractor.area;
            //printf("usedArea is %.6e\n",usedArea);
            //printf("Array area is %.6e\n",areaArray);
            //printf("slSwitchMatrix.area is %.6e\n",slSwitchMatrix.area);
            //printf("mux.area is %.6e\n",mux.area);
            //printf("readCircuit.area %.6e\n",readCircuit.area);
  					emptyArea = area - usedArea;
          }
          else
          {
  					wlDecoder.CalculateArea(heightArray, NULL, NONE);
  					colDecoder.CalculateArea(NULL, widthArray, NONE);
  					colDecoderDriver.CalculateArea(NULL, widthArray, NONE);
                      
  					// Get Mux height, compare it with Mux decoder height, and select whichever is larger for Mux
  					mux.CalculateArea(NULL, widthArray, NONE);
  					muxDecoder.CalculateArea(NULL, NULL, NONE);
  					double minMuxHeight = MAX(muxDecoder.height, mux.height);
  					mux.CalculateArea(minMuxHeight, widthArray, OVERRIDE);
  
  					voltageSenseAmp.CalculateUnitArea();
  					voltageSenseAmp.CalculateArea(mux.widthTgShared);
  					adder.CalculateArea(NULL, widthArray, NONE);
  					dff.CalculateArea(NULL, widthArray, NONE);
  					subtractor.CalculateArea(NULL, widthArray, NONE);
  					if (shiftAddEnable) {
  						shiftAdd.CalculateArea(NULL, widthArray, NONE);
  					}
  					height = colDecoder.height + colDecoderDriver.height + heightArray + mux.height + voltageSenseAmp.height + adder.height + dff.height + shiftAdd.height + subtractor.height;
  					width = MAX(wlDecoder.width, muxDecoder.width) + widthArray;
  					area = height * width;
  					usedArea = areaArray + wlDecoder.area + colDecoder.area + colDecoderDriver.area + mux.area + voltageSenseAmp.area + muxDecoder.area + adder.area + dff.area + shiftAdd.area + subtractor.area;
  					emptyArea = area - usedArea;
                   }

				} else {  //analog mode 1T1R RRAM
					wlDecoder.CalculateArea(heightArray, NULL, NONE);
					wlDecoderOutput.CalculateArea(heightArray, NULL, NONE);
					slSwitchMatrix.CalculateArea(NULL, widthArray, NONE);
					blSwitchMatrix.CalculateArea(heightArray, NULL, NONE);

					// Get Mux height, compare it with Mux decoder height, and select whichever is larger for Mux
					mux.CalculateArea(NULL, widthArray, NONE);
					muxDecoder.CalculateArea(NULL, NULL, NONE);
					double minMuxHeight = MAX(muxDecoder.height, mux.height);
					mux.CalculateArea(minMuxHeight, widthArray, OVERRIDE);

					readCircuit.CalculateUnitArea();
					readCircuit.CalculateArea(mux.width);
					subtractor.CalculateArea(NULL, widthArray, NONE);
					if (shiftAddEnable) {
						shiftAdd.CalculateArea(NULL, mux.width, NONE);
					}

					height = slSwitchMatrix.height + heightArray + mux.height + readCircuit.height + shiftAdd.height + subtractor.height;
					width = MAX(wlDecoder.width+wlDecoderOutput.width, muxDecoder.width) + widthArray + blSwitchMatrix.width;
					area = height * width;
					usedArea = areaArray + wlDecoder.area + wlDecoderOutput.area + slSwitchMatrix.area + blSwitchMatrix.area + mux.area + readCircuit.area + voltageSenseAmp.area + muxDecoder.area + shiftAdd.area + subtractor.area;
					emptyArea = area - usedArea;
				}

			} else {        // Cross-point
				
				// Array only
				heightArray = lengthCol;
				widthArray = lengthRow;
				areaArray = heightArray * widthArray;
				
				if (digitalModeNeuro) {  //digital cross-point 
					wlDecoder.CalculateArea(heightArray, NULL, NONE);
					wlDecoderDriver.CalculateArea(heightArray, NULL, NONE);
					colDecoder.CalculateArea(NULL, widthArray, NONE);
					colDecoderDriver.CalculateArea(NULL, widthArray, NONE);

					// Get Mux height, compare it with Mux decoder height, and select whichever is larger for Mux
					mux.CalculateArea(NULL, widthArray, NONE);
					muxDecoder.CalculateArea(NULL, NULL, NONE);
					double minMuxHeight = MAX(muxDecoder.height, mux.height);
					mux.CalculateArea(minMuxHeight, widthArray, OVERRIDE);

					voltageSenseAmp.CalculateUnitArea();
					voltageSenseAmp.CalculateArea(mux.widthTgShared);
					subtractor.CalculateArea(NULL, widthArray, NONE);

					height = colDecoder.height + colDecoderDriver.height + heightArray + mux.height + voltageSenseAmp.height + adder.height + dff.height + shiftAdd.height + subtractor.height;
					width = MAX(wlDecoder.width+wlDecoderDriver.width, muxDecoder.width) + widthArray;
					area = height * width;
					usedArea = areaArray + wlDecoder.area + wlDecoderDriver.area + colDecoder.area + colDecoderDriver.area + mux.area + voltageSenseAmp.area + muxDecoder.area + adder.area + dff.area + shiftAdd.area + subtractor.area;
					emptyArea = area - usedArea;

				} else {  //analog cross-point
					wlSwitchMatrix.CalculateArea(heightArray, NULL, NONE);
					blSwitchMatrix.CalculateArea(NULL, widthArray, NONE);
					// Get Mux height, compare it with Mux decoder height, and select whichever is larger for Mux
					mux.CalculateArea(NULL, widthArray, NONE);
					muxDecoder.CalculateArea(NULL, NULL, NONE);
					double minMuxHeight = MAX(muxDecoder.height, mux.height);
					mux.CalculateArea(minMuxHeight, widthArray, OVERRIDE);
					
					readCircuit.CalculateUnitArea();
					readCircuit.CalculateArea(mux.width);
					subtractor.CalculateArea(NULL, widthArray, NONE);
					if (shiftAddEnable) {
						shiftAdd.CalculateArea(NULL, mux.width, NONE);
					}
					height = blSwitchMatrix.height + heightArray + mux.height + readCircuit.height + shiftAdd.height + subtractor.height;
					width = MAX(wlSwitchMatrix.width, muxDecoder.width) + widthArray;
					area = height * width;
					usedArea = areaArray + wlSwitchMatrix.area + blSwitchMatrix.area + mux.area + readCircuit.area + voltageSenseAmp.area + muxDecoder.area + shiftAdd.area + subtractor.area;
					emptyArea = area - usedArea;
				}
			}
		}
	}
}


void SubArray::CalculateLatency(double _rampInput) {   //calculate latency for different mode 
	if (!initialized) {
		cout << "[Subarray] Error: Require initialization first!" << endl;
	} else {
		
		if (cell.memCellType == Type::SRAM) {
			
			int numReadOperationPerRow = (int)ceil((double)numCol/numReadCellPerOperationNeuro);
			int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);
			wlDecoder.CalculateLatency(1e20, capRow1, NULL, numRow*numReadPulse*activityRowRead, numRow*activityRowWrite);
			precharger.CalculateLatency(1e20, capCol, numReadOperationPerRow*numRow*numReadPulse*activityRowRead, numWriteOperationPerRow*numRow*activityRowWrite);
			sramWriteDriver.CalculateLatency(1e20, capCol, resCol, numWriteOperationPerRow*numRow*activityRowWrite);
			senseAmp.CalculateLatency(numReadOperationPerRow*numRow*numReadPulse*activityRowRead);
			adder.CalculateLatency(1e20, dff.capTgDrain, numReadOperationPerRow*numRow*numReadPulse*activityRowRead);
			dff.CalculateLatency(1e20, numReadOperationPerRow*numRow*numReadPulse*activityRowRead);
			subtractor.CalculateLatency(1e20, dff.capTgDrain, numReadOperationPerRow*numReadPulse);
			if (shiftAddEnable) {
				shiftAdd.CalculateLatency(numReadPulse);	// There are numReadPulse times of shift-and-add
			}

			// Read
			double resPullDown = CalculateOnResistance(cell.widthSRAMCellNMOS * tech.featureSize, NMOS, inputParameter.temperature, tech);
			double tau = (resCellAccess + resPullDown) * (capCellAccess + capCol) + resCol * capCol / 2;
			tau *= log(tech.vdd / (tech.vdd - cell.minSenseVoltage / 2));   /* one signal raises and the other drops, so cell.minSenseVoltage/2 is enough */
			double gm = CalculateTransconductance(cell.widthAccessCMOS * tech.featureSize, NMOS, tech);
			double beta = 1 / (resPullDown * gm);
			double colRamp = 0;
			colDelay = horowitz(tau, beta, wlDecoder.rampOutput, &colRamp) * numReadOperationPerRow * numRow * numReadPulse * activityRowRead;

			readLatency += wlDecoder.readLatency;
			readLatency += precharger.readLatency;
			readLatency += colDelay;
			readLatency += senseAmp.readLatency;
			readLatency += subtractor.readLatency;
			readLatency += adder.readLatency;
			readLatency += dff.readLatency;
			readLatency += shiftAdd.readLatency;

			// Write (assume the average delay of pullup and pulldown inverter in SRAM cell)
			double resPull;
			resPull = (CalculateOnResistance(cell.widthSRAMCellNMOS * tech.featureSize, NMOS, inputParameter.temperature, tech) + CalculateOnResistance(cell.widthSRAMCellPMOS * tech.featureSize, PMOS, inputParameter.temperature, tech)) / 2;    // take average
			tau = resPull * cell.capSRAMCell;
			gm = (CalculateTransconductance(cell.widthSRAMCellNMOS * tech.featureSize, NMOS, tech) + CalculateTransconductance(cell.widthSRAMCellPMOS * tech.featureSize, PMOS, tech)) / 2;   // take average
			beta = 1 / (resPull * gm);

			writeLatency += horowitz(tau, beta, 1e20, NULL) * numWriteOperationPerRow * numRow * activityRowWrite;

			writeLatency += wlDecoder.writeLatency;
			writeLatency += precharger.writeLatency;
			writeLatency += sramWriteDriver.writeLatency;

	    } else if (cell.memCellType == Type::RRAM) {
			if (cell.accessType == CMOS_access) {   // 1T1R
				
				if (digitalModeNeuro) { //1T1R for digital neuro
					double capBL = lengthCol * 0.2e-15/1e-6;
					wlDecoder.CalculateLatency(1e20, capRow2, NULL, numRow*activityRowRead*numReadPulse*numColMuxed, numRow*activityRowWrite);
					double colDecoderLoad = (colDecoderDriver.capInvInput + colDecoderDriver.capTgGateN * 2 + colDecoderDriver.capTgGateP) * numWriteCellPerOperationNeuro;
					colDecoder.CalculateLatency(1e20, colDecoderLoad, NULL, 1, numRow*activityRowWrite*numWriteColMuxed);	// Doesn't matter for read
					if (colDecoder.rampOutput > 0) {
						colDecoderDriver.CalculateLatency(colDecoder.rampOutput, capCol, capBL, resCol, 1, numRow*activityRowWrite*numWriteColMuxed*2);	// Doesn't matter for read. *2 means 2-step write
					} else {    // The case where column decoder doesn't exist
						colDecoderDriver.CalculateLatency(1e20, capCol, capBL, resCol, 1, numRow*activityRowWrite*numWriteColMuxed*2);  // Doesn't matter for read. *2 means 2-step write
					}

					// Calculate column latency
					double colRamp = 0;
					double tau = resCol * capBL / 2 * (cell.resMemCellOff + resCol / 3) / (cell.resMemCellOff + resCol);
					colDelay = horowitz(tau, 0, 1e20, &colRamp);	// Just to generate colRamp
					
					// Read
					mux.CalculateLatency(colRamp, 0, 1);
					// Here numColMuxed can mean how many synapses share 1 adder or how many columns share 1 S/A
					int numAdder = (int)ceil((numCol/numCellPerSynapse)/numColMuxed);   // numCol is divisible by numCellPerSynapse
					int numInput = numAdder * numCellPerSynapse;
					muxDecoder.CalculateLatency(1e20, mux.capTgGateN*numInput, mux.capTgGateP*numInput, 1, 1);
					double capInputLoad = capBL + mux.capTgDrain * (2 + numColMuxed - 1);
					voltageSenseAmp.CalculateLatency(capInputLoad, numColMuxed*numRow*numReadPulse*activityRowRead);
					adder.CalculateLatency(1e20, dff.capTgDrain, numColMuxed*numRow*numReadPulse*activityRowRead);
					dff.CalculateLatency(1e20, numColMuxed*numRow*numReadPulse*activityRowRead);
					subtractor.CalculateLatency(1e20, dff.capTgDrain, numColMuxed*numReadPulse);
					if (shiftAddEnable) {
						shiftAdd.CalculateLatency(numReadPulse);	// There are numReadPulse times of shift-and-add
					}


					readLatency += MAX(wlDecoder.readLatency, muxDecoder.readLatency + mux.readLatency);
					readLatency += voltageSenseAmp.readLatency;
					readLatency += subtractor.readLatency;
					readLatency += adder.readLatency;
					readLatency += dff.readLatency;
					readLatency += shiftAdd.readLatency;

					// Write
					writeLatency += MAX(wlDecoder.writeLatency, colDecoder.writeLatency + colDecoderDriver.writeLatency);
				} else {
					int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);
					wlDecoder.CalculateLatency(1e20, wlDecoderOutput.capNorInput, NULL, 1, numRow*activityRowWrite);
					wlDecoderOutput.CalculateLatency(wlDecoder.rampOutput, capRow2, resRow, 1, numRow*activityRowWrite);
					slSwitchMatrix.CalculateLatency(1e20, capCol, resCol, 1, 1);
					blSwitchMatrix.CalculateLatency(1e20, capRow1, resRow, numReadPulse, maxNumWritePulse*2*numWriteOperationPerRow*numRow*activityRowWrite);	// *2 means 2-step write

					// Calculate column latency
					double colRamp = 0;
					double tau = resCol * capCol / 2 * (cell.resMemCellOff + resCol / 3) / (cell.resMemCellOff + resCol);
					colDelay = horowitz(tau, 0, blSwitchMatrix.rampOutput, &colRamp);	// Just to generate colRamp

					mux.CalculateLatency(colRamp, 0, 1);
					int numInput = (int)ceil((double)numCol/numColMuxed);
					muxDecoder.CalculateLatency(1e20, mux.capTgGateN*numInput, mux.capTgGateP*numInput, 1, 1);
					
					// Read
					if (readCircuit.mode == CMOS) {
						double Cin = capCol + mux.capTgDrain * (2 + numColMuxed - 1) + readCircuit.capTgDrain + readCircuit.capPmosGate;
						double Imax = numRow * cell.readVoltage / cell.resMemCellOn;
						cell.readPulseWidth = Cin * readCircuit.voltageIntThreshold / Imax * readCircuit.maxNumIntPerCycle;
					} else {	// mode==OSCILLATION
						double Cin = capCol + mux.capTgDrain * (2 + numColMuxed - 1) + readCircuit.capInvInput;
						double Rmin = cell.resMemCellOn / numRow;
						double Rp = 1 / (1/Rmin + 1/readCircuit.R_OSC_OFF);
						double t_rise = -Rp * Cin * log((readCircuit.Vth - readCircuit.Vrow*Rp/Rmin) / (readCircuit.Vhold - readCircuit.Vrow*Rp/Rmin));
						cell.readPulseWidth = t_rise * readCircuit.maxNumIntPerCycle;
					}
					readCircuit.CalculateLatency(numColMuxed*numReadPulse);
					subtractor.CalculateLatency(1e20, 0, numColMuxed*numReadPulse);
					if (shiftAddEnable) {
						shiftAdd.CalculateLatency(numReadPulse);
					}
					
					readLatency += wlDecoderOutput.readLatency;
					readLatency += blSwitchMatrix.readLatency;
					readLatency += readCircuit.readLatency;
					readLatency += subtractor.readLatency;
					readLatency += shiftAdd.readLatency;
					
					// Write
					writeLatency += wlDecoder.writeLatency;
					writeLatency += wlDecoderOutput.writeLatency;
					writeLatency += blSwitchMatrix.writeLatency;
				}

			} else {	// Cross-point
				if (digitalModeNeuro) {
					double numReadCells = (int)ceil((double)numCol/numColMuxed);    // similar parameter as numReadCellPerOperationNeuro, which is for SRAM
					double wlDecoderLoad = colDecoderDriver.capInvInput + colDecoderDriver.capTgGateN + colDecoderDriver.capTgGateP;
					wlDecoder.CalculateLatency(1e20, wlDecoderLoad, NULL, numRow*activityRowRead*numReadPulse*numColMuxed, numRow*activityRowWrite);
					wlDecoderDriver.CalculateLatency(wlDecoder.rampOutput, capRow1, capRow1, resRow, numRow*activityRowRead*numReadPulse*numColMuxed, numRow*activityRowWrite);
					double colDecoderLoad = (colDecoderDriver.capInvInput + colDecoderDriver.capTgGateN + colDecoderDriver.capTgGateP) * numWriteCellPerOperationNeuro;
					colDecoder.CalculateLatency(1e20, colDecoderLoad, NULL, 1, numRow*activityRowWrite*numWriteColMuxed);  // Doesn't matter for read
					if (colDecoder.rampOutput > 0) {
						colDecoderDriver.CalculateLatency(colDecoder.rampOutput, capCol, capCol, resCol, 1, numRow*activityRowWrite*numWriteColMuxed*2); // Doesn't matter for read. *2 means 2-step write
					} else {    // The case where column decoder doesn't exist
						colDecoderDriver.CalculateLatency(1e20, capCol, capCol, resCol, 1, numRow*activityRowWrite*numWriteColMuxed*2);  // Doesn't matter for read. *2 means 2-step write
					}
					
					// Calculate column latency
					double colRamp = 0;
					double tau = resCol * capCol / 2 * (cell.resMemCellOff + resCol / 3) / (cell.resMemCellOff + resCol);
					colDelay = horowitz(tau, 0, 1e20, &colRamp);	// Just to generate colRamp
					
					// Read
					mux.CalculateLatency(colRamp, 0, 1);
					// Here numColMuxed can mean how many synapses share 1 adder or how many columns share 1 S/A
					int numAdder = (int)ceil(((double)numCol/numCellPerSynapse)/numColMuxed);   // numCol is divisible by numCellPerSynapse
					int numInput = numAdder * numCellPerSynapse;
					muxDecoder.CalculateLatency(1e20, mux.capTgGateN*numInput, mux.capTgGateP*numInput, 1, 1);
					double capInputLoad = capCol + mux.capTgDrain * (2 + numColMuxed - 1);
					voltageSenseAmp.CalculateLatency(capInputLoad, numColMuxed*numRow*numReadPulse*activityRowRead);
					adder.CalculateLatency(1e20, dff.capTgDrain, numColMuxed*numRow*numReadPulse*activityRowRead);
					dff.CalculateLatency(1e20, numColMuxed*numRow*numReadPulse*activityRowRead);
					subtractor.CalculateLatency(1e20, dff.capTgDrain, numColMuxed*numReadPulse);
					if (shiftAddEnable) {
						shiftAdd.CalculateLatency(numReadPulse);	// There are numReadPulse times of shift-and-add
					}

					readLatency += MAX(wlDecoder.readLatency + wlDecoderDriver.readLatency, muxDecoder.readLatency + muxDecoder.readLatency);
					readLatency += voltageSenseAmp.readLatency;
					readLatency += subtractor.readLatency;
					readLatency += adder.readLatency;
					readLatency += dff.readLatency;
					readLatency += shiftAdd.readLatency;

					// Write
					writeLatency += MAX(wlDecoder.writeLatency + wlDecoderDriver.writeLatency, colDecoder.writeLatency + colDecoderDriver.writeLatency);
				} else {
					int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);
					wlSwitchMatrix.CalculateLatency(1e20, capRow1, resRow, numReadPulse, maxNumWritePulse*2*numWriteOperationPerRow*numRow*activityRowWrite);	// *2 means 2-step write
					blSwitchMatrix.CalculateLatency(1e20, capCol, resCol, 1, 1);

					// Calculate column latency
					double colRamp = 0;
					double tau = resCol * capCol / 2 * (cell.resMemCellOff + resCol / 3) / (cell.resMemCellOff + resCol);
					colDelay = horowitz(tau, 0, wlSwitchMatrix.rampOutput, &colRamp);	// Just to generate colRamp

					mux.CalculateLatency(colRamp, 0, 1);
					int numInput = (int)ceil((double)numCol/numColMuxed);
					muxDecoder.CalculateLatency(1e20, mux.capTgGateN*numInput, mux.capTgGateP*numInput, 1, 1);

					// Read
					if (readCircuit.mode == CMOS) {
						double Cin = capCol + mux.capTgDrain * (2 + numColMuxed - 1) + readCircuit.capTgDrain + readCircuit.capPmosGate;
						double Imax = numRow * cell.readVoltage / cell.resMemCellOn;
						cell.readPulseWidth = Cin * readCircuit.voltageIntThreshold / Imax * readCircuit.maxNumIntPerCycle;
					} else {    // mode==OSCILLATION
						double Cin = capCol + mux.capTgDrain * (2 + numColMuxed - 1) + readCircuit.capInvInput;
						double Rmin = cell.resMemCellOn / numRow;
						double Rp = 1 / (1/Rmin + 1/readCircuit.R_OSC_OFF);
						double t_rise = -Rp * Cin * log((readCircuit.Vth - readCircuit.Vrow*Rp/Rmin) / (readCircuit.Vhold - readCircuit.Vrow*Rp/Rmin));
						cell.readPulseWidth = t_rise * readCircuit.maxNumIntPerCycle;
					}
					readCircuit.CalculateLatency(numColMuxed*numReadPulse);
					subtractor.CalculateLatency(1e20, 0, numColMuxed*numReadPulse);
					if (shiftAddEnable) {
						shiftAdd.CalculateLatency(numReadPulse);
					}
					
					readLatency += wlSwitchMatrix.readLatency;
					readLatency += readCircuit.readLatency;
					readLatency += subtractor.readLatency;
					readLatency += shiftAdd.readLatency;
					
					// Write
					writeLatency += wlSwitchMatrix.writeLatency;
				}
			}
		}
	}
}

void SubArray::CalculatePower() {
	if (!initialized) {
		cout << "[Subarray] Error: Require initialization first!" << endl;
	} else {
		if (cell.memCellType == Type::SRAM) {
			
			double numReadOperationPerRow;   // average value (can be non-integer for energy calculation)
			if (numCol > numReadCellPerOperationNeuro)
				numReadOperationPerRow = numCol / numReadCellPerOperationNeuro;
			else
				numReadOperationPerRow = 1;

			double numWriteOperationPerRow;   // average value (can be non-integer for energy calculation)
			if (numCol * activityColWrite > numWriteCellPerOperationNeuro)
				numWriteOperationPerRow = numCol * activityColWrite / numWriteCellPerOperationNeuro;
			else
				numWriteOperationPerRow = 1;

			wlDecoder.CalculatePower(numRow*numReadPulse*activityRowRead, numRow*activityRowWrite);
			precharger.CalculatePower(numReadOperationPerRow*numRow*numReadPulse*activityRowRead, numWriteOperationPerRow*numRow*activityRowWrite);
			sramWriteDriver.CalculatePower(numWriteOperationPerRow*numRow*activityRowWrite);
			senseAmp.CalculatePower(numReadOperationPerRow*numRow*numReadPulse*activityRowRead);
			adder.CalculatePower(numReadOperationPerRow*numRow*numReadPulse*activityRowRead, numReadCellPerOperationNeuro/numCellPerSynapse);				
			dff.CalculatePower(numReadOperationPerRow*numRow*numReadPulse*activityRowRead, numReadCellPerOperationNeuro/numCellPerSynapse*(adder.numBit+1));	// +1 because the adder output is 1 bit more than the input
			subtractor.CalculatePower(numReadOperationPerRow*numReadPulse, numReadCellPerOperationNeuro/numCellPerSynapse);	
			if (shiftAddEnable) {
				shiftAdd.CalculatePower(numReadPulse);
			}

			// Array
			readDynamicEnergyArray = 0; // Just BL discharging
			writeDynamicEnergyArray = cell.capSRAMCell * tech.vdd * tech.vdd * 2 * numCol * activityColWrite * numRow * activityRowWrite;    // flip Q and Q_bar

			// Read
			readDynamicEnergy += wlDecoder.readDynamicEnergy;
			readDynamicEnergy += precharger.readDynamicEnergy;
			readDynamicEnergy += readDynamicEnergyArray;
			readDynamicEnergy += senseAmp.readDynamicEnergy;
			readDynamicEnergy += subtractor.readDynamicEnergy;
			readDynamicEnergy += adder.readDynamicEnergy;
			readDynamicEnergy += dff.readDynamicEnergy;
			readDynamicEnergy += shiftAdd.readDynamicEnergy;

			// Write
			writeDynamicEnergy += wlDecoder.writeDynamicEnergy;
			writeDynamicEnergy += precharger.writeDynamicEnergy;
			writeDynamicEnergy += sramWriteDriver.writeDynamicEnergy;
			writeDynamicEnergy += writeDynamicEnergyArray;
			
			// Array leakage (assume 2 INV)
			leakage += CalculateGateLeakage(INV, 1, cell.widthSRAMCellNMOS * tech.featureSize,
					cell.widthSRAMCellPMOS * tech.featureSize, inputParameter.temperature, tech) * tech.vdd * 2;
			leakage *= numRow * numCol;

			// Leakage
			leakage += wlDecoder.leakage;
			leakage += precharger.leakage;
			leakage += sramWriteDriver.leakage;
			leakage += senseAmp.leakage;
			leakage += mux.leakage;
			leakage += muxDecoder.leakage;
			leakage += dff.leakage;
			leakage += adder.leakage;
			leakage += shiftAdd.leakage;
			leakage += subtractor.leakage;
			
	    } else if (cell.memCellType == Type::RRAM) {
			if (cell.accessType == CMOS_access) {   // 1T1R
							
				if (digitalModeNeuro) {
					double numReadCells = (int)ceil((double)numCol/numColMuxed);    // similar parameter as numReadCellPerOperationNeuro, which is for SRAM
					int numWriteCells = (int)ceil((double)numCol/numWriteColMuxed);
					wlDecoder.CalculatePower(numRow*activityRowRead*numReadPulse*numColMuxed, numRow*activityRowWrite);
					colDecoder.CalculatePower(1, numRow*activityRowWrite*numWriteColMuxed);  // Doesn't matter for read
					colDecoderDriver.CalculatePower(numReadCells, numWriteCells, 1, numRow*activityRowWrite*numWriteColMuxed*2); // Doesn't matter for read. *2 means 2-step write
					mux.CalculatePower(numRow*activityRowRead*numColMuxed*numReadPulse);	// Mux still consumes energy during row-by-row read
					muxDecoder.CalculatePower(numColMuxed*numReadPulse, 1); //
					voltageSenseAmp.CalculatePower(numRow*activityRowRead*numColMuxed*numReadPulse);
					adder.CalculatePower(numColMuxed*numRow*numReadPulse*activityRowRead, numReadCells/numCellPerSynapse);
					dff.CalculatePower(numColMuxed*numRow*numReadPulse*activityRowRead, numReadCells/numCellPerSynapse*(adder.numBit+1));    // +1 because the adder output is 1 bit more than the input
					subtractor.CalculatePower(numColMuxed*numReadPulse, numReadCells/numCellPerSynapse);
					if (shiftAddEnable) {
						shiftAdd.CalculatePower(numReadPulse);
					}

					// Array (here it is regular 1T1R array)
					double capBL = lengthCol * 0.2e-15/1e-6;
					readDynamicEnergyArray += capBL * cell.readVoltage * cell.readVoltage * numReadCells; // Selected BLs
					readDynamicEnergyArray += capRow2 * tech.vdd * tech.vdd; // Selected WL
					// Use average case in write energy calculation: half SET and half RESET, and Ron and Roff are also half and half in each operation
					readDynamicEnergyArray *= numRow * activityRowRead * numReadPulse * numColMuxed;
					// SET
					writeDynamicEnergyArray += capBL * cell.writeVoltage * cell.writeVoltage * MIN(numWriteCells, numCol)/2;
					writeDynamicEnergyArray += cell.writeVoltage * cell.writeVoltage * (1/cell.resMemCellOn + 1/cell.resMemCellOff) / 2 * MIN(numWriteCells, numCol) / 2 * cell.writePulseWidth; // half SET with Ron and Roff 50/50
					// RESET
					writeDynamicEnergyArray += capCol * cell.writeVoltage * cell.writeVoltage * MIN(numWriteCells, numCol)/2;
					writeDynamicEnergyArray += cell.writeVoltage * cell.writeVoltage * (1/cell.resMemCellOn + 1/cell.resMemCellOff) / 2 * MIN(numWriteCells, numCol) / 2 * cell.writePulseWidth;    // half RESET with Ron and Roff 50/50
					writeDynamicEnergyArray *= numRow * activityRowWrite * numWriteColMuxed;

					// Read
					readDynamicEnergy += wlDecoder.readDynamicEnergy;
					readDynamicEnergy += readDynamicEnergyArray;
					readDynamicEnergy += mux.readDynamicEnergy;
					readDynamicEnergy += muxDecoder.readDynamicEnergy;
					readDynamicEnergy += voltageSenseAmp.readDynamicEnergy;
					readDynamicEnergy += subtractor.readDynamicEnergy;
					readDynamicEnergy += adder.readDynamicEnergy;
					readDynamicEnergy += dff.readDynamicEnergy;
					readDynamicEnergy += shiftAdd.readDynamicEnergy;

					// Write
					writeDynamicEnergy += wlDecoder.writeDynamicEnergy;
					writeDynamicEnergy += colDecoder.writeDynamicEnergy;
					writeDynamicEnergy += colDecoderDriver.writeDynamicEnergy;
					writeDynamicEnergy += writeDynamicEnergyArray;

				} else {
					double numWriteOperationPerRow;   // average value (can be non-integer for energy calculation)
					if (numCol * activityColWrite > numWriteCellPerOperationNeuro)
						numWriteOperationPerRow = numCol * activityColWrite / numWriteCellPerOperationNeuro;
					else
						numWriteOperationPerRow = 1;
					wlDecoder.CalculatePower(numReadPulse, numRow*activityRowWrite);
					wlDecoderOutput.CalculatePower(numReadPulse, numRow*activityRowWrite);
					blSwitchMatrix.CalculatePower(numReadPulse, numRow*activityRowWrite);
					slSwitchMatrix.CalculatePower(1, numWriteOperationPerRow*numRow*activityRowWrite);
					mux.CalculatePower(numColMuxed*numReadPulse);
					muxDecoder.CalculatePower(numColMuxed*numReadPulse, 1);
					readCircuit.CalculatePower(numColMuxed*numReadPulse);
					subtractor.CalculatePower(numColMuxed*numReadPulse, numCol/numColMuxed);
					if (shiftAddEnable) {
						shiftAdd.CalculatePower(numReadPulse);
					}
					
					// Array
					readDynamicEnergyArray += capRow1 * readCircuit.Vrow * readCircuit.Vrow * numRow * activityRowRead;   // Selected BLs
					readDynamicEnergyArray += capCol * readCircuit.Vcol * readCircuit.Vcol * numCol;	// Read all columns in total
					readDynamicEnergyArray += capRow2 * tech.vdd * tech.vdd * numRow; // All WLs open
					readDynamicEnergyArray += cell.readVoltage * cell.readVoltage / cell.resMemCellAvg * cell.readPulseWidth * numRow * activityRowRead * numCol; // Unselected SLs are floating (assume the read integration threshold is small enough)
					readDynamicEnergyArray *= numReadPulse;
					// Use average case in write energy calculation: half LTP and half LTD with average resistance
					// LTP
					writeDynamicEnergyArray += capCol * cell.writeVoltage * cell.writeVoltage * MIN(numWriteCellPerOperationNeuro, numCol * activityColWrite) / 2 * numWritePulse;	// Selected SLs
					writeDynamicEnergyArray += capCol * cell.writeVoltage * cell.writeVoltage * (numCol - MIN(numWriteCellPerOperationNeuro, numCol*activityColWrite)/2);	// Unselected SLs
					writeDynamicEnergyArray += capRow1 * cell.writeVoltage * cell.writeVoltage;	// Selected BL
					writeDynamicEnergyArray += capRow2 * tech.vdd * tech.vdd;	// Selected WL
					writeDynamicEnergyArray += cell.writeVoltage * cell.writeVoltage / cell.resMemCellAvg * MIN(numWriteCellPerOperationNeuro, numCol*activityColWrite) / 2 * numWritePulse * cell.writePulseWidth;	// LTP
					// LTD
					writeDynamicEnergyArray += capCol * cell.writeVoltage * cell.writeVoltage * MIN(numWriteCellPerOperationNeuro, numCol * activityColWrite) / 2 * numWritePulse;    // Selected SLs
					writeDynamicEnergyArray += capRow2 * tech.vdd * tech.vdd; // Selected WL
					writeDynamicEnergyArray += cell.writeVoltage * cell.writeVoltage / cell.resMemCellAvg * MIN(numWriteCellPerOperationNeuro, numCol*activityColWrite) / 2 * numWritePulse * cell.writePulseWidth;	// LTD
					writeDynamicEnergyArray *= numWriteOperationPerRow * numRow * activityRowWrite;

					// Read
					readDynamicEnergy += wlDecoder.readDynamicEnergy;
					readDynamicEnergy += wlDecoderOutput.readDynamicEnergy;
					readDynamicEnergy += blSwitchMatrix.readDynamicEnergy;
					readDynamicEnergy += readDynamicEnergyArray;
					readDynamicEnergy += mux.readDynamicEnergy;
					readDynamicEnergy += muxDecoder.readDynamicEnergy;
					readDynamicEnergy += readCircuit.readDynamicEnergy;
					readDynamicEnergy += subtractor.readDynamicEnergy;
					readDynamicEnergy += shiftAdd.readDynamicEnergy;
					
					// Write
					writeDynamicEnergy += wlDecoder.writeDynamicEnergy;
					writeDynamicEnergy += wlDecoderOutput.writeDynamicEnergy;
					writeDynamicEnergy += blSwitchMatrix.writeDynamicEnergy;
					writeDynamicEnergy += slSwitchMatrix.writeDynamicEnergy;
					writeDynamicEnergy += writeDynamicEnergyArray;
				}
				
				// Leakage
				leakage += wlDecoder.leakage;
				leakage += wlDecoderOutput.leakage;
				leakage += colDecoder.leakage;
				leakage += colDecoderDriver.leakage;
				leakage += blSwitchMatrix.leakage;
				leakage += slSwitchMatrix.leakage;
				leakage += mux.leakage;
				leakage += muxDecoder.leakage;
				leakage += readCircuit.leakage;
				leakage += voltageSenseAmp.leakage;
				leakage += shiftAdd.leakage;
				leakage += subtractor.leakage;

			} else {	// Cross-point

				if (digitalModeNeuro) {
					double numReadCells = (int)ceil((double)numCol/numColMuxed);    // similar parameter as numReadCellPerOperationNeuro, which is for SRAM
					int numWriteCells = (int)ceil((double)numCol/numWriteColMuxed);
					wlDecoder.CalculatePower(numRow*activityRowRead*numReadPulse*numColMuxed, numRow*activityRowWrite);
					wlDecoderDriver.CalculatePower(numReadCells, numWriteCells, numRow*activityRowRead*numReadPulse*numColMuxed, numRow*activityRowWrite*numWriteColMuxed);
					colDecoder.CalculatePower(1, numRow*activityRowWrite*numWriteColMuxed);	// Doesn't matter for read
					colDecoderDriver.CalculatePower(numReadCells, numWriteCells, 1, numRow*activityRowWrite*numWriteColMuxed*2);	// Doesn't matter for read. *2 means 2-step write
					mux.CalculatePower(numRow*activityRowRead*numColMuxed*numReadPulse);	// Mux still consumes energy during row-by-row read
					muxDecoder.CalculatePower(numColMuxed*numReadPulse, 1);
					voltageSenseAmp.CalculatePower(numRow*activityRowRead*numColMuxed*numReadPulse);
					adder.CalculatePower(numColMuxed*numRow*numReadPulse*activityRowRead, numReadCells/numCellPerSynapse);
					dff.CalculatePower(numColMuxed*numRow*numReadPulse*activityRowRead, numReadCells/numCellPerSynapse*(adder.numBit+1));    // +1 because the adder output is 1 bit more than the input
					subtractor.CalculatePower(numColMuxed*numReadPulse, numReadCells/numCellPerSynapse);
					if (shiftAddEnable) {
						shiftAdd.CalculatePower(numReadPulse);
					}

					// Array
					readDynamicEnergyArray += capRow1 * cell.readVoltage * cell.readVoltage * (numRow-1);   // All WLs except the selected one
					readDynamicEnergyArray += capCol * cell.readVoltage * cell.readVoltage * numReadCells;   // Selected BLs
					readDynamicEnergyArray *= numRow * activityRowWrite * numColMuxed * numReadPulse;
					// Use average case in write energy calculation: half SET and half RESET, and Ron and Roff are also half and half in each operation
					// SET
					writeDynamicEnergyArray += capRow1 * cell.writeVoltage * cell.writeVoltage;   // Selected WL
					writeDynamicEnergyArray += cell.writeVoltage * cell.writeVoltage * (1/cell.resMemCellOnAtVw + 1/cell.resMemCellOffAtVw) / 2 * MIN(numWriteCells, numCol) / 2 * cell.writePulseWidth;  // half SET with Ron and Roff 50/50
					writeDynamicEnergyArray += cell.writeVoltage/2 * cell.writeVoltage/2 * (1/cell.resMemCellOnAtHalfVw + 1/cell.resMemCellOffAtHalfVw) / 2 * (numCol - MIN(numWriteCells, numCol)/2) * cell.writePulseWidth;    // Half-selected cells on the row
					writeDynamicEnergyArray += cell.writeVoltage/2 * cell.writeVoltage/2 * (1/cell.resMemCellOnAtHalfVw + 1/cell.resMemCellOffAtHalfVw) / 2 * (numRow-1) * MIN(numWriteCells, numCol) / 2 * cell.writePulseWidth;   // Half-selected cells on the selected columns
					// RESET
					writeDynamicEnergyArray += cell.writeVoltage * cell.writeVoltage * (1/cell.resMemCellOnAtVw + 1/cell.resMemCellOffAtVw) / 2 * MIN(numWriteCells, numCol) / 2 * cell.writePulseWidth;  // half RESET with Ron and Roff 50/50
					writeDynamicEnergyArray += cell.writeVoltage/2 * cell.writeVoltage/2 * (1/cell.resMemCellOnAtHalfVw + 1/cell.resMemCellOffAtHalfVw) / 2 * (numCol - MIN(numWriteCells, numCol)/2) * cell.writePulseWidth;    // Half-selected cells on the row
					writeDynamicEnergyArray += cell.writeVoltage/2 * cell.writeVoltage/2 * (1/cell.resMemCellOnAtHalfVw + 1/cell.resMemCellOffAtHalfVw) / 2 * (numRow-1) * MIN(numWriteCells, numCol) / 2 * cell.writePulseWidth;   // Half-selected cells on the selected columns
					// Both SET and RESET
					writeDynamicEnergyArray += capCol * cell.writeVoltage/2 * cell.writeVoltage/2 * numCol; // Unselected BLs (every BL has one time to charge to V/2 within the 2-step write)
					writeDynamicEnergyArray += capRow1 * cell.writeVoltage/2 * cell.writeVoltage/2 * (numRow-1);  // Unselected WLs
					writeDynamicEnergyArray *= numRow * activityRowWrite * numWriteColMuxed;

					// Read
					readDynamicEnergy += wlDecoder.readDynamicEnergy;
					readDynamicEnergy += wlDecoderDriver.readDynamicEnergy;
					readDynamicEnergy += readDynamicEnergyArray;
					readDynamicEnergy += mux.readDynamicEnergy;
					readDynamicEnergy += muxDecoder.readDynamicEnergy;
					readDynamicEnergy += voltageSenseAmp.readDynamicEnergy;
					readDynamicEnergy += subtractor.readDynamicEnergy;
					readDynamicEnergy += adder.readDynamicEnergy;
					readDynamicEnergy += dff.readDynamicEnergy;
					readDynamicEnergy += shiftAdd.readDynamicEnergy;

					// Write
					writeDynamicEnergy += wlDecoder.writeDynamicEnergy;
					writeDynamicEnergy += wlDecoderDriver.writeDynamicEnergy;
					writeDynamicEnergy += colDecoder.writeDynamicEnergy;
					writeDynamicEnergy += colDecoderDriver.writeDynamicEnergy;
					writeDynamicEnergy += writeDynamicEnergyArray;
				} else {

					double numWriteOperationPerRow;   // average value (can be non-integer for energy calculation)
					if (numCol * activityColWrite > numWriteCellPerOperationNeuro)
						numWriteOperationPerRow = numCol * activityColWrite / numWriteCellPerOperationNeuro;
					else
						numWriteOperationPerRow = 1;
					wlSwitchMatrix.CalculatePower(numReadPulse, numRow*activityRowWrite);
					blSwitchMatrix.CalculatePower(1, numWriteOperationPerRow*numRow*activityRowWrite);
					mux.CalculatePower(numColMuxed*numReadPulse);
					muxDecoder.CalculatePower(numColMuxed*numReadPulse, 1);
					readCircuit.CalculatePower(numColMuxed*numReadPulse);
					subtractor.CalculatePower(numColMuxed*numReadPulse, numCol/numColMuxed);
					if (shiftAddEnable) {
						shiftAdd.CalculatePower(numReadPulse);
					}
					
					// Array
					readDynamicEnergyArray += capRow1 * readCircuit.Vrow * readCircuit.Vrow * numRow * activityRowRead;   // Selected BLs
					readDynamicEnergyArray += capCol * readCircuit.Vcol * readCircuit.Vcol * numCol;    // Read all columns in total
					readDynamicEnergyArray += cell.readVoltage * cell.readVoltage / cell.resMemCellAvg * cell.readPulseWidth * numRow * activityRowRead * numCol; // Unselected SLs are floating (assume the read integration threshold is small enough)
					readDynamicEnergyArray *= numReadPulse;
					// Use average case in write energy calculation: half LTP and half LTD with average resistance
					double totalWriteTime = cell.writePulseWidth * maxNumWritePulse;
					// LTP
					writeDynamicEnergyArray += capRow1 * cell.writeVoltage * cell.writeVoltage;	// Selected WL
					writeDynamicEnergyArray += capCol * cell.writeVoltage * cell.writeVoltage * MIN(numWriteCellPerOperationNeuro, numCol*activityColWrite) / 2 * numWritePulse;	// Selected BLs
					writeDynamicEnergyArray += cell.writeVoltage * cell.writeVoltage / cell.resMemCellAvgAtVw * MIN(numWriteCellPerOperationNeuro, numCol*activityColWrite) / 2 * numWritePulse * cell.writePulseWidth;	// LTP
					writeDynamicEnergyArray += cell.writeVoltage/2 * cell.writeVoltage/2 / cell.resMemCellAvgAtHalfVw * (numCol - MIN(numWriteCellPerOperationNeuro, numCol*activityColWrite)/2) * totalWriteTime;    // Half-selected cells on the row
					writeDynamicEnergyArray += cell.writeVoltage/2 * cell.writeVoltage/2 / cell.resMemCellAvgAtHalfVw * (numRow-1) * MIN(numWriteCellPerOperationNeuro, numCol*activityColWrite) / 2 * totalWriteTime;   // Half-selected cells on the selected columns
					// LTD
					writeDynamicEnergyArray += capCol * cell.writeVoltage * cell.writeVoltage * MIN(numWriteCellPerOperationNeuro, numCol*activityColWrite) / 2 * numWritePulse;  // Selected BLs
					writeDynamicEnergyArray += cell.writeVoltage * cell.writeVoltage / cell.resMemCellAvgAtVw * MIN(numWriteCellPerOperationNeuro, numCol*activityColWrite) / 2 * numWritePulse * cell.writePulseWidth; // LTD
					writeDynamicEnergyArray += cell.writeVoltage/2 * cell.writeVoltage/2 / cell.resMemCellAvgAtHalfVw * (numCol - MIN(numWriteCellPerOperationNeuro, numCol*activityColWrite)/2) * totalWriteTime;    // Half-selected cells on the row
					writeDynamicEnergyArray += cell.writeVoltage/2 * cell.writeVoltage/2 / cell.resMemCellAvgAtHalfVw * (numRow-1) * MIN(numWriteCellPerOperationNeuro, numCol*activityColWrite) / 2 * totalWriteTime;   // Half-selected cells on the selected columns
					// Both SET and RESET
					writeDynamicEnergyArray += capCol * cell.writeVoltage/2 * cell.writeVoltage/2 * numCol; // Unselected BLs (every BL has one time to charge to V/2 within the 2-step write)
					writeDynamicEnergyArray += capRow1 * cell.writeVoltage/2 * cell.writeVoltage/2 * (numRow-1);  // Unselected WLs
					
					writeDynamicEnergyArray *= numWriteOperationPerRow * numRow * activityRowWrite;
					
					// Read
					readDynamicEnergy += wlSwitchMatrix.readDynamicEnergy;
					readDynamicEnergy += readDynamicEnergyArray;
					readDynamicEnergy += mux.readDynamicEnergy;
					readDynamicEnergy += muxDecoder.readDynamicEnergy;
					readDynamicEnergy += readCircuit.readDynamicEnergy;
					readDynamicEnergy += subtractor.readDynamicEnergy;
					readDynamicEnergy += shiftAdd.readDynamicEnergy;

					// Write
					writeDynamicEnergy += wlSwitchMatrix.writeDynamicEnergy;
					writeDynamicEnergy += blSwitchMatrix.writeDynamicEnergy;
					writeDynamicEnergy += writeDynamicEnergyArray;
				}
				
				// Leakage
				leakage += wlDecoder.leakage;
				leakage += wlDecoderDriver.leakage;
				leakage += colDecoder.leakage;
				leakage += colDecoderDriver.leakage;
				leakage += wlSwitchMatrix.leakage;
				leakage += blSwitchMatrix.leakage;
				leakage += mux.leakage;
				leakage += muxDecoder.leakage;
				leakage += readCircuit.leakage;
				leakage += voltageSenseAmp.leakage;
				leakage += shiftAdd.leakage;
				leakage += subtractor.leakage;

			}
		}

		if (!readLatency) {
			cout << "[SubArray] Error: Need to calculate read latency first" << endl;
		} else {
			readPower = readDynamicEnergy/readLatency + leakage;
		}
		if (!writeLatency) {
			cout << "[SubArray] Error: Need to calculate write latency first" << endl;
		} else {
			writePower = writeDynamicEnergy/writeLatency + leakage;
		}

	}
}

void SubArray::PrintProperty() {
	cout << endl << endl;
	cout << "Array:" << endl;
	cout << "Area = " << heightArray*1e6 << "um x " << widthArray*1e6 << "um = " << areaArray*1e12 << "um^2" << endl;
	cout << "Read Dynamic Energy = " << readDynamicEnergyArray*1e12 << "pJ" << endl;
	cout << "Write Dynamic Energy = " << writeDynamicEnergyArray*1e12 << "pJ" << endl;
	if (cell.memCellType == Type::SRAM) {
		wlDecoder.PrintProperty("wlDecoder");
		precharger.PrintProperty("precharger");
		sramWriteDriver.PrintProperty("sramWriteDriver");
		adder.PrintProperty("adder");
		dff.PrintProperty("dff");
		subtractor.PrintProperty("subtractor");
		if (shiftAddEnable) {
			shiftAdd.PrintProperty("shiftAdd");
		}
	} else if (cell.memCellType == Type::RRAM) {
		if (cell.accessType == CMOS_access) {   // 1T1R
			if (digitalModeNeuro) {
				wlDecoder.PrintProperty("wlDecoder");
				colDecoder.PrintProperty("colDecoder");
				colDecoderDriver.PrintProperty("colDecoderDriver");
				mux.PrintProperty("mux");
				muxDecoder.PrintProperty("muxDecoder");
				voltageSenseAmp.PrintProperty("voltageSenseAmp");
				adder.PrintProperty("adder");
				dff.PrintProperty("dff");
				subtractor.PrintProperty("subtractor");
				if (shiftAddEnable) {
					shiftAdd.PrintProperty("shiftAdd");
				}
			} else {
				wlDecoderOutput.PrintProperty("wlDecoderOutput");
				wlDecoder.PrintProperty("wlDecoder");
				slSwitchMatrix.PrintProperty("slSwitchMatrix");
				blSwitchMatrix.PrintProperty("blSwitchMatrix");
				mux.PrintProperty("mux");
				muxDecoder.PrintProperty("muxDecoder");
				readCircuit.PrintProperty("readCircuit");
				subtractor.PrintProperty("subtractor");
				if (shiftAddEnable) {
					shiftAdd.PrintProperty("shiftAdd");
				}
			}
		} else {	// Crosspoint
			if (digitalModeNeuro) {
				wlDecoder.PrintProperty("wlDecoder");
				wlDecoderDriver.PrintProperty("wlDecoderDriver");
				colDecoder.PrintProperty("colDecoder");
				colDecoderDriver.PrintProperty("colDecoderDriver");
				mux.PrintProperty("mux");
				muxDecoder.PrintProperty("muxDecoder");
				voltageSenseAmp.PrintProperty("voltageSenseAmp");
				adder.PrintProperty("adder");
				dff.PrintProperty("dff");
				subtractor.PrintProperty("subtractor");
				if (shiftAddEnable) {
					shiftAdd.PrintProperty("shiftAdd");
				}
			} else {
				wlSwitchMatrix.PrintProperty("wlSwitchMatrix");
				blSwitchMatrix.PrintProperty("blSwitchMatrix");
				mux.PrintProperty("mux");
				muxDecoder.PrintProperty("muxDecoder");
				readCircuit.PrintProperty("readCircuit");
				subtractor.PrintProperty("subtractor");
				if (shiftAddEnable) {
					shiftAdd.PrintProperty("shiftAdd");
				}
			}
		}
	}
	FunctionUnit::PrintProperty("SubArray");
	cout << "Used Area = " << usedArea*1e12 << "um^2" << endl;
	cout << "Empty Area = " << emptyArea*1e12 << "um^2" << endl;
}

