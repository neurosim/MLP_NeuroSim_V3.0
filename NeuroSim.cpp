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
*   Pai-Yu Chen     Email: pchen72 at asu dot edu 
*                     
*   Xiaochen Peng   Email: xpeng15 at asu dot edu
********************************************************************************/

#include <cmath>
#include <iostream>
#include "NeuroSim.h"
#include "NeuroSim/constant.h"
#include "NeuroSim/formula.h"
#include "Cell.h"
#include "Param.h"

using namespace std;

extern Param *param;

// NeuroSim.cpp is the interface between the subarray and MLP
// it assigns the value from cell objects in Arrays to the Memcell objects in subarrays
// subarrays is used only for area calculation
void NeuroSimSubArrayInitialize(SubArray *& subArray, Array *array, InputParameter& inputParameter, Technology& tech, MemCell& cell) {

	/* Create SubArray object and link the required global objects (not initialization) */
	subArray = new SubArray(inputParameter, tech, cell);

	inputParameter.deviceRoadmap = HP;	// HP: high performance, LSTP: low power
	inputParameter.temperature = 301;	// Temperature (K)
	inputParameter.processNode = param->processNode;	// Technology node
	
	tech.Initialize(inputParameter.processNode, inputParameter.deviceRoadmap);
	
	subArray->activityRowWrite = (double)1/2;  // Dynamic parameter (to be determined)
	subArray->activityColWrite = (double)1/2;	   // Dynamic parameter (to be determined)
	subArray->activityRowRead = (double)1/2;  // Dynamic parameter (to be determined)
	subArray->spikingMode = NONSPIKING;    // NONSPIKING: input data using pulses in binary representation
										                                  
	if (subArray->spikingMode == SPIKING) {  // SPIKING: input data using # of pulses
		subArray->numReadPulse = pow(2, param->numBitInput); // use a set of pulse train to represent the input
	} else {
		subArray->numReadPulse = param->numBitInput;
	}
	subArray->numWritePulse = 8;	// Dynamic parameter (to be determined)
    
	if (DigitalNVM *temp = dynamic_cast<DigitalNVM*>(array->cell[0][0])) {
		subArray->digitalModeNeuro = 1;	// Use digital RRAM for Neuromorphic mode
	} else {
		subArray->digitalModeNeuro = 0;
	}
    
	subArray->clkFreq = param->clkFreq;		                                // Clock frequency
	subArray->numCellPerSynapse = array->numCellPerSynapse;	// # of cells per synapse
	subArray->numColMuxed = param->numColMuxed;              // How many columns share 1 read circuit (for analog RRAM) or 1 S/A (for digital RRAM)
	subArray->numWriteColMuxed = param->numWriteColMuxed;    // How many columns share 1 write column decoder driver (for digital RRAM)
	
    if (subArray->spikingMode == NONSPIKING && subArray->numReadPulse > 1) {
		subArray->shiftAddEnable = true;	// Need to shift & add the partial weighted sum
	} else {
		subArray->shiftAddEnable = false;
	}
	
    subArray->relaxArrayCellHeight = param->relaxArrayCellHeight;
	subArray->relaxArrayCellWidth = param->relaxArrayCellWidth;

	cell.heightInFeatureSize = (array->cell[0][0])->heightInFeatureSize;	// Cell height in feature size
	cell.widthInFeatureSize = (array->cell[0][0])->widthInFeatureSize;		// Cell width in feature size
   
	if (SRAM *temp = dynamic_cast<SRAM*>(array->cell[0][0])) {	// SRAM
		/* Transfer the cell properties from MLP simulator to NeuroSim */
		cell.memCellType = Type::SRAM;
		cell.widthSRAMCellNMOS = static_cast<SRAM*>(array->cell[0][0])->widthSRAMCellNMOS;
		cell.widthSRAMCellPMOS = static_cast<SRAM*>(array->cell[0][0])->widthSRAMCellPMOS;
		cell.widthAccessCMOS = static_cast<SRAM*>(array->cell[0][0])->widthAccessCMOS;
		cell.minSenseVoltage = static_cast<SRAM*>(array->cell[0][0])->minSenseVoltage;	// The minimum voltage difference for sensing
		subArray->avgWeightBit = subArray->numCellPerSynapse;   // Average weight for each synapse (value can range from 0 to numCellPerSynapse)
		subArray->numReadCellPerOperationNeuro = (int)ceil((double)array->arrayColSize / param->numColMuxed) * subArray->numCellPerSynapse;	// # of SRAM read cells in neuromorphic mode 

	}  
    else {	// eNVM (RRAM, PCM, STT-MRAM)
        cell.memCellType = Type::RRAM; // Type if RRAM if it is not 2T1F cell
		subArray->readCircuitMode  = CMOS;	// CMOS implementation for integrate-and-fire neuron
		subArray->maxNumIntBit = param->numBitPartialSum;	// Max # bits for the integrate-and-fire neuron
		
		if (subArray->digitalModeNeuro) {
			subArray->avgWeightBit = subArray->numCellPerSynapse;   // Average weight for each synapse (value can range from 0 to numCellPerSynapse)
            // code added for STT-MRAM
            if ( static_cast<DigitalNVM*>(array->cell[0][0])->parallelRead ==true) 
            {
                subArray->parallelRead =true;
            }
            else {
                 subArray->parallelRead=false;
            }
           // code added for STT-MRAM ends
           
		} else { //Analog mode
			int maxNumLevelLTP = static_cast<AnalogNVM*>(array->cell[0][0])->maxNumLevelLTP;
			int maxNumLevelLTD = static_cast<AnalogNVM*>(array->cell[0][0])->maxNumLevelLTD;
			subArray->maxNumWritePulse = (maxNumLevelLTP > maxNumLevelLTD)? maxNumLevelLTP : maxNumLevelLTD;
		}
		
		cell.accessType = (static_cast<eNVM*>(array->cell[0][0])->cmosAccess)? CMOS_access : none_access;	// CMOS_access: 1T1R (pseudo-crossbar), none_access: crossbar

		cell.resistanceOn = 1/static_cast<eNVM*>(array->cell[0][0])->avgMaxConductance;	// Ron resistance at Vr in the reported measurement data (need to recalculate below if considering the nonlinearity)
		cell.resistanceOff = 1/static_cast<eNVM*>(array->cell[0][0])->avgMinConductance;	// Roff resistance at Vr in the reported measurement dat (need to recalculate below if considering the nonlinearity)
		cell.resistanceAvg = (cell.resistanceOn + cell.resistanceOff)/2;	// Average resistance (used for energy estimation)
		cell.resCellAccess = static_cast<eNVM*>(array->cell[0][0])->resistanceAccess;   // Access transistor resistance
		cell.readVoltage = static_cast<eNVM*>(array->cell[0][0])->readVoltage;	// On-chip read voltage for memory cell
		double writeVoltageLTP = static_cast<eNVM*>(array->cell[0][0])->writeVoltageLTP;
		double writeVoltageLTD = static_cast<eNVM*>(array->cell[0][0])->writeVoltageLTD;
		cell.writeVoltage = sqrt(writeVoltageLTP * writeVoltageLTP + writeVoltageLTD * writeVoltageLTD);	// Use an average value of write voltage for NeuroSim
		cell.readPulseWidth = static_cast<eNVM*>(array->cell[0][0])->readPulseWidth;
		double writePulseWidthLTP = static_cast<eNVM*>(array->cell[0][0])->writePulseWidthLTP;
		double writePulseWidthLTD = static_cast<eNVM*>(array->cell[0][0])->writePulseWidthLTD;
		cell.writePulseWidth = (writePulseWidthLTP + writePulseWidthLTD) / 2;
		cell.nonlinearIV = static_cast<eNVM*>(array->cell[0][0])->nonlinearIV; // This option is to consider I-V nonlinearity in cross-point array or not
		cell.nonlinearity = (cell.nonlinearIV)? 10 : 2;	// This is the nonlinearity for the current ratio at Vw and Vw/2
		if (cell.nonlinearIV) {
			double Vr_exp = 1;  // XXX: Modify this to Vr in the reported measurement data (can be different than cell.readVoltage)
			// Calculation of resistance at on-chip Vr
			cell.resistanceOn = NonlinearResistance(cell.resistanceOn, cell.nonlinearity, cell.writeVoltage, Vr_exp, cell.readVoltage);
			cell.resistanceOff = NonlinearResistance(cell.resistanceOff, cell.nonlinearity, cell.writeVoltage, Vr_exp, cell.readVoltage);
			cell.resistanceAvg = (cell.resistanceOn + cell.resistanceOff)/2;    // Average resistance (for energy estimation)
		}
		    cell.accessVoltage = 1.1;	// Gate voltage for the transistor in 1T1R
	}


	int numRow = array->arrayRowSize;	// Transfer the parameter of # of rows from the MLP simulator to NeuroSim
	int numCol = array->arrayColSize * subArray->numCellPerSynapse;	// Transfer the parameter of # of columns from the MLP simulator to NeuroSim (times the # of cells per synapse for SRAM)
    
    // code added for STT-MRAM
    if (  (subArray->digitalModeNeuro && subArray->parallelRead == true)) {
         numCol += 1;  // need reference column for parallelRead
                                 // 1 more column for the 3T1C cell
     }
    // code added for STT-MRAM ends
	if (param->numColMuxed > numCol) {	// Set the upperbound of param->numColMuxed
		param->numColMuxed = numCol;
	}

	cell.featureSize = array->wireWidth * 1e-9;
	if(cell.featureSize <= 0) {
		puts("NeuroSim does not take ideal array. Has Assigned the width to 200nm.");
		cell.featureSize = 200e-9;
	}

	subArray->numWriteCellPerOperationNeuro = (int)ceil((double)numCol / subArray->numWriteColMuxed);
	
	/* NeuroSim SubArray Initialization */
	double unitLengthWireResistance = array->unitLengthWireResistance;
	subArray->Initialize(numRow, numCol, unitLengthWireResistance);
	/* Recalculate wire resistance after possible layout adjustment by NeuroSim */
	array->wireResistanceRow = subArray->lengthRow / numCol * unitLengthWireResistance; // the wire resistance of each cell along row direction
	array->wireResistanceCol = subArray->lengthCol / numRow * unitLengthWireResistance; // the wire resistance of each cell along column direction
	/* Transfer the wire capacitances from NeuroSim to MLP simulator */
	array->wireCapRow = subArray->capRow1;
	array->wireCapCol = subArray->capCol;
	array->wireGateCapRow = subArray->capRow2;
	array->wireCapBLCol = subArray->lengthCol * 0.2e-15/1e-6;	// For BL cap of digital eNVM in 1T1R
	/* Transfer the write energy of SRAM cell from NeuroSim to MLP simulator */
	array->writeEnergySRAMCell = subArray->cell.capSRAMCell * subArray->tech.vdd * subArray->tech.vdd * 2;	// flip Q and Q_bar
}

void NeuroSimSubArrayArea(SubArray *subArray) { // calculate the area from Subarray class
	subArray->CalculateArea();
}

double NeuroSimSubArrayReadLatency(SubArray *subArray) {	// For 1 weighted sum task on selected columns
	if (!param->NeuroSimDynamicPerformance) { return 0; }	// Skip this function if param->NeuroSimDynamicPerformance is false
	if (subArray->cell.memCellType == Type::SRAM) {   // SRAM
		subArray->wlDecoder.CalculateLatency(1e20, subArray->capRow1, NULL, subArray->numRow * subArray->numReadPulse * subArray->activityRowRead, 1);	// Don't care write
		subArray->precharger.CalculateLatency(1e20, subArray->capCol, subArray->numRow * subArray->numReadPulse * subArray->activityRowRead, 1);	// Don't care write
		subArray->senseAmp.CalculateLatency(subArray->numRow * subArray->numReadPulse * subArray->activityRowRead);
		subArray->adder.CalculateLatency(1e20, subArray->dff.capTgDrain, subArray->numRow * subArray->numReadPulse * subArray->activityRowRead);
		subArray->dff.CalculateLatency(1e20, subArray->numRow * subArray->numReadPulse * subArray->activityRowRead);
		subArray->subtractor.CalculateLatency(1e20, subArray->dff.capTgDrain, subArray->numReadPulse);
		if (subArray->shiftAddEnable) {
			subArray->shiftAdd.CalculateLatency(subArray->numReadPulse);	// There are numReadPulse times of shift-and-add
		}
		double resPullDown = CalculateOnResistance(subArray->cell.widthSRAMCellNMOS * subArray->tech.featureSize, NMOS, subArray->inputParameter.temperature, subArray->tech);
		double tau = (subArray->resCellAccess + resPullDown) * (subArray->capCellAccess + subArray->capCol) + subArray->resCol * subArray->capCol / 2;
		tau *= log(subArray->tech.vdd / (subArray->tech.vdd - subArray->cell.minSenseVoltage / 2));   /* one signal raises and the other drops, so cell.minSenseVoltage/2 is enough */
		double gm = CalculateTransconductance(subArray->cell.widthAccessCMOS * subArray->tech.featureSize, NMOS, subArray->tech);
		double beta = 1 / (resPullDown * gm);
		double colRamp = 0;
		subArray->colDelay = horowitz(tau, beta, subArray->wlDecoder.rampOutput, &colRamp) * subArray->numRow * subArray->numReadPulse * subArray->activityRowRead;

		return 	subArray->wlDecoder.readLatency +
				subArray->precharger.readLatency +
				subArray->colDelay +
				subArray->senseAmp.readLatency +
				subArray->adder.readLatency +
				subArray->dff.readLatency +
				subArray->subtractor.readLatency +
				subArray->shiftAdd.readLatency;

	} 
    else {	// eNVM
		if (subArray->digitalModeNeuro) {	// Digital eNVM, row by row operation
			if (subArray->cell.accessType == CMOS_access) {   // 1T1R
				if(subArray->parallelRead == true) 
                {   // for the parallel readout
                    //void NewSwitchMatrix::CalculateLatency(double _rampInput, double _capLoad, double _resLoad, double numRead, double numWrite) {	// For simplicity, assume shift register is ideal
                    double capToDrive=MAX(subArray->capRow2,subArray->capRow1);
                    double resToDrive=subArray->resRow;
                    subArray->wlBlSwitchMatrix.CalculateLatency(1e20, capToDrive, resToDrive, subArray->numReadPulse * subArray->activityRowRead, 1);
                   
                    // only need the wl-bl decoder
                    double capBL = subArray->lengthCol * 0.2e-15 / 1e-6;
                    subArray->colDelay = 2.3 * subArray->resCol * capBL; //column delay

 
                    // the read circuit
                    // The input capacitance of the read circuit
                    double Cin_ReadCircuit = subArray->capCol + subArray->mux.capTgDrain * (2 + subArray->numColMuxed - 1) + subArray->readCircuit.capTgDrain + subArray->readCircuit.capPmosGate;

                    // Use the maximum read current to determine the readpulsewidth
                    double Imax = subArray->numRow * subArray->cell.readVoltage / subArray->cell.resMemCellOn;
                    subArray->cell.readPulseWidth = Cin_ReadCircuit * subArray->readCircuit.voltageIntThreshold / Imax * subArray->readCircuit.maxNumIntPerCycle;

                    // Delay at the Mux the mux is driving the read circuit
                    double colRamp=0;
                    subArray->mux.CalculateLatency(colRamp, Cin_ReadCircuit, 1); // the drive resistance should be the input resistance of the read circuit, the cap is the cap of

                    // Here numColMuxed can mean how many synapses share 1 adder or how many columns share 1 S/A
                    int numAdder = (int)ceil(((double)subArray->numCol / subArray->numCellPerSynapse) / subArray->numColMuxed);   // numCol is divisible by numCellPerSynapse
                    int numInput = numAdder * subArray->numCellPerSynapse; // number of input of the mux
                    subArray->muxDecoder.CalculateLatency(1e20, subArray->mux.capTgGateN * numInput, subArray->mux.capTgGateP * numInput, 1, 1);
                    subArray->readCircuit.CalculateLatency(subArray->numReadPulse);
                    subArray->subtractor.CalculateLatency(1e20, 0, subArray->numReadPulse);
                   if (subArray->shiftAddEnable) {
                       // two shift adders are needed. one to add
                       subArray->shiftAdd.CalculateLatency(subArray->numReadPulse);
                   }
                                
                   return  MAX(subArray->wlBlSwitchMatrix.readLatency, subArray->muxDecoder.readLatency + subArray->mux.readLatency)+
                           subArray->readCircuit.readLatency +
                           subArray->subtractor.readLatency +
                           subArray->colDelay+ // need furthercheck
                           subArray->shiftAdd.readLatency;                
                  }
          else {
                   double capBL = subArray->lengthCol * 0.2e-15 / 1e-6;
				   subArray->wlDecoder.CalculateLatency(1e20, subArray->capRow2, NULL, subArray->numRow * subArray->numReadPulse * subArray->activityRowRead, 1);	// Don't care write
				   double colRamp = 0;
				   double tau = subArray->resCol * capBL / 2 * (subArray->cell.resMemCellOff + subArray->resCol / 3) / (subArray->cell.resMemCellOff + subArray->resCol);
				   subArray->colDelay = horowitz(tau, 0, 1e20, &colRamp);
				   subArray->colDelay = 2.3 * subArray->resCol * capBL;
				   subArray->mux.CalculateLatency(colRamp, 0, 1);
				   // Here numColMuxed can mean how many synapses share 1 adder or how many columns share 1 S/A
				   int numAdder = (int)ceil(((double)subArray->numCol / subArray->numCellPerSynapse) / subArray->numColMuxed);   // numCol is divisible by numCellPerSynapse
				   int numInput = numAdder * subArray->numCellPerSynapse; // number of input of the mux
				   subArray->muxDecoder.CalculateLatency(1e20, subArray->mux.capTgGateN * numInput, subArray->mux.capTgGateP * numInput, 1, 1);
				   double capInputLoad = capBL + subArray->mux.capTgDrain * (2 + subArray->numColMuxed - 1);
				   subArray->voltageSenseAmp.CalculateLatency(capInputLoad, subArray->numRow * subArray->numReadPulse * subArray->activityRowRead);
				   subArray->adder.CalculateLatency(1e20, subArray->dff.capTgDrain, subArray->numRow * subArray->numReadPulse * subArray->activityRowRead);
				   subArray->dff.CalculateLatency(1e20, subArray->numRow * subArray->numReadPulse * subArray->activityRowRead);
				   subArray->subtractor.CalculateLatency(1e20, subArray->dff.capTgDrain, subArray->numReadPulse);
				   if (subArray->shiftAddEnable) {
					   subArray->shiftAdd.CalculateLatency(subArray->numReadPulse);    // There are numReadPulse times of shift-and-add
				}

				return  MAX(subArray->wlDecoder.readLatency, subArray->muxDecoder.readLatency + subArray->mux.readLatency)+
						subArray->voltageSenseAmp.readLatency +
						subArray->adder.readLatency +
						subArray->dff.readLatency +
						subArray->subtractor.readLatency +
                        // subArray->colDelay+
						subArray->shiftAdd.readLatency;
                }
                        
			} else {        // Cross-point
				double wlDecoderLoad = subArray->colDecoderDriver.capInvInput + subArray->colDecoderDriver.capTgGateN + subArray->colDecoderDriver.capTgGateP;
				subArray->wlDecoder.CalculateLatency(1e20, wlDecoderLoad, NULL, subArray->numRow * subArray->activityRowRead * subArray->numReadPulse, 1);	// Don't care write
				subArray->wlDecoderDriver.CalculateLatency(subArray->wlDecoder.rampOutput, subArray->capRow1, subArray->capRow1, subArray->resRow, subArray->numRow * subArray->activityRowRead * subArray->numReadPulse, 1);	// Don't care write
				double colRamp = 0;
				double tau = subArray->resCol * subArray->capCol / 2 * (subArray->cell.resMemCellOff + subArray->resCol / 3) / (subArray->cell.resMemCellOff + subArray->resCol);
				subArray->colDelay = horowitz(tau, 0, 1e20, &colRamp);
				subArray->colDelay = 2.3 * subArray->resCol * subArray->capCol;
				subArray->mux.CalculateLatency(colRamp, 0, 1);
				// Here numColMuxed can mean how many synapses share 1 adder or how many columns share 1 S/A
				int numAdder = (int)ceil(((double)subArray->numCol / subArray->numCellPerSynapse) / subArray->numColMuxed);   // numCol is divisible by numCellPerSynapse
				int numInput = numAdder * subArray->numCellPerSynapse;
				subArray->muxDecoder.CalculateLatency(1e20, subArray->mux.capTgGateN * numInput, subArray->mux.capTgGateP*numInput, 1, 1);
				double capInputLoad = subArray->capCol + subArray->mux.capTgDrain * (2 + subArray->numColMuxed - 1);
				subArray->voltageSenseAmp.CalculateLatency(capInputLoad, subArray->numRow * subArray->numReadPulse * subArray->activityRowRead);
				subArray->adder.CalculateLatency(1e20, subArray->dff.capTgDrain, subArray->numRow * subArray->numReadPulse * subArray->activityRowRead);
				subArray->dff.CalculateLatency(1e20, subArray->numRow * subArray->numReadPulse * subArray->activityRowRead);
				subArray->subtractor.CalculateLatency(1e20, subArray->dff.capTgDrain, subArray->numReadPulse);
				if (subArray->shiftAddEnable) {
					subArray->shiftAdd.CalculateLatency(subArray->numReadPulse);    // There are numReadPulse times of shift-and-add
				}
				return  MAX(subArray->wlDecoder.readLatency + subArray->wlDecoderDriver.readLatency, subArray->muxDecoder.readLatency + subArray->mux.readLatency);
						subArray->voltageSenseAmp.readLatency +
						subArray->adder.readLatency +
						subArray->dff.readLatency +
						subArray->subtractor.readLatency +
						subArray->shiftAdd.readLatency;
			}
		} else {	// Analog eNVM
			if (subArray->cell.accessType == CMOS_access) {   // 1T1R
				subArray->wlDecoder.CalculateLatency(1e20, subArray->wlDecoderOutput.capNorInput, NULL, 1, 1);	// Don't care write
				subArray->wlDecoderOutput.CalculateLatency(subArray->wlDecoder.rampOutput, subArray->capRow2, subArray->resRow, 1, 1);	// Don't care write
				subArray->blSwitchMatrix.CalculateLatency(1e20, subArray->capRow1, subArray->resRow, subArray->numReadPulse, 1);    // Don't care write
				if (subArray->readCircuit.mode == CMOS) {
                    // Cin is the capacitance to collect the charge
					double Cin = subArray->capCol + subArray->mux.capTgDrain * (2 + subArray->numColMuxed - 1) + subArray->readCircuit.capTgDrain + subArray->readCircuit.capPmosGate;
					// the maximum read current
                    double Imax = subArray->numRow * subArray->cell.readVoltage / subArray->cell.resMemCellOn;
					subArray->cell.readPulseWidth = Cin * subArray->readCircuit.voltageIntThreshold / Imax * subArray->readCircuit.maxNumIntPerCycle;
				} else {    // mode==OSCILLATION
					double Cin = subArray->capCol + subArray->mux.capTgDrain * (2 + subArray->numColMuxed - 1) + subArray->readCircuit.capInvInput;
					double Rmin = subArray->cell.resMemCellOn / subArray->numRow;
					double Rp = 1 / (1/Rmin + 1/subArray->readCircuit.R_OSC_OFF);
					double t_rise = -Rp * Cin * log((subArray->readCircuit.Vth - subArray->readCircuit.Vrow * Rp / Rmin) / (subArray->readCircuit.Vhold - subArray->readCircuit.Vrow * Rp / Rmin));
					subArray->cell.readPulseWidth = t_rise * subArray->readCircuit.maxNumIntPerCycle;
				}
				subArray->readCircuit.CalculateLatency(subArray->numReadPulse);
				subArray->subtractor.CalculateLatency(1e20, 0, subArray->numReadPulse);
				if (subArray->shiftAddEnable) {
					subArray->shiftAdd.CalculateLatency(subArray->numReadPulse);
				}
				return 	subArray->wlDecoderOutput.readLatency +
						subArray->blSwitchMatrix.readLatency +
						subArray->readCircuit.readLatency +
						subArray->subtractor.readLatency +
						subArray->shiftAdd.readLatency;

			} else {		// Cross-point
				subArray->wlSwitchMatrix.CalculateLatency(1e20, subArray->capRow1, subArray->resRow, subArray->numReadPulse, 1);	// Don't care write
				if (subArray->readCircuit.mode == CMOS) {
					double Cin = subArray->capCol + subArray->mux.capTgDrain * (2 + subArray->numColMuxed - 1) + subArray->readCircuit.capTgDrain + subArray->readCircuit.capPmosGate;
					double Imax = subArray->numRow * subArray->cell.readVoltage / subArray->cell.resMemCellOn;
					subArray->cell.readPulseWidth = Cin * subArray->readCircuit.voltageIntThreshold / Imax * subArray->readCircuit.maxNumIntPerCycle;
				} else {    // mode==OSCILLATION
					double Cin = subArray->capCol + subArray->mux.capTgDrain * (2 + subArray->numColMuxed - 1) + subArray->readCircuit.capInvInput;
					double Rmin = subArray->cell.resMemCellOn / subArray->numRow;
					double Rp = 1 / (1/Rmin + 1/subArray->readCircuit.R_OSC_OFF);
					double t_rise = -Rp * Cin * log((subArray->readCircuit.Vth - subArray->readCircuit.Vrow * Rp / Rmin) / (subArray->readCircuit.Vhold - subArray->readCircuit.Vrow * Rp / Rmin));
					subArray->cell.readPulseWidth = t_rise * subArray->readCircuit.maxNumIntPerCycle;
				}
				subArray->readCircuit.CalculateLatency(subArray->numReadPulse);
				subArray->subtractor.CalculateLatency(1e20, 0, subArray->numReadPulse);
				if (subArray->shiftAddEnable) {
					subArray->shiftAdd.CalculateLatency(subArray->numReadPulse);
				}
				return 	subArray->wlSwitchMatrix.readLatency +
						subArray->readCircuit.readLatency +
						subArray->subtractor.readLatency +
						subArray->shiftAdd.readLatency;
			}
		}
	}
}

double NeuroSimSubArrayWriteLatency(SubArray *subArray, int numWriteOperationPerRow, double sumWriteLatencyAnalogNVM) {	// For 1 weight update task of whole array
	if (!param->NeuroSimDynamicPerformance) { return 0; }	// Skip this function if param->NeuroSimDynamicPerformance is false
	subArray->activityRowWrite = 1;
	subArray->activityColWrite = 1;
	if (subArray->cell.memCellType == Type::SRAM) {	// SRAM
		subArray->wlDecoder.CalculateLatency(1e20, subArray->capRow1, NULL, 1, subArray->numRow * subArray->activityRowWrite);	// Don't care read
		subArray->precharger.CalculateLatency(1e20, subArray->capCol, 1, numWriteOperationPerRow * subArray->numRow * subArray->activityRowWrite);	// Don't care read
		subArray->sramWriteDriver.CalculateLatency(1e20, subArray->capCol, subArray->resCol, numWriteOperationPerRow * subArray->numRow * subArray->activityRowWrite);	// Don't care read
		
		// Write (assume the average delay of pullup and pulldown inverter in SRAM cell)
		double resPull = (CalculateOnResistance(subArray->cell.widthSRAMCellNMOS * subArray->tech.featureSize, NMOS, subArray->inputParameter.temperature, subArray->tech) + CalculateOnResistance(subArray->cell.widthSRAMCellPMOS * subArray->tech.featureSize, PMOS, subArray->inputParameter.temperature, subArray->tech)) / 2;    // take average
		double tau = resPull * subArray->cell.capSRAMCell;
		double gm = (CalculateTransconductance(subArray->cell.widthSRAMCellNMOS * subArray->tech.featureSize, NMOS, subArray->tech) + CalculateTransconductance(subArray->cell.widthSRAMCellPMOS * subArray->tech.featureSize, PMOS, subArray->tech)) / 2;   // take average
		double beta = 1 / (resPull * gm);

		return	subArray->wlDecoder.writeLatency +
				subArray->precharger.writeLatency +
				subArray->sramWriteDriver.writeLatency +
				horowitz(tau, beta, 1e20, NULL) * numWriteOperationPerRow * subArray->numRow * subArray->activityRowWrite;

	}
    else {	// eNVM
		if (subArray->digitalModeNeuro) {   // Digital eNVM
			if (subArray->cell.accessType == CMOS_access) {   // 1T1R
				  if(subArray->parallelRead==true) {// parallel read for digital neuron
                     double capToDrive=MAX(subArray->capRow2,subArray->capRow1);
                     double resToDrive=subArray->resRow;
                     subArray->wlBlSwitchMatrix.CalculateLatency(1e20, capToDrive, resToDrive, 1, subArray->numRow * subArray->activityRowWrite);
                     return subArray->wlBlSwitchMatrix.writeLatency;
                  } 
                  else {
                       double capBL = subArray->lengthCol * 0.2e-15 / 1e-6;
                       subArray->wlDecoder.CalculateLatency(1e20, subArray->capRow2, NULL, 1, subArray->numRow * subArray->activityRowWrite);	// Don't care read
                       double colDecoderLoad = (subArray->colDecoderDriver.capInvInput + subArray->colDecoderDriver.capTgGateN * 2 + subArray->colDecoderDriver.capTgGateP) * subArray->numWriteCellPerOperationNeuro;
                       subArray->colDecoder.CalculateLatency(1e20, colDecoderLoad, NULL, 1, subArray->numRow * subArray->activityRowWrite * numWriteOperationPerRow);	// Doesn't matter for read
                       if (subArray->colDecoder.rampOutput > 0) {
                                  subArray->colDecoderDriver.CalculateLatency(subArray->colDecoder.rampOutput, subArray->capCol, capBL, subArray->resCol, 1, subArray->numRow * subArray->activityRowWrite * numWriteOperationPerRow * 2);	// Doesn't matter for read. *2 means 2-step write
                        } 
                        else {    // The case where column decoder doesn't exist
                                  subArray->colDecoderDriver.CalculateLatency(1e20, subArray->capCol, capBL, subArray->resCol, 1, subArray->numRow * subArray->activityRowWrite * numWriteOperationPerRow * 2);	// Doesn't matter for read. *2 means 2-step write
                        }
                       return 	MAX(subArray->wlDecoder.writeLatency, subArray->colDecoder.writeLatency + subArray->colDecoderDriver.writeLatency);
			      }
            } 
            else {	// Cross-point
				double wlDecoderLoad = subArray->colDecoderDriver.capInvInput + subArray->colDecoderDriver.capTgGateN + subArray->colDecoderDriver.capTgGateP;
				subArray->wlDecoder.CalculateLatency(1e20, wlDecoderLoad, NULL, 1, subArray->numRow * subArray->activityRowWrite);
				subArray->wlDecoderDriver.CalculateLatency(subArray->wlDecoder.rampOutput, subArray->capRow1, subArray->capRow1, subArray->resRow, 1, subArray->numRow * subArray->activityRowWrite);
				double colDecoderLoad = (subArray->colDecoderDriver.capInvInput + subArray->colDecoderDriver.capTgGateN + subArray->colDecoderDriver.capTgGateP) * subArray->numWriteCellPerOperationNeuro;
				subArray->colDecoder.CalculateLatency(1e20, colDecoderLoad, NULL, 1, subArray->numRow * subArray->activityRowWrite * numWriteOperationPerRow);  // Doesn't matter for read
				if (subArray->colDecoder.rampOutput > 0) {
					subArray->colDecoderDriver.CalculateLatency(subArray->colDecoder.rampOutput, subArray->capCol, subArray->capCol, subArray->resCol, 1, subArray->numRow * subArray->activityRowWrite * numWriteOperationPerRow * 2); // Doesn't matter for read. *2 means 2-step write
				} else {    // The case where column decoder doesn't exist
					subArray->colDecoderDriver.CalculateLatency(1e20, subArray->capCol, subArray->capCol, subArray->resCol, 1, subArray->numRow * subArray->activityRowWrite * numWriteOperationPerRow * 2);  // Doesn't matter for read. *2 means 2-step write
				}
				return	MAX(subArray->wlDecoder.writeLatency + subArray->wlDecoderDriver.writeLatency, subArray->colDecoder.writeLatency + subArray->colDecoderDriver.writeLatency);
			}
		} else {	// Analog eNVM
			if (subArray->cell.accessType == CMOS_access) {	// 1T1R
				subArray->wlDecoder.CalculateLatency(1e20, subArray->wlDecoderOutput.capNorInput, NULL, 1, subArray->numRow * subArray->activityRowWrite);	// Don't care read
				subArray->wlDecoderOutput.CalculateLatency(subArray->wlDecoder.rampOutput, subArray->capRow2, subArray->resRow, 1, subArray->numRow * subArray->activityRowWrite);	// Don't care read
				subArray->blSwitchMatrix.writeLatency = sumWriteLatencyAnalogNVM;

				return 	subArray->wlDecoder.writeLatency +
						subArray->wlDecoderOutput.writeLatency +
						subArray->blSwitchMatrix.writeLatency;

			} else {	// Cross-point
				subArray->wlSwitchMatrix.writeLatency = sumWriteLatencyAnalogNVM;

				return subArray->wlSwitchMatrix.writeLatency;
			}
		}
	}
}

double NeuroSimSubArrayReadEnergy(SubArray *subArray) {	// For 1 weighted sum task on selected columns
	if (!param->NeuroSimDynamicPerformance) { return 0; }	// Skip this function if param->NeuroSimDynamicPerformance is false
	if (subArray->cell.memCellType == Type::SRAM) {   // SRAM
		subArray->wlDecoder.CalculatePower(subArray->numRow * subArray->numReadPulse * subArray->activityRowRead, 1);	// Don't care write
		subArray->precharger.CalculatePower(subArray->numRow * subArray->numReadPulse * subArray->activityRowRead, 1);	// Don't care write
		subArray->senseAmp.CalculatePower(subArray->numRow * subArray->numReadPulse * subArray->activityRowRead);
		subArray->adder.CalculatePower(subArray->numRow * subArray->numReadPulse * subArray->activityRowRead, subArray->numReadCellPerOperationNeuro / subArray->numCellPerSynapse);
		subArray->dff.CalculatePower(subArray->numRow * subArray->numReadPulse * subArray->activityRowRead, subArray->numReadCellPerOperationNeuro / subArray->numCellPerSynapse * (subArray->adder.numBit+1));
		subArray->subtractor.CalculatePower(subArray->numReadPulse, subArray->numReadCellPerOperationNeuro / subArray->numCellPerSynapse);
		if (subArray->shiftAddEnable) {
			subArray->shiftAdd.CalculatePower(subArray->numReadPulse);
		}

		return 	subArray->wlDecoder.readDynamicEnergy +
				subArray->precharger.readDynamicEnergy +
				subArray->senseAmp.readDynamicEnergy +
				subArray->adder.readDynamicEnergy +
				subArray->dff.readDynamicEnergy +
				subArray->subtractor.readDynamicEnergy +
				subArray->shiftAdd.readDynamicEnergy;

	} 
    else {    // eNVM
		if (subArray->digitalModeNeuro) {   // Digital eNVM
			if (subArray->cell.accessType == CMOS_access) {   // 1T1R
                   if (subArray->parallelRead == true) 
                   {      
                     double numReadCells = (int)ceil((double)subArray->numCol / subArray->numColMuxed);
                     subArray->wlBlSwitchMatrix.CalculatePower(subArray->activityRowRead * subArray->numReadPulse, 1);	// Don't care write
                     subArray->slSwitchMatrix.CalculatePower(subArray->numReadPulse, 1);	// Don't care write
                     subArray->mux.CalculatePower(subArray->activityRowRead * subArray->numReadPulse);	// Mux still consumes energy during row-by-row read
                     subArray->muxDecoder.CalculatePower(subArray->numReadPulse, 1);	// Don't care write
                     subArray->readCircuit.CalculatePower(subArray->numReadPulse);
                     // Need to chack the second parameter of subtractorpower
                     subArray->subtractor.CalculatePower(subArray->numReadPulse, subArray->numReadCellPerOperationNeuro / subArray->numCellPerSynapse);
                     if (subArray->shiftAddEnable) 
                        {
                                subArray->shiftAdd.CalculatePower(subArray->numReadPulse);
                         }
                     return  subArray->wlBlSwitchMatrix.readDynamicEnergy+
                                 subArray->slSwitchMatrix.readDynamicEnergy+
                                 subArray->mux.readDynamicEnergy +
                                 subArray->muxDecoder.readDynamicEnergy +
                                 subArray->readCircuit.readDynamicEnergy+
                                 subArray->subtractor.readDynamicEnergy +
                                 subArray->shiftAdd.readDynamicEnergy;
                    } 
                    else 
                    {      // row-by-row readout
                         double numReadCells = (int)ceil((double)subArray->numCol / subArray->numColMuxed);
                         subArray->wlDecoder.CalculatePower(subArray->numRow * subArray->activityRowRead * subArray->numReadPulse, 1);	// Don't care write
                         subArray->mux.CalculatePower(subArray->numRow * subArray->activityRowRead * subArray->numReadPulse);	// Mux still consumes energy during row-by-row read
                         subArray->muxDecoder.CalculatePower(subArray->numReadPulse, 1);	// Don't care write
                         subArray->voltageSenseAmp.CalculatePower(subArray->numRow * subArray->activityRowRead * subArray->numReadPulse);
                         subArray->adder.CalculatePower(subArray->numRow * subArray->numReadPulse * subArray->activityRowRead, numReadCells / subArray->numCellPerSynapse);
                         subArray->dff.CalculatePower(subArray->numRow * subArray->numReadPulse * subArray->activityRowRead, numReadCells / subArray->numCellPerSynapse*(subArray->adder.numBit+1));	// +1 because the adder output is 1 bit more than the input
                         subArray->subtractor.CalculatePower(subArray->numReadPulse, subArray->numReadCellPerOperationNeuro / subArray->numCellPerSynapse);
                         if (subArray->shiftAddEnable) 
                         {
                                  subArray->shiftAdd.CalculatePower(subArray->numReadPulse);
                          }
                         return  subArray->wlDecoder.readDynamicEnergy +
                                     subArray->mux.readDynamicEnergy +
                                     subArray->muxDecoder.readDynamicEnergy +
                                     subArray->voltageSenseAmp.readDynamicEnergy +
                                     subArray->adder.readDynamicEnergy +
                                     subArray->dff.readDynamicEnergy +
                                     subArray->subtractor.readDynamicEnergy +
                                     subArray->shiftAdd.readDynamicEnergy;
                     }
			} else {    // Cross-point
				double numReadCells = (int)ceil((double)subArray->numCol / subArray->numColMuxed);
				subArray->wlDecoder.CalculatePower(subArray->numRow * subArray->activityRowRead * subArray->numReadPulse, 1);	// Don't care write
				subArray->wlDecoderDriver.CalculatePower(numReadCells, 1, subArray->numRow * subArray->activityRowRead * subArray->numReadPulse, 1);
				subArray->mux.CalculatePower(subArray->numRow * subArray->activityRowRead * subArray->numReadPulse);	// Mux still consumes energy during row-by-row read
				subArray->muxDecoder.CalculatePower(subArray->numReadPulse, 1);
				subArray->voltageSenseAmp.CalculatePower(subArray->numRow * subArray->activityRowRead * subArray->numReadPulse);
				subArray->adder.CalculatePower(subArray->numRow * subArray->numReadPulse * subArray->activityRowRead, numReadCells / subArray->numCellPerSynapse);
				subArray->dff.CalculatePower(subArray->numRow * subArray->numReadPulse * subArray->activityRowRead, numReadCells / subArray->numCellPerSynapse*(subArray->adder.numBit+1));	// +1 because the adder output is 1 bit more than the input
				subArray->subtractor.CalculatePower(subArray->numReadPulse, subArray->numReadCellPerOperationNeuro / subArray->numCellPerSynapse);
				if (subArray->shiftAddEnable) {
					subArray->shiftAdd.CalculatePower(subArray->numReadPulse);
				}
				return	subArray->wlDecoder.readDynamicEnergy +
						subArray->wlDecoderDriver.readDynamicEnergy +
						subArray->mux.readDynamicEnergy +
						subArray->muxDecoder.readDynamicEnergy +
						subArray->voltageSenseAmp.readDynamicEnergy +
						subArray->adder.readDynamicEnergy +
						subArray->dff.readDynamicEnergy +
						subArray->subtractor.readDynamicEnergy +
						subArray->shiftAdd.readDynamicEnergy;
			}
		} else {	// Analog eNVM
			if (subArray->cell.accessType == CMOS_access) {   // 1T1R
				subArray->wlDecoder.CalculatePower(subArray->numReadPulse, 1);	// Don't care write
				subArray->wlDecoderOutput.CalculatePower(subArray->numReadPulse, 1);	// Don't care write
				subArray->blSwitchMatrix.CalculatePower(subArray->numReadPulse, 1);	// Don't care write
				subArray->mux.CalculatePower(subArray->numReadPulse);
				subArray->muxDecoder.CalculatePower(subArray->numReadPulse, 1);	// Don't care write
				subArray->readCircuit.CalculatePower(subArray->numReadPulse);
				subArray->subtractor.CalculatePower(subArray->numReadPulse, subArray->numReadCellPerOperationNeuro / subArray->numCellPerSynapse);
				if (subArray->shiftAddEnable) {
					subArray->shiftAdd.CalculatePower(subArray->numReadPulse);
				}
				return	subArray->wlDecoder.readDynamicEnergy +
						    subArray->wlDecoderOutput.readDynamicEnergy +
						    subArray->blSwitchMatrix.readDynamicEnergy +
						    subArray->mux.readDynamicEnergy +
						    subArray->muxDecoder.readDynamicEnergy +
						    subArray->readCircuit.readDynamicEnergy +
						    subArray->subtractor.readDynamicEnergy +
						    subArray->shiftAdd.readDynamicEnergy;
			} else {	// Cross-point
				subArray->wlSwitchMatrix.CalculatePower(subArray->numReadPulse, 1);	// Don't care write
				subArray->mux.CalculatePower(subArray->numReadPulse);
				subArray->muxDecoder.CalculatePower(subArray->numReadPulse, 1);	// Don't care write
				subArray->readCircuit.CalculatePower(subArray->numReadPulse);
				subArray->subtractor.CalculatePower(subArray->numReadPulse, subArray->numReadCellPerOperationNeuro / subArray->numCellPerSynapse);
				if (subArray->shiftAddEnable) {
					subArray->shiftAdd.CalculatePower(subArray->numReadPulse);
				}
				return 	subArray->wlSwitchMatrix.readDynamicEnergy +
						subArray->mux.readDynamicEnergy +
						subArray->muxDecoder.readDynamicEnergy +
						subArray->readCircuit.readDynamicEnergy +
						subArray->subtractor.readDynamicEnergy +
						subArray->shiftAdd.readDynamicEnergy;
			}
		}
	}
}

double NeuroSimSubArrayWriteEnergy(SubArray *subArray, int numWriteOperationPerRow, double numWriteCellPerOperation) {	// For 1 weight update task of one row
	if (!param->NeuroSimDynamicPerformance) { return 0; }	// Skip this function if param->NeuroSimDynamicPerformance is false
	subArray->activityRowWrite = 1;
	subArray->activityColWrite = 1;
	if (subArray->cell.memCellType == Type::SRAM) {   // SRAM
		subArray->wlDecoder.CalculatePower(1, 1);	// Don't care read, should be different for parallel read
		subArray->precharger.CalculatePower(1, numWriteOperationPerRow);	// Don't care read
		subArray->sramWriteDriver.CalculatePower(numWriteOperationPerRow);

		return	subArray->wlDecoder.writeDynamicEnergy +
				subArray->precharger.writeDynamicEnergy +
				subArray->sramWriteDriver.writeDynamicEnergy;

	} 
    else {    // eNVM
		if (subArray->digitalModeNeuro) {   // Digital eNVM
			if (subArray->cell.accessType == CMOS_access) {   // 1T1R
                 if(subArray->parallelRead==true)
                 {  //parallel readout
                      subArray->wlBlSwitchMatrix.CalculatePower(1,1);
                      subArray->slSwitchMatrix.CalculatePower(1, numWriteOperationPerRow);	// Don't care read
                                    
                      return subArray->wlBlSwitchMatrix.writeDynamicEnergy+
                             subArray->slSwitchMatrix.writeDynamicEnergy;
            
                  }
                  else
                  {
                      subArray->wlDecoder.CalculatePower(1, 1);	// Don't care read
                      subArray->colDecoder.CalculatePower(1, numWriteOperationPerRow);  // Doesn't matter for read
                      subArray->colDecoderDriver.CalculatePower(1, numWriteCellPerOperation, 1, numWriteOperationPerRow * 2); // Doesn't matter for read. *2 means 2-step write
                      return	subArray->wlDecoder.writeDynamicEnergy +
                                subArray->colDecoder.writeDynamicEnergy +
                                subArray->colDecoderDriver.writeDynamicEnergy;
                  }
			} else {    // Cross-point
				subArray->wlDecoder.CalculatePower(1, 1);	// Don't care read
				subArray->wlDecoderDriver.CalculatePower(1, numWriteCellPerOperation, 1, numWriteOperationPerRow);
				subArray->colDecoder.CalculatePower(1, numWriteOperationPerRow);  // Doesn't matter for read
				subArray->colDecoderDriver.CalculatePower(1, numWriteCellPerOperation, 1, numWriteOperationPerRow * 2); // Doesn't matter for read. *2 means 2-step write
				return	subArray->wlDecoder.writeDynamicEnergy +
						subArray->wlDecoderDriver.writeDynamicEnergy +
						subArray->colDecoder.writeDynamicEnergy +
						subArray->colDecoderDriver.writeDynamicEnergy;
			}
		} else {    // Analog eNVM
			if (subArray->cell.accessType == CMOS_access) {   // 1T1R
				subArray->blSwitchMatrix.numWritePulse = 1;	// Does not matter
				subArray->slSwitchMatrix.numWritePulse = subArray->numWritePulse;
				subArray->wlDecoder.CalculatePower(1, 1);	// Don't care read
				subArray->wlDecoderOutput.CalculatePower(1, 1);	// Don't care read
				subArray->blSwitchMatrix.CalculatePower(1, 1);	// Don't care read
				subArray->slSwitchMatrix.CalculatePower(1, numWriteOperationPerRow);	// Don't care read
				return 	subArray->wlDecoder.writeDynamicEnergy +
						subArray->wlDecoderOutput.writeDynamicEnergy +
						subArray->blSwitchMatrix.writeDynamicEnergy +
						subArray->slSwitchMatrix.writeDynamicEnergy;
			} else {    // Cross-point
				subArray->wlSwitchMatrix.numWritePulse = subArray->numWritePulse;
				subArray->blSwitchMatrix.numWritePulse = subArray->numWritePulse;
				subArray->wlSwitchMatrix.CalculatePower(1, 1);	// Don't care read
				subArray->blSwitchMatrix.CalculatePower(1, numWriteOperationPerRow);	// Don't care read
				return  subArray->wlSwitchMatrix.writeDynamicEnergy +
						subArray->blSwitchMatrix.writeDynamicEnergy;
			}
		}
	}
}

double NeuroSimSubArrayLeakagePower(SubArray *subArray) {
    //printf("calculating the array leakage power");
	if (subArray->cell.memCellType == Type::SRAM) {	// SRAM
		subArray->wlDecoder.CalculatePower(1, 1);
		subArray->precharger.CalculatePower(1, 1);
		subArray->sramWriteDriver.CalculatePower(1);
		subArray->senseAmp.CalculatePower(1);
		subArray->adder.CalculatePower(1, 1);
		subArray->dff.CalculatePower(1, 1);
		subArray->subtractor.CalculatePower(1, 1);
		if (subArray->shiftAddEnable) {
			subArray->shiftAdd.CalculatePower(1);
		}
		
		// Array leakage (assume 2 INV)
		subArray->leakage += CalculateGateLeakage(INV, 1, subArray->cell.widthSRAMCellNMOS * subArray->tech.featureSize, subArray->cell.widthSRAMCellPMOS * subArray->tech.featureSize, subArray->inputParameter.temperature, subArray->tech) * subArray->tech.vdd * 2;
		subArray->leakage *= subArray->numRow * subArray->numCol;

		subArray->leakage += subArray->wlDecoder.leakage;
		subArray->leakage += subArray->precharger.leakage;
		subArray->leakage += subArray->sramWriteDriver.leakage;
		subArray->leakage += subArray->senseAmp.leakage;
		subArray->leakage += subArray->adder.leakage;
		subArray->leakage += subArray->dff.leakage;
		subArray->leakage += subArray->subtractor.leakage;
		subArray->leakage += subArray->shiftAdd.leakage;

	} 
    else if(subArray->cell.memCellType==Type::_2T1F){
             // leakage from the 2 access transistor??
            subArray->wlSwitchMatrix.CalculatePower(1, 1);	// Don't care numRead and numWrite
            subArray->blSwitchMatrix.CalculatePower(1, 1);	// Don't care numRead and numWrite
            subArray->slSwitchMatrix.CalculatePower(1, 1);	// Don't care numRead and numWrite
            subArray->plSwitchMatrix.CalculatePower(1, 1);	// Don't care numRead and numWrite
            subArray->mux.CalculatePower(1);	// Don't care numRead
            subArray->muxDecoder.CalculatePower(1, 1);	// Don't care numRead and numWrite
            subArray->readCircuit.CalculatePower(1);	// Don't care numRead
            subArray->subtractor.CalculatePower(1, 1);
            if (subArray->shiftAddEnable) {
                subArray->shiftAdd.CalculatePower(1);	// Don't care numRead
            }
		    // width of the PMOS is 0;
            // leakage power of the two access transistor of the 3T1C cell
            // should use the access transistor's width
            subArray->leakage += CalculateGateLeakage(INV, 1, subArray->cell.widthAccessNMOS * subArray->tech.featureSize, subArray->cell.widthAccessPMOS * subArray->tech.featureSize , subArray->inputParameter.temperature, subArray->tech) * subArray->tech.vdd * 2;
		    subArray->leakage *= subArray->numRow * subArray->numCol;            
            subArray->leakage += subArray->wlSwitchMatrix.leakage;
            subArray->leakage += subArray->blSwitchMatrix.leakage;
            subArray->leakage += subArray->slSwitchMatrix.leakage;
            subArray->leakage += subArray->plSwitchMatrix.leakage;
            subArray->leakage += subArray->mux.leakage;
            subArray->leakage += subArray->muxDecoder.leakage;
            subArray->leakage += subArray->readCircuit.leakage;
            subArray->leakage += subArray->subtractor.leakage;
            subArray->leakage += subArray->shiftAdd.leakage;       
    }
    else {	// eNVM
		if (subArray->digitalModeNeuro) {   // Digital eNVM
			if (subArray->cell.accessType == CMOS_access) {   // 1T1R
                if (subArray->parallelRead==true)
                {
				    /*subArray->wlDecoder.CalculatePower(1, 1);
                    subArray->wlDecoderOutput.CalculatePower(1,1);
                    subArray->blSwitchMatrix.CalculatePower(1, 1);	// Don't care numRead and numWrite */
                    
                    subArray->wlBlSwitchMatrix.CalculatePower(1, 1);
                    subArray->slSwitchMatrix.CalculatePower(1, 1);
				    subArray->mux.CalculatePower(1);
				    subArray->muxDecoder.CalculatePower(1, 1);
                    subArray->readCircuit.CalculatePower(1);
				    subArray->subtractor.CalculatePower(1, 1);
				    if (subArray->shiftAddEnable) {
					    subArray->shiftAdd.CalculatePower(1);
				     }
				    /*
                    subArray->leakage += subArray->wlDecoder.leakage;
				    subArray->leakage += subArray->wlDecoderOutput.leakage;
                    subArray->leakage += subArray->blSwitchMatrix.leakage; 
                     */
                    subArray->leakage += subArray->wlBlSwitchMatrix.leakage;
                    subArray->leakage += subArray->mux.leakage;
				    subArray->leakage += subArray->muxDecoder.leakage;
                    
                    subArray->leakage += subArray->slSwitchMatrix.leakage;
                    subArray->leakage += subArray->readCircuit.leakage;
				    subArray->leakage += subArray->subtractor.leakage;
				    subArray->leakage += subArray->shiftAdd.leakage;
                }
                else {
                    subArray->wlDecoder.CalculatePower(1, 1);
				    subArray->colDecoder.CalculatePower(1, 1);
				    subArray->colDecoderDriver.CalculatePower(1, 1, 1, 1);
				    subArray->mux.CalculatePower(1);
				    subArray->muxDecoder.CalculatePower(1, 1);
				    subArray->voltageSenseAmp.CalculatePower(1);
				    subArray->adder.CalculatePower(1, 1);
				    subArray->dff.CalculatePower(1, 1);
				    subArray->subtractor.CalculatePower(1, 1);
				    if (subArray->shiftAddEnable) {
					    subArray->shiftAdd.CalculatePower(1);
				    }
				    subArray->leakage += subArray->wlDecoder.leakage;
				    subArray->leakage += subArray->colDecoder.leakage;
				    subArray->leakage += subArray->colDecoderDriver.leakage;
				    subArray->leakage += subArray->mux.leakage;
				    subArray->leakage += subArray->muxDecoder.leakage;
				    subArray->leakage += subArray->voltageSenseAmp.leakage;
				    subArray->leakage += subArray->adder.leakage;
				    subArray->leakage += subArray->dff.leakage;
				    subArray->leakage += subArray->subtractor.leakage;
				    subArray->leakage += subArray->shiftAdd.leakage;
                }
			} else {    // Cross-point
				subArray->wlDecoder.CalculatePower(1, 1);
				subArray->wlDecoderDriver.CalculatePower(1, 1, 1, 1);
				subArray->colDecoder.CalculatePower(1, 1);
				subArray->colDecoderDriver.CalculatePower(1, 1, 1, 1);
				subArray->mux.CalculatePower(1);
				subArray->muxDecoder.CalculatePower(1, 1);
				subArray->voltageSenseAmp.CalculatePower(1);
				subArray->adder.CalculatePower(1, 1);
				subArray->dff.CalculatePower(1, 1);
				subArray->subtractor.CalculatePower(1, 1);
				if (subArray->shiftAddEnable) {
					subArray->shiftAdd.CalculatePower(1);
				}
				subArray->leakage += subArray->wlDecoder.leakage;
				subArray->leakage += subArray->wlDecoderDriver.leakage;
				subArray->leakage += subArray->colDecoder.leakage;
				subArray->leakage += subArray->colDecoderDriver.leakage;
				subArray->leakage += subArray->mux.leakage;
				subArray->leakage += subArray->muxDecoder.leakage;
				subArray->leakage += subArray->voltageSenseAmp.leakage;
				subArray->leakage += subArray->adder.leakage;
				subArray->leakage += subArray->dff.leakage;
				subArray->leakage += subArray->subtractor.leakage;
				subArray->leakage += subArray->shiftAdd.leakage;
			}
		} else {    // Analog eNVM
			if (subArray->cell.accessType == CMOS_access) {	// 1T1R
				subArray->wlDecoder.CalculatePower(1, 1);	// Don't care numRead and numWrite
				subArray->wlDecoderOutput.CalculatePower(1, 1);	// Don't care numRead and numWrite
				subArray->blSwitchMatrix.CalculatePower(1, 1);	// Don't care numRead and numWrite
				subArray->slSwitchMatrix.CalculatePower(1, 1);	// Don't care numRead and numWrite
				subArray->mux.CalculatePower(1);	// Don't care numRead
				subArray->muxDecoder.CalculatePower(1, 1);	// Don't care numRead and numWrite
				subArray->readCircuit.CalculatePower(1);	// Don't care numRead
				subArray->subtractor.CalculatePower(1, 1);
				if (subArray->shiftAddEnable) {
					subArray->shiftAdd.CalculatePower(1);	// Don't care numRead
				}
				subArray->leakage += subArray->wlDecoder.leakage;
				subArray->leakage += subArray->wlDecoderOutput.leakage;
				subArray->leakage += subArray->blSwitchMatrix.leakage;
				subArray->leakage += subArray->slSwitchMatrix.leakage;
				subArray->leakage += subArray->mux.leakage;
				subArray->leakage += subArray->muxDecoder.leakage;
				subArray->leakage += subArray->readCircuit.leakage;
				subArray->leakage += subArray->subtractor.leakage;
				subArray->leakage += subArray->shiftAdd.leakage;
			} else {	// Cross-point
				subArray->wlSwitchMatrix.CalculatePower(1, 1);	// Don't care numRead and numWrite
				subArray->blSwitchMatrix.CalculatePower(1, 1);	// Don't care numRead and numWrite
				subArray->mux.CalculatePower(1);	// Don't care numRead
				subArray->muxDecoder.CalculatePower(1, 1);	// Don't care numRead and numWrite
				subArray->readCircuit.CalculatePower(1);	// Don't care numRead
				subArray->subtractor.CalculatePower(1, 1);
				if (subArray->shiftAddEnable) {
					subArray->shiftAdd.CalculatePower(1);	// Don't care numRead
				}
				subArray->leakage += subArray->wlSwitchMatrix.leakage;
				subArray->leakage += subArray->blSwitchMatrix.leakage;
				subArray->leakage += subArray->mux.leakage;
				subArray->leakage += subArray->muxDecoder.leakage;
				subArray->leakage += subArray->readCircuit.leakage;
				subArray->leakage += subArray->subtractor.leakage;
				subArray->leakage += subArray->shiftAdd.leakage;
			}
		}
	}
}

// Neuron refers to the periphery circuit
void NeuroSimNeuronInitialize(SubArray *& subArray, InputParameter& inputParameter, Technology& tech, MemCell& cell, Adder& adder, Mux& mux, RowDecoder& muxDecoder, DFF& dff, Subtractor& subtractor) {
    int numAdderBit; 
	if (subArray->shiftAddEnable) 
     {	// Here we only support adder in non-spiking fashion
    // numReadPulse: controls the precise of the input vectors
		 numAdderBit = subArray->shiftAdd.numAdderBit + 1 + subArray->shiftAdd.numReadPulse - 1;
	  } 
    else 
    { //No need to use shift adder
		   if (cell.memCellType == Type::SRAM) 
            {	// SRAM
			     numAdderBit = subArray->adder.numBit + 1;
		    } 
            else 
            {	// eNVM
                 numAdderBit = param->numBitPartialSum;
		    }
	}
	numAdderBit = numAdderBit + 2;	// Need 1 more bit for *2 in 2w'-1 of MLP algorithm, and 1 more bit for 2's complement implementation
	
	int numAdder = (int)ceil((double)subArray->numCol/subArray->numCellPerSynapse/subArray->numColMuxed);

	// Only need the MSB of adder output to determine it is positive or negative in 2's complement
	dff.Initialize(subArray->numCol/subArray->numCellPerSynapse, subArray->clkFreq);
	if (subArray->numColMuxed > 1) {
		mux.Initialize(numAdder*param->numBitInput, subArray->numColMuxed, NULL, true);    // Digital Mux
		muxDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(subArray->numColMuxed)), true);
		adder.Initialize(numAdderBit, numAdder);
		subtractor.Initialize(numAdderBit, numAdder);
	} else {	// No need for Mux and Mux decoder
		adder.Initialize(numAdderBit, numAdder);
		subtractor.Initialize(numAdderBit, numAdder);
	}
}

void NeuroSimNeuronArea(SubArray *subArray, Adder& adder, Mux& mux, RowDecoder& muxDecoder, DFF& dff, Subtractor& subtractor, double *height, double *width) {
	adder.CalculateArea(NULL, subArray->widthArray, NONE);
	if (subArray->numColMuxed > 1) {
		mux.CalculateArea(NULL, subArray->widthArray, NONE);	// Digital Mux
		muxDecoder.CalculateArea(mux.height * 4, NULL, NONE);	// Set muxDecoder height to be 4 times (for example) mux height to avoid large muxDecoder width
	}
	dff.CalculateArea(NULL, subArray->widthArray, NONE);
	subtractor.CalculateArea(NULL, subArray->widthArray, NONE);
		
	*height = MAX(adder.height + mux.height + dff.height + subtractor.height, muxDecoder.height);
	*width = subArray->widthArray + muxDecoder.width;
}

double NeuroSimNeuronReadLatency(SubArray *subArray, Adder& adder, Mux& mux, RowDecoder& muxDecoder, DFF& dff, Subtractor& subtractor) {	// For 1 weighted sum task on selected columns
	if (!param->NeuroSimDynamicPerformance) { return 0; }	// Skip this function if param->NeuroSimDynamicPerformance is false
	if (subArray->numColMuxed > 1) 
    {
		adder.CalculateLatency(1e20, mux.capTgDrain, 1);
		mux.CalculateLatency(adder.rampOutput, dff.capTgDrain, 1);
		muxDecoder.CalculateLatency(1e20, mux.capTgGateN * adder.numAdder, mux.capTgGateP * adder.numAdder, 1, 1);	// Don't care write
		subtractor.CalculateLatency(1e20, mux.capTgDrain, 1);
	} 
    else 
    {	// No need for Mux and Mux decoder
		adder.CalculateLatency(1e20, dff.capTgDrain, 1);
		subtractor.CalculateLatency(1e20, dff.capTgDrain, 1);
	}
	dff.CalculateLatency(1e20, 1);
    if (subArray-> parallelRead== true)
    {
       return mux.readLatency + subtractor.readLatency;     
    }
    else
    {
	    return adder.readLatency + mux.readLatency + dff.readLatency + subtractor.readLatency;
    }
}

double NeuroSimNeuronReadEnergy(SubArray *subArray, Adder& adder, Mux& mux, RowDecoder& muxDecoder, DFF& dff, Subtractor& subtractor) {	// For 1 weighted sum task on selected columns
	if (!param->NeuroSimDynamicPerformance) { return 0; }	// Skip this function if param->NeuroSimDynamicPerformance is false
	adder.CalculatePower(1, adder.numAdder);
	if (subArray->numColMuxed > 1) 
    {
		mux.CalculatePower(1);
		muxDecoder.CalculatePower(1, 1);	// Don't care write
	}
	dff.CalculatePower(1, adder.numAdder);
	subtractor.CalculatePower(1, adder.numAdder);
    if(subArray->parallelRead==true)
    {
        return mux.readDynamicEnergy + muxDecoder.readDynamicEnergy+subtractor.readDynamicEnergy;
    }
    else 
    {
        return adder.readDynamicEnergy + mux.readDynamicEnergy + muxDecoder.readDynamicEnergy + dff.readDynamicEnergy + subtractor.readDynamicEnergy;
    }
    
}

double NeuroSimNeuronLeakagePower(SubArray *subArray, Adder& adder, Mux& mux, RowDecoder& muxDecoder, DFF& dff, Subtractor& subtractor) { // Same as NeuroSimNeuronReadEnergy
    adder.CalculatePower(1, adder.numAdder);
	if (subArray->numColMuxed > 1) 
    {
		mux.CalculatePower(1);
		muxDecoder.CalculatePower(1, 1);	// Don't care write
	}
	dff.CalculatePower(1, adder.numAdder);
	subtractor.CalculatePower(1, adder.numAdder);
    if(subArray->parallelRead==true)
    {
	    return  mux.leakage + muxDecoder.leakage  + subtractor.leakage;
    }
    else
    {
        return adder.leakage + mux.leakage + muxDecoder.leakage + dff.leakage + subtractor.leakage;

    }
}

double NeuroSimNeuronTransferEnergy(SubArray *subArray, int numWriteOperationPerRow, double numWriteCellPerOperation)
{
        double energyReadLSB=0, energyWriteMSB=0;
        subArray->wlSwitchMatrix_LSB.CalculatePower(subArray->numReadPulse, 1);	// Don't care write
        subArray->blSwitchMatrix_LSB.CalculatePower(subArray->numReadPulse, 1);	// Don't care write
        subArray->mux.CalculatePower(subArray->numReadPulse);
        subArray->mux_2to1.CalculatePower(subArray->numReadPulse);
        subArray->muxDecoder.CalculatePower(subArray->numReadPulse, 1);	// Don't care write
        subArray->muxDecoder_2to1.CalculatePower(subArray->numReadPulse,1);
        subArray->readCircuit.CalculatePower(subArray->numReadPulse);
        subArray->subtractor.CalculatePower(subArray->numReadPulse, subArray->numReadCellPerOperationNeuro / subArray->numCellPerSynapse);
        if (subArray->shiftAddEnable) {
            subArray->shiftAdd.CalculatePower(subArray->numReadPulse);
        }
        energyReadLSB = subArray->wlSwitchMatrix_LSB.readDynamicEnergy +
        subArray->blSwitchMatrix_LSB.readDynamicEnergy +
        subArray->mux.readDynamicEnergy +
        subArray->numColMuxed*subArray->mux_2to1.readDynamicEnergy+
        subArray->muxDecoder.readDynamicEnergy +
        subArray->numColMuxed*subArray->muxDecoder_2to1.readDynamicEnergy+
        subArray->readCircuit.readDynamicEnergy +
        subArray->subtractor.readDynamicEnergy +
        subArray->shiftAdd.readDynamicEnergy; 
        subArray->blSwitchMatrix.numWritePulse = 1;	// Does not matter
        subArray->slSwitchMatrix.numWritePulse = subArray->numWritePulse;
        subArray->wlSwitchMatrix.CalculatePower(1, 1);	// Don't care read
        subArray->blSwitchMatrix.CalculatePower(1, 1);	// Don't care read
        subArray->slSwitchMatrix.CalculatePower(1, numWriteOperationPerRow);	// Don't care read
        energyWriteMSB = subArray->wlSwitchMatrix.writeDynamicEnergy +
                                      subArray->blSwitchMatrix.writeDynamicEnergy +
                                      subArray->slSwitchMatrix.writeDynamicEnergy; 
        return energyReadLSB + energyWriteMSB; 
}


