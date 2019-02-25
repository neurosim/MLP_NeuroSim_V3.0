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

#include <ctime>
#include <iostream>
#include <math.h>
#include "formula.h"
#include "Array.h"
#include "Cell.h"


/* General eNVM */
void AnalogNVM::WriteEnergyCalculation(double wireCapCol) {
    //printf("calculating write energy consumption\n");
	if (nonlinearIV) {  // Currently only for cross-point array
		/* I-V nonlinearity */
		double conductancePrevAtVwLTP = NonlinearConductance(conductancePrev, NL, writeVoltageLTP, readVoltage, writeVoltageLTP);
		double conductancePrevAtHalfVwLTP = NonlinearConductance(conductancePrev, NL, writeVoltageLTP, readVoltage, writeVoltageLTP/2);
		double conductancePrevAtVwLTD = NonlinearConductance(conductancePrev, NL, writeVoltageLTD, readVoltage, writeVoltageLTD);
		double conductancePrevAtHalfVwLTD = NonlinearConductance(conductancePrev, NL, writeVoltageLTD, readVoltage, writeVoltageLTD/2);
		conductanceAtVwLTP = NonlinearConductance(conductance, NL, writeVoltageLTP, readVoltage, writeVoltageLTP);
		conductanceAtHalfVwLTP = NonlinearConductance(conductance, NL, writeVoltageLTP, readVoltage, writeVoltageLTP/2);
		conductanceAtVwLTD = NonlinearConductance(conductance, NL, writeVoltageLTD, readVoltage, writeVoltageLTD);
		conductanceAtHalfVwLTD = NonlinearConductance(conductance, NL, writeVoltageLTD, readVoltage, writeVoltageLTD/2);
		if (numPulse > 0) { // If the cell needs LTP pulses
			writeEnergy = writeVoltageLTP * writeVoltageLTP * (conductancePrevAtVwLTP+conductanceAtVwLTP)/2 * writePulseWidthLTP * numPulse;
			writeEnergy += writeVoltageLTP * writeVoltageLTP * wireCapCol * numPulse;
			if (nonIdenticalPulse) {
				writeVoltageLTD = VinitLTD + (VinitLTD + VstepLTD * maxNumLevelLTD);
			}
			writeEnergy += writeVoltageLTD/2 * writeVoltageLTD/2 * conductanceAtHalfVwLTD * writeLatencyLTD;    // Half-selected during LTD phase (use the new conductance value if LTP phase is before LTD phase)
			writeEnergy += writeVoltageLTD/2 * writeVoltageLTD/2 * wireCapCol;
		} else if (numPulse < 0) {  // If the cell needs LTD pulses
			if (nonIdenticalPulse) {
				writeVoltageLTP = VinitLTP + (VinitLTP + VstepLTP * maxNumLevelLTP);
			}
			writeEnergy = writeVoltageLTP/2 * writeVoltageLTP/2 * conductancePrevAtHalfVwLTP * writeLatencyLTP;    // Half-selected during LTP phase (use the old conductance value if LTP phase is before LTD phase)
			writeEnergy += writeVoltageLTP/2 * writeVoltageLTP/2 * wireCapCol;
			writeEnergy += writeVoltageLTD * writeVoltageLTD * wireCapCol * (-numPulse);
			writeEnergy += writeVoltageLTD * writeVoltageLTD * (conductancePrevAtVwLTD+conductanceAtVwLTD)/2 * writePulseWidthLTD * (-numPulse);
		} else {    // Half-selected during both LTP and LTD phases
			if (nonIdenticalPulse) {
				writeVoltageLTP = VinitLTP + (VinitLTP + VstepLTP * maxNumLevelLTP);
				writeVoltageLTD = VinitLTD + (VinitLTD + VstepLTD * maxNumLevelLTD);
			}
			writeEnergy = writeVoltageLTP/2 * writeVoltageLTP/2 * conductancePrevAtHalfVwLTP * writeLatencyLTP;
			writeEnergy += writeVoltageLTP/2 * writeVoltageLTP/2 * wireCapCol;
			writeEnergy += writeVoltageLTD/2 * writeVoltageLTD/2 * conductancePrevAtHalfVwLTD * writeLatencyLTD;
			writeEnergy += writeVoltageLTD/2 * writeVoltageLTD/2 * wireCapCol;
		}
	} else {    // If not cross-point array or not considering I-V nonlinearity
		if (FeFET) {	// FeFET structure
			if (cmosAccess) {
				if (numPulse > 0) { // If the cell needs LTP pulses
					writeEnergy = writeVoltageLTP * writeVoltageLTP * (gateCapFeFET + wireCapCol) * numPulse;
					if (nonIdenticalPulse) {
						writeVoltageLTD = VinitLTD + VstepLTD * maxNumLevelLTD;
					}
					writeEnergy += writeVoltageLTD * writeVoltageLTD * (gateCapFeFET + wireCapCol);
				} else if (numPulse < 0) {  // If the cell needs LTD pulses
					writeEnergy = writeVoltageLTD * writeVoltageLTD * (gateCapFeFET + wireCapCol) * (-numPulse);
				} else {    // Half-selected during both LTP and LTD phases
					if (nonIdenticalPulse) {
						writeVoltageLTD = VinitLTD + VstepLTD * maxNumLevelLTD;
					}
					writeEnergy = writeVoltageLTD * writeVoltageLTD * (gateCapFeFET + wireCapCol);
				}
			} else {
				puts("FeFET structure is not compatible with crossbar");
				exit(-1);
			}
		} else {
			if (numPulse > 0) { // If the cell needs LTP pulses
				writeEnergy = writeVoltageLTP * writeVoltageLTP * (conductancePrev+conductance)/2 * writePulseWidthLTP * numPulse;
				writeEnergy += writeVoltageLTP * writeVoltageLTP * wireCapCol * numPulse;
				if (!cmosAccess) {	// Crossbar
					if (nonIdenticalPulse) {
						writeVoltageLTD = VinitLTD + (VinitLTD + VstepLTD * maxNumLevelLTD);
					}
					writeEnergy += writeVoltageLTD/2 * writeVoltageLTD/2 * conductance * writeLatencyLTD;    // Half-selected during LTD phase (use the new conductance value if LTP phase is before LTD phase)
					writeEnergy += writeVoltageLTD/2 * writeVoltageLTD/2 * wireCapCol;
				}
			} else if (numPulse < 0) {  // If the cell needs LTD pulses
				if (!cmosAccess) {	// Crossbar
					if (nonIdenticalPulse) {
						writeVoltageLTP = VinitLTP + (VinitLTP + VstepLTP * maxNumLevelLTP);
					}
					writeEnergy = writeVoltageLTP/2 * writeVoltageLTP/2 * conductancePrev * writeLatencyLTP;    // Half-selected during LTP phase (use the old conductance value if LTP phase is before LTD phase)
					writeEnergy += writeVoltageLTP/2 * writeVoltageLTP/2 * wireCapCol;
				} else {	// 1T1R
					if (nonIdenticalPulse) {
						writeVoltageLTP = VinitLTP + VstepLTP * maxNumLevelLTP;
					}
					writeEnergy = writeVoltageLTP * writeVoltageLTP * wireCapCol;
				}
				writeEnergy += writeVoltageLTD * writeVoltageLTD * wireCapCol * (-numPulse);
				writeEnergy += writeVoltageLTD * writeVoltageLTD * (conductancePrev+conductance)/2 * writePulseWidthLTD * (-numPulse);
			} else {    // Half-selected during both LTP and LTD phases
				if (!cmosAccess) {	// Crossbar
					if (nonIdenticalPulse) {
						writeVoltageLTP = VinitLTP + (VinitLTP + VstepLTP * maxNumLevelLTP);
						writeVoltageLTD = VinitLTD + (VinitLTD + VstepLTD * maxNumLevelLTD);
					}
					writeEnergy = writeVoltageLTP/2 * writeVoltageLTP/2 * conductancePrev * writeLatencyLTP;
					writeEnergy += writeVoltageLTP/2 * writeVoltageLTP/2 * wireCapCol;
					writeEnergy += writeVoltageLTD/2 * writeVoltageLTD/2 * conductancePrev * writeLatencyLTD;
					writeEnergy += writeVoltageLTD/2 * writeVoltageLTD/2 * wireCapCol;
				} else {	// 1T1R
					if (nonIdenticalPulse) {
						writeVoltageLTP = VinitLTP + VstepLTP * maxNumLevelLTP;
					}
					writeEnergy = writeVoltageLTP * writeVoltageLTP * wireCapCol;
				}
			}
		}
	}
    //printf("writeEnergy is %.4e\n", writeEnergy);
}

/* Ideal device (no weight update nonlinearity) */
IdealDevice::IdealDevice(int x, int y) {
	this->x = x; this->y = y;	// Cell location: x (column) and y (row) start from index 0
	maxConductance = 5e-6;		// Maximum cell conductance (S)
	minConductance = 100e-9;	    // Minimum cell conductance (S)
	avgMaxConductance = maxConductance; // Average maximum cell conductance (S)
	avgMinConductance = minConductance; // Average minimum cell conductance (S)
	conductance = minConductance;	// Current conductance (S) (dynamic variable)
	conductancePrev = conductance;	// Previous conductance (S) (dynamic variable)
	readVoltage = 0.5;	// On-chip read voltage (Vr) (V)
	readPulseWidth = 5e-9;	// Read pulse width (s) (will be determined by ADC)
	writeVoltageLTP = 2;	// Write voltage (V) for LTP or weight increase
	writeVoltageLTD = 2;	// Write voltage (V) for LTD or weight decrease
	writePulseWidthLTP = 10e-9;	// Write pulse width (s) for LTP or weight increase
	writePulseWidthLTD = 10e-9;	// Write pulse width (s) for LTD or weight decrease
	writeEnergy = 0;	// Dynamic variable for calculation of write energy (J)
	maxNumLevelLTP = 63;	// Maximum number of conductance states during LTP or weight increase
	maxNumLevelLTD = 63;	// Maximum number of conductance states during LTD or weight decrease
	numPulse = 0;	// Number of write pulses used in the most recent write operation (dynamic variable)
	cmosAccess = true;	// True: Pseudo-crossbar (1T1R), false: cross-point
	FeFET = false;		// True: FeFET structure (Pseudo-crossbar only, should be cmosAccess=1)
	gateCapFeFET = 2.1717e-18;	// Gate capacitance of FeFET (F)
	resistanceAccess = 15e3;	// The resistance of transistor (Ohm) in Pseudo-crossbar array when turned ON
	nonlinearIV = false;	// Consider I-V nonlinearity or not (Currently for cross-point array only)
	nonIdenticalPulse = false;	// Use non-identical pulse scheme in weight update or not (should be false here)
								// Don't care other non-identical pulse parameters
	NL = 10;	// Nonlinearity in write scheme (the current ratio between Vw and Vw/2), assuming for the LTP side
	if (nonlinearIV) {	// Currently for cross-point array only
		double Vr_exp = readVoltage;  // XXX: Modify this value to Vr in the reported measurement data (can be different than readVoltage)
		// Calculation of conductance at on-chip Vr
		maxConductance = NonlinearConductance(maxConductance, NL, writeVoltageLTP, Vr_exp, readVoltage);
		minConductance = NonlinearConductance(minConductance, NL, writeVoltageLTP, Vr_exp, readVoltage);
	}
	readNoise = false;	// Consider read noise or not
	sigmaReadNoise = 0.25;	// Sigma of read noise in gaussian distribution
	gaussian_dist = new std::normal_distribution<double>(0, sigmaReadNoise);	// Set up mean and stddev for read noise
	
	/* Conductance range variation */	
	conductanceRangeVar = false;	// Consider variation of conductance range or not
	maxConductanceVar = 0;	// Sigma of maxConductance variation (S)
	minConductanceVar = 0;	// Sigma of minConductance variation (S)
	std::mt19937 localGen;
	localGen.seed(std::time(0));
	gaussian_dist_maxConductance = new std::normal_distribution<double>(0, maxConductanceVar);
	gaussian_dist_minConductance = new std::normal_distribution<double>(0, minConductanceVar);
	if (conductanceRangeVar) {
		maxConductance += (*gaussian_dist_maxConductance)(localGen);
		minConductance += (*gaussian_dist_minConductance)(localGen);
		if (minConductance >= maxConductance || maxConductance < 0 || minConductance < 0 ) {	// Conductance variation check
			puts("[Error] Conductance variation check not passed. The variation may be too large.");
			exit(-1);
		}
		// Use the code below instead for re-choosing the variation if the check is not passed
		//do {
		//	maxConductance = avgMaxConductance + (*gaussian_dist_maxConductance)(localGen);
		//	minConductance = avgMinConductance + (*gaussian_dist_minConductance)(localGen);
		//} while (minConductance >= maxConductance || maxConductance < 0 || minConductance < 0);
	}
	
	heightInFeatureSize = cmosAccess? 4 : 2;	// Cell height = 4F (Pseudo-crossbar) or 2F (cross-point)
	widthInFeatureSize = cmosAccess? (FeFET? 6 : 4) : 2;	// Cell width = 6F (FeFET) or 4F (Pseudo-crossbar) or 2F (cross-point)
}

double IdealDevice::Read(double voltage) {
	extern std::mt19937 gen;
	// TODO: nonlinear read
	if (readNoise) {
		return voltage * conductance * (1 + (*gaussian_dist)(gen));
	} else {
		return voltage * conductance;
	}
}

void IdealDevice::Write(double deltaWeightNormalized, double weight, double minWeight, double maxWeight) {
	extern std::mt19937 gen;
	if (deltaWeightNormalized >= 0) {
		deltaWeightNormalized = deltaWeightNormalized/(maxWeight-minWeight);
		deltaWeightNormalized = truncate(deltaWeightNormalized, maxNumLevelLTP);
		numPulse = deltaWeightNormalized * maxNumLevelLTP;
	} else {
		deltaWeightNormalized = deltaWeightNormalized/(maxWeight-minWeight);
		deltaWeightNormalized = truncate(deltaWeightNormalized, maxNumLevelLTD);
		numPulse = deltaWeightNormalized * maxNumLevelLTD;	                          // will be a negative number
	}
	double conductanceNew = conductance + deltaWeightNormalized * (maxConductance - minConductance);
	if (conductanceNew > maxConductance) {
		conductanceNew = maxConductance;
	} else if (conductanceNew < minConductance) {
		conductanceNew = minConductance;
	}

	/* Write latency calculation */
	if (numPulse > 0) {	// LTP
		writeLatencyLTP = numPulse * writePulseWidthLTP;
		writeLatencyLTD = 0;
	} else {	// LTD
		writeLatencyLTP = 0;
		writeLatencyLTD = -numPulse * writePulseWidthLTD;
	}
	
	conductancePrev = conductance;
	conductance = conductanceNew;
}

/* Real Device */
RealDevice::RealDevice(int x, int y) {
	this->x = x; this->y = y;	// Cell location: x (column) and y (row) start from index 0
	maxConductance = 3.8462e-8;		// Maximum cell conductance (S)
	minConductance = 3.0769e-9;	// Minimum cell conductance (S)
	avgMaxConductance = maxConductance; // Average maximum cell conductance (S)
	avgMinConductance = minConductance; // Average minimum cell conductance (S)
	conductance = minConductance;	// Current conductance (S) (dynamic variable)
	conductancePrev = conductance;	// Previous conductance (S) (dynamic variable)
	readVoltage = 0.5;	// On-chip read voltage (Vr) (V)
	readPulseWidth = 5e-9;	// Read pulse width (s) (will be determined by ADC)
	writeVoltageLTP = 3.2;	// Write voltage (V) for LTP or weight increase
	writeVoltageLTD = 2.8;	// Write voltage (V) for LTD or weight decrease
	writePulseWidthLTP = 300e-6;	// Write pulse width (s) for LTP or weight increase
	writePulseWidthLTD = 300e-6;	// Write pulse width (s) for LTD or weight decrease
	writeEnergy = 0;	// Dynamic variable for calculation of write energy (J)
	maxNumLevelLTP = 97;	// Maximum number of conductance states during LTP or weight increase
	maxNumLevelLTD = 100;	// Maximum number of conductance states during LTD or weight decrease
	numPulse = 0;	// Number of write pulses used in the most recent write operation (dynamic variable)
	cmosAccess = true;	// True: Pseudo-crossbar (1T1R), false: cross-point
    FeFET = false;		// True: FeFET structure (Pseudo-crossbar only, should be cmosAccess=1)
	gateCapFeFET = 2.1717e-18;	// Gate capacitance of FeFET (F)
	resistanceAccess = 15e3;	// The resistance of transistor (Ohm) in Pseudo-crossbar array when turned ON
	nonlinearIV = false;	// Consider I-V nonlinearity or not (Currently for cross-point array only)
	NL = 10;    // I-V nonlinearity in write scheme (the current ratio between Vw and Vw/2), assuming for the LTP side
	if (nonlinearIV) {  // Currently for cross-point array only
		double Vr_exp = readVoltage;  // XXX: Modify this value to Vr in the reported measurement data (can be different than readVoltage)
		// Calculation of conductance at on-chip Vr
		maxConductance = NonlinearConductance(maxConductance, NL, writeVoltageLTP, Vr_exp, readVoltage);
		minConductance = NonlinearConductance(minConductance, NL, writeVoltageLTP, Vr_exp, readVoltage);
	}
	nonlinearWrite = true;	// Consider weight update nonlinearity or not
	nonIdenticalPulse = false;	// Use non-identical pulse scheme in weight update or not
	if (nonIdenticalPulse) {
		VinitLTP = 2.85;	// Initial write voltage for LTP or weight increase (V)
		VstepLTP = 0.05;	// Write voltage step for LTP or weight increase (V)
		VinitLTD = 2.1;		// Initial write voltage for LTD or weight decrease (V)
		VstepLTD = 0.05; 	// Write voltage step for LTD or weight decrease (V)
		PWinitLTP = 75e-9;	// Initial write pulse width for LTP or weight increase (s)
		PWstepLTP = 5e-9;	// Write pulse width for LTP or weight increase (s)
		PWinitLTD = 75e-9;	// Initial write pulse width for LTD or weight decrease (s)
		PWstepLTD = 5e-9;	// Write pulse width for LTD or weight decrease (s)
		writeVoltageSquareSum = 0;	// Sum of V^2 of non-identical pulses (dynamic variable)
	}
	readNoise = false;		// Consider read noise or not
	sigmaReadNoise = 0;		// Sigma of read noise in gaussian distribution
	gaussian_dist = new std::normal_distribution<double>(0, sigmaReadNoise);	// Set up mean and stddev for read noise

	std::mt19937 localGen;	// It's OK not to use the external gen, since here the device-to-device vairation is a one-time deal
	localGen.seed(std::time(0));
	
	/* Device-to-device weight update variation */
	NL_LTP = 2.4;	// LTP nonlinearity
	NL_LTD = -4.88;	// LTD nonlinearity
	sigmaDtoD = 0;	// Sigma of device-to-device weight update vairation in gaussian distribution
	gaussian_dist2 = new std::normal_distribution<double>(0, sigmaDtoD);	// Set up mean and stddev for device-to-device weight update vairation
	paramALTP = getParamA(NL_LTP + (*gaussian_dist2)(localGen)) * maxNumLevelLTP;	// Parameter A for LTP nonlinearity
	paramALTD = getParamA(NL_LTD + (*gaussian_dist2)(localGen)) * maxNumLevelLTD;	// Parameter A for LTD nonlinearity

	/* Cycle-to-cycle weight update variation */
	sigmaCtoC = 0.035* (maxConductance - minConductance);	// Sigma of cycle-to-cycle weight update vairation: defined as the percentage of conductance range
	gaussian_dist3 = new std::normal_distribution<double>(0, sigmaCtoC);    // Set up mean and stddev for cycle-to-cycle weight update vairation

	/* Conductance range variation */
	conductanceRangeVar = false;    // Consider variation of conductance range or not
	maxConductanceVar = 0;  // Sigma of maxConductance variation (S)
	minConductanceVar = 0;  // Sigma of minConductance variation (S)
	gaussian_dist_maxConductance = new std::normal_distribution<double>(0, maxConductanceVar);
	gaussian_dist_minConductance = new std::normal_distribution<double>(0, minConductanceVar);
	if (conductanceRangeVar) {
		maxConductance += (*gaussian_dist_maxConductance)(localGen);
		minConductance += (*gaussian_dist_minConductance)(localGen);
		if (minConductance >= maxConductance || maxConductance < 0 || minConductance < 0 ) {    // Conductance variation check
			puts("[Error] Conductance variation check not passed. The variation may be too large.");
			exit(-1);
		}
		// Use the code below instead for re-choosing the variation if the check is not passed
		//do {
		//  maxConductance = avgMaxConductance + (*gaussian_dist_maxConductance)(localGen);
		//  minConductance = avgMinConductance + (*gaussian_dist_minConductance)(localGen);
		//} while (minConductance >= maxConductance || maxConductance < 0 || minConductance < 0);
	}
 
        heightInFeatureSize = cmosAccess? 4 : 2; // Cell height = 4F (Pseudo-crossbar) or 2F (cross-point)
        widthInFeatureSize = cmosAccess? (FeFET? 6 : 4) : 2; //// Cell width = 6F (FeFET) or 4F (Pseudo-crossbar) or 2F (cross-point)
}

double RealDevice::Read(double voltage) {	// Return read current (A)
	extern std::mt19937 gen;
	if (nonlinearIV) {
		// TODO: nonlinear read
		if (readNoise) {
			return voltage * conductance * (1 + (*gaussian_dist)(gen));
		} else {
			return voltage * conductance;
		}
	} else {
		if (readNoise) {
			return voltage * conductance * (1 + (*gaussian_dist)(gen));
		} else {
			return voltage * conductance;
		}
	}
}

void RealDevice::Write(double deltaWeightNormalized, double weight, double minWeight, double maxWeight) {
	double conductanceNew = conductance;	// =conductance if no update
	if (deltaWeightNormalized > 0) {	// LTP
		deltaWeightNormalized = deltaWeightNormalized/(maxWeight-minWeight);
		deltaWeightNormalized = truncate(deltaWeightNormalized, maxNumLevelLTP);
		numPulse = deltaWeightNormalized * maxNumLevelLTP;
		if (nonlinearWrite) {
			paramBLTP = (maxConductance - minConductance) / (1 - exp(-maxNumLevelLTP/paramALTP));
			xPulse = InvNonlinearWeight(conductance, maxNumLevelLTP, paramALTP, paramBLTP, minConductance);
			conductanceNew = NonlinearWeight(xPulse+numPulse, maxNumLevelLTP, paramALTP, paramBLTP, minConductance);
		} else {
			xPulse = (conductance - minConductance) / (maxConductance - minConductance) * maxNumLevelLTP;
			conductanceNew = (xPulse+numPulse) / maxNumLevelLTP * (maxConductance - minConductance) + minConductance;
		}
	} else {	// LTD
		deltaWeightNormalized = deltaWeightNormalized/(maxWeight-minWeight);
		deltaWeightNormalized = truncate(deltaWeightNormalized, maxNumLevelLTD);
		numPulse = deltaWeightNormalized * maxNumLevelLTD;
		if (nonlinearWrite) {
			paramBLTD = (maxConductance - minConductance) / (1 - exp(-maxNumLevelLTD/paramALTD));
			xPulse = InvNonlinearWeight(conductance, maxNumLevelLTD, paramALTD, paramBLTD, minConductance);
			conductanceNew = NonlinearWeight(xPulse+numPulse, maxNumLevelLTD, paramALTD, paramBLTD, minConductance);
		} else {
			xPulse = (conductance - minConductance) / (maxConductance - minConductance) * maxNumLevelLTD;
			conductanceNew = (xPulse+numPulse) / maxNumLevelLTD * (maxConductance - minConductance) + minConductance;
		}
	}

	/* Cycle-to-cycle variation */
	extern std::mt19937 gen;
	if (sigmaCtoC && numPulse != 0) {
		conductanceNew += (*gaussian_dist3)(gen) * sqrt(abs(numPulse));	// Absolute variation
	}
	
	if (conductanceNew > maxConductance) {
		conductanceNew = maxConductance;
	} else if (conductanceNew < minConductance) {
		conductanceNew = minConductance;
	}

	/* Write latency calculation */
	if (!nonIdenticalPulse) {	// Identical write pulse scheme
		if (numPulse > 0) { // LTP
			writeLatencyLTP = numPulse * writePulseWidthLTP;
			writeLatencyLTD = 0;
		} else {    // LTD
			writeLatencyLTP = 0;
			writeLatencyLTD = -numPulse * writePulseWidthLTD;
		}
	} else {	// Non-identical write pulse scheme
		writeLatencyLTP = 0;
		writeLatencyLTD = 0;
		writeVoltageSquareSum = 0;
		double V = 0;
		double PW = 0;
		if (numPulse > 0) { // LTP
			for (int i=0; i<numPulse; i++) {
				V = VinitLTP + (xPulse+i) * VstepLTP;
				PW = PWinitLTP + (xPulse+i) * PWstepLTP;
				writeLatencyLTP += PW;
				writeVoltageSquareSum += V * V;
			}
			writePulseWidthLTP = writeLatencyLTP / numPulse;
		} else {    // LTD
			for (int i=0; i<(-numPulse); i++) {
				V = VinitLTD + (maxNumLevelLTD-xPulse+i) * VstepLTD;
				PW = PWinitLTD + (maxNumLevelLTD-xPulse+i) * PWstepLTD;
				writeLatencyLTD += PW;
				writeVoltageSquareSum += V * V;
			}
			writePulseWidthLTD = writeLatencyLTD / (-numPulse);
		}
	}
	
	conductancePrev = conductance;
	conductance = conductanceNew;
}

/* Measured device */
MeasuredDevice::MeasuredDevice(int x, int y) {
	this->x = x; this->y = y;	// Cell location: x (column) and y (row) start from index 0
	readVoltage = 0.5;	// On-chip read voltage (Vr) (V)
	readPulseWidth = 5e-9;	// Read pulse width (s) (will be determined by ADC)
	writeVoltageLTP = 2;	// Write voltage (V) for LTP or weight increase
	writeVoltageLTD = 2;	// Write voltage (V) for LTD or weight decrease
	writePulseWidthLTP = 100e-9;	// Write pulse width (s) for LTP or weight increase
	writePulseWidthLTD = 100e-9;	// Write pulse width (s) for LTD or weight decrease
	writeEnergy = 0;	// Dynamic variable for calculation of write energy (J)
	numPulse = 0;	// Number of write pulses used in the most recent write operation (dynamic variable)
	cmosAccess = true;	// True: Pseudo-crossbar (1T1R), false: cross-point
	FeFET = false;		// True: FeFET structure (Pseudo-crossbar only, should be cmosAccess=1)
	gateCapFeFET = 2.1717e-18;	// Gate capacitance of FeFET (F)
	resistanceAccess = 15e3;	// The resistance of transistor (Ohm) in Pseudo-crossbar array when turned ON
	nonlinearIV = false;	// Currently for cross-point array only
	nonlinearWrite = false;	// Consider weight update nonlinearity or not
	nonIdenticalPulse = false;	// Use non-identical pulse scheme in weight update or not
	if (nonIdenticalPulse) {
		VinitLTP = 2.85;    // Initial write voltage for LTP or weight increase (V)
		VstepLTP = 0.05;    // Write voltage step for LTP or weight increase (V)
		VinitLTD = 2.1;     // Initial write voltage for LTD or weight decrease (V)
		VstepLTD = 0.05;    // Write voltage step for LTD or weight decrease (V)
		PWinitLTP = 75e-9;  // Initial write pulse width for LTP or weight increase (s)
		PWstepLTP = 5e-9;   // Write pulse width for LTP or weight increase (s)
		PWinitLTD = 75e-9;  // Initial write pulse width for LTD or weight decrease (s)
		PWstepLTD = 5e-9;   // Write pulse width for LTD or weight decrease (s)
		writeVoltageSquareSum = 0;  // Sum of V^2 of non-identical pulses (dynamic variable)
	}
	readNoise = false;		// Consider read noise or not
	sigmaReadNoise = 0.0289;	// Sigma of read noise in gaussian distribution
	NL = 10;	// Nonlinearity in write scheme (the current ratio between Vw and Vw/2), assuming for the LTP side
	gaussian_dist = new std::normal_distribution<double>(0, sigmaReadNoise);    // Set up mean and stddev for read noise
	symLTPandLTD = false;	// True: use LTP conductance data for LTD

	/* LTP */
	double rawDataConductanceLTP[] = {0,1.00e-09,2.00e-09,3.00e-09,4.00e-09,5.00e-09,6.00e-09,7.00e-09,8.00e-09,9.00e-09,1.00e-08,1.10e-08,1.20e-08,1.30e-08,1.40e-08,1.50e-08,1.60e-08,1.70e-08,1.80e-08,1.90e-08,2.00e-08,2.10e-08,2.20e-08,2.30e-08,2.40e-08,2.50e-08,2.60e-08,2.70e-08,2.80e-08,2.90e-08,3.00e-08,3.10e-08,3.20e-08,3.30e-08,3.40e-08,3.50e-08,3.60e-08,3.70e-08,3.80e-08,3.90e-08,4.00e-08,4.10e-08,4.20e-08,4.30e-08,4.40e-08,4.50e-08,4.60e-08,4.70e-08,4.80e-08,4.90e-08,5.00e-08,5.10e-08,5.20e-08,5.30e-08,5.40e-08,5.50e-08,5.60e-08,5.70e-08,5.80e-08,5.90e-08,6.00e-08,6.10e-08,6.20e-08,6.30e-08};
	dataConductanceLTP.insert(dataConductanceLTP.begin(), rawDataConductanceLTP, rawDataConductanceLTP + sizeof(rawDataConductanceLTP)/sizeof(rawDataConductanceLTP[0]));	// Put the raw data into a member variable of vector
	maxNumLevelLTP = dataConductanceLTP.size() - 1;
	/* LTD */
	if (symLTPandLTD) {	// Use LTP conductance data for LTD
		for (int i=maxNumLevelLTP; i>=0; i--) {
			dataConductanceLTD.push_back(dataConductanceLTP[i]);
		}
		maxNumLevelLTD = dataConductanceLTD.size() - 1;
	} else {	// Use provided LTD conductance data
		double rawDataConductanceLTD[] = {6.30e-08,6.20e-08,6.10e-08,6.00e-08,5.90e-08,5.80e-08,5.70e-08,5.60e-08,5.50e-08,5.40e-08,5.30e-08,5.20e-08,5.10e-08,5.00e-08,4.90e-08,4.80e-08,4.70e-08,4.60e-08,4.50e-08,4.40e-08,4.30e-08,4.20e-08,4.10e-08,4.00e-08,3.90e-08,3.80e-08,3.70e-08,3.60e-08,3.50e-08,3.40e-08,3.30e-08,3.20e-08,3.10e-08,3.00e-08,2.90e-08,2.80e-08,2.70e-08,2.60e-08,2.50e-08,2.40e-08,2.30e-08,2.20e-08,2.10e-08,2.00e-08,1.90e-08,1.80e-08,1.70e-08,1.60e-08,1.50e-08,1.40e-08,1.30e-08,1.20e-08,1.10e-08,1.00e-08,9.00e-09,8.00e-09,7.00e-09,6.00e-09,5.00e-09,4.00e-09,3.00e-09,2.00e-09,1.00e-09,0};
		dataConductanceLTD.insert(dataConductanceLTD.begin(), rawDataConductanceLTD, rawDataConductanceLTD + sizeof(rawDataConductanceLTD)/sizeof(rawDataConductanceLTD[0]));	// Put the raw data into a member variable of vector
		maxNumLevelLTD = dataConductanceLTD.size() - 1;
	}
	/* Define max/min/initial conductance */
	maxConductance = (dataConductanceLTP.back() > dataConductanceLTD.front())? dataConductanceLTD.front() : dataConductanceLTP.back();      // The last conductance point of LTP or the first conductance point of LTD, depending on which one is smaller
	minConductance = (dataConductanceLTP.front() > dataConductanceLTD.back())? dataConductanceLTP.front() : dataConductanceLTD.back();  // The first conductance point of LTP or the last conductance point of LTD, depending on which one is larger
	avgMaxConductance = maxConductance; // Average maximum cell conductance (S)
	avgMinConductance = minConductance; // Average minimum cell conductance (S)
	conductance = minConductance;
	conductancePrev = conductance;

	// Data check
	/* Check if the conductance range of LTP and LTD are consistent */
	if (dataConductanceLTP.back() != dataConductanceLTD.front() || dataConductanceLTP.front() != dataConductanceLTD.back()) {
		puts("[Error] Conductance range of LTP and LTD are not consistent");
		exit(-1);
	}
	/* Check if LTP conductance is monotonically increasing */
	for (int i=1; i<=maxNumLevelLTP; i++) {
		if (dataConductanceLTP[i] - dataConductanceLTP[i-1] <= 0) {
			puts("[Error] LTP conductance should be monotonically increasing");
			exit(-1);
		}
	}
	/* Check if LTD conductance is monotonically decreasing */
	for (int i=1; i<=maxNumLevelLTD; i++) {
		if (dataConductanceLTD[i] - dataConductanceLTD[i-1] >= 0) {
			puts("[Error] LTD conductance should be monotonically decreasing");
			exit(-1);
		}
	}

	heightInFeatureSize = cmosAccess? 4 : 2;	// Cell height = 4F (Pseudo-crossbar) or 2F (cross-point)
	widthInFeatureSize = cmosAccess? (FeFET? 6 : 4) : 2;	// Cell width = 6F (FeFET) or 4F (Pseudo-crossbar) or 2F (cross-point)
}

double MeasuredDevice::Read(double voltage) {	// Return read current (A)
	extern std::mt19937 gen;
	if (nonlinearIV) {
		// TODO: nonlinear read
		if (readNoise) {
			return voltage * conductance * (1 + (*gaussian_dist)(gen));
		} else {
			return voltage * conductance;
		}
	} else {
		if (readNoise) {
			return voltage * conductance * (1 + (*gaussian_dist)(gen));
		} else {
			return voltage * conductance;
		}
	}
}

void MeasuredDevice::Write(double deltaWeightNormalized, double weight, double minWeight, double maxWeight) {
	double conductanceNew;
	if (deltaWeightNormalized > 0) {    // LTP
		deltaWeightNormalized = deltaWeightNormalized/(maxWeight-minWeight);
		deltaWeightNormalized = truncate(deltaWeightNormalized, maxNumLevelLTP);
		numPulse = deltaWeightNormalized * maxNumLevelLTP;
		if (nonlinearWrite) {
			xPulse = InvMeasuredLTP(conductance, maxNumLevelLTP, dataConductanceLTP);
			conductanceNew = MeasuredLTP(xPulse+numPulse, maxNumLevelLTP, dataConductanceLTP);
		} else {
			xPulse = (conductance - minConductance) / (maxConductance - minConductance) * maxNumLevelLTP;
			conductanceNew = (weight-minWeight)/(maxWeight-minWeight) * (maxConductance - minConductance) + minConductance;
			if (conductanceNew > maxConductance) {
				conductanceNew = maxConductance;
			}
		}
	} else {    // LTD
		deltaWeightNormalized = deltaWeightNormalized/(maxWeight-minWeight);
		deltaWeightNormalized = truncate(deltaWeightNormalized, maxNumLevelLTD);
		numPulse = deltaWeightNormalized * maxNumLevelLTD;
		if (nonlinearWrite) {
			xPulse = InvMeasuredLTP(conductance, maxNumLevelLTP, dataConductanceLTP);
			conductanceNew = MeasuredLTP(xPulse+numPulse, maxNumLevelLTP, dataConductanceLTP);	// Use xPulse-numPulse here because the conductance will decrease with larger pulse position in dataConductanceLTD
		} else {
			xPulse = (conductance - minConductance) / (maxConductance - minConductance) * maxNumLevelLTD;
			conductanceNew = (weight-minWeight)/(maxWeight-minWeight) * (maxConductance - minConductance) + minConductance;
			if (conductanceNew < minConductance) {
				conductanceNew = minConductance;
			}
		}
	}

	/* Write latency calculation */
	if (!nonIdenticalPulse) {   // Identical write pulse scheme
		if (numPulse > 0) { // LTP
			writeLatencyLTP = numPulse * writePulseWidthLTP;
			writeLatencyLTD = 0;
		} else {    // LTD
			writeLatencyLTP = 0;
			writeLatencyLTD = -numPulse * writePulseWidthLTD;
		}
	} else {    // Non-identical write pulse scheme
		writeLatencyLTP = 0;
		writeLatencyLTD = 0;
		writeVoltageSquareSum = 0;
		double V = 0;
		double PW = 0;
		if (numPulse > 0) { // LTP
			for (int i=0; i<numPulse; i++) {
				V = VinitLTP + (xPulse+i) * VstepLTP;
				PW = PWinitLTP + (xPulse+i) * PWstepLTP;
				writeLatencyLTP += PW;
				writeVoltageSquareSum += V * V;
			}
			writePulseWidthLTP = writeLatencyLTP / numPulse;
		} else {    // LTD
			for (int i=0; i<(-numPulse); i++) {
				V = VinitLTD + (maxNumLevelLTD-xPulse+i) * VstepLTD;
				PW = PWinitLTD + (maxNumLevelLTD-xPulse+i) * PWstepLTD;
				writeLatencyLTD += PW;
				writeVoltageSquareSum += V * V;
			}
			writePulseWidthLTD = writeLatencyLTD / (-numPulse);
		}
	}

	conductancePrev = conductance;
	conductance = conductanceNew;
}

/* SRAM */
SRAM::SRAM(int x, int y) {
	this->x = x; this->y = y;
	bit = 0;	// Stored bit (1 or 0) (dynamic variable)
	bitPrev = 0;	// Previous bit
	heightInFeatureSize = 14.6;	// Cell height in terms of feature size (F)
	widthInFeatureSize = 10;	// Cell width in terms of feature size (F)
	widthSRAMCellNMOS = 2.08;	// Pull-down NMOS width in terms of feature size (F)
	widthSRAMCellPMOS = 1.23;	// Pull-up PMOS width in terms of feature size (F)
	widthAccessCMOS = 1.31;		// Access transistor width in terms of feature size (F)
	minSenseVoltage = 0.1;		// Minimum voltage difference (V) for sensing
	readEnergy = 0;				// Dynamic variable for calculation of read energy (J)
	writeEnergy = 0;			// Dynamic variable for calculation of write energy (J)
	readEnergySRAMCell = 0;		// Read energy (J) per SRAM cell (currently not used, it is included in the peripheral circuits of SRAM array in NeuroSim)
	writeEnergySRAMCell = 0;	// Write energy (J) per SRAM cell (will be obtained from NeuroSim)
}

/* Digital eNVM */
DigitalNVM::DigitalNVM(int x, int y) {
	this->x = x; this->y = y;	// Cell location: x (column) and y (row) start from index 0	
	bit = 0;	// Stored bit (1 or 0) (dynamic variable), for internel check only and not be used for read
	bitPrev = 0;	// Previous bit
	maxConductance = 5e-6;		// Maximum cell conductance (S)
	minConductance = 100e-9;	// Minimum cell conductance (S)
	avgMaxConductance = maxConductance; // Average maximum cell conductance (S)
	avgMinConductance = minConductance; // Average minimum cell conductance (S)
	conductance = minConductance;	// Current conductance (S) (dynamic variable)
	conductancePrev = conductance;	// Previous conductance (S) (dynamic variable)
	readVoltage = 0.5;	// On-chip read voltage (Vr) (V)
	readPulseWidth = 5e-9;	// Read pulse width (s) (will be determined by S/A)
	writeVoltageLTP = 2.5;	// Write voltage (V) for LTP or weight increase
	writeVoltageLTD = 2.5;	// Write voltage (V) for LTD or weight decrease
	writePulseWidthLTP = 10e-9;	// Write pulse width (s) for LTP or weight increase
	writePulseWidthLTD = 10e-9;	// Write pulse width (s) for LTD or weight decrease
	readEnergy = 0;		// Read pulse width (s) (currently not used)
	writeEnergy = 0;    // Dynamic variable for calculation of write energy (J)
	cmosAccess = true;	// True: Pseudo-crossbar (1T1R), false: cross-point
    parallelRead = false; // if it is a parallel readout scheme
	resistanceAccess = 15e3;	// The resistance of transistor (Ohm) in Pseudo-crossbar array when turned ON
	nonlinearIV = false;	// Consider I-V nonlinearity or not (Currently for cross-point array only)
	NL = 10;    // Nonlinearity in write scheme (the current ratio between Vw and Vw/2), assuming for the LTP side
	if (nonlinearIV) {  // Currently for cross-point array only
		double Vr_exp = readVoltage;  // XXX: Modify this value to Vr in the reported measurement data (can be different than readVoltage)
		// Calculation of conductance at on-chip Vr
		maxConductance = NonlinearConductance(maxConductance, NL, writeVoltageLTP, Vr_exp, readVoltage);
		minConductance = NonlinearConductance(minConductance, NL, writeVoltageLTP, Vr_exp, readVoltage);
	}
	readNoise = false;		// Consider read noise or not
	sigmaReadNoise = 0.25;	// Sigma of read noise in gaussian distribution
	gaussian_dist = new std::normal_distribution<double>(0, sigmaReadNoise);    // Set up mean and stddev for read noise
  if(cmosAccess){ // the reference current for 1T1R cell, should include the resistance
      double Rmax=1/maxConductance;
      double Rmin=1/minConductance;
      refCurrent = readVoltage/(0.5*(Rmax+Rmin+2*resistanceAccess));
  }
  else { // the reference current for cross-point array
      refCurrent = readVoltage * (avgMaxConductance + avgMinConductance) / 2;	// Set up reference current for sensing       
  }

	/* Conductance range variation */
	conductanceRangeVar =false;    // Consider variation of conductance range or not
	maxConductanceVar = 0.01*maxConductance;  // Sigma of maxConductance variation (S)
	minConductanceVar = 0.01*minConductance;  // Sigma of minConductance variation (S)
	std::mt19937 localGen;
	localGen.seed(std::time(0));
	gaussian_dist_maxConductance = new std::normal_distribution<double>(0, maxConductanceVar);
	gaussian_dist_minConductance = new std::normal_distribution<double>(0, minConductanceVar);
	if (conductanceRangeVar) {
		maxConductance += (*gaussian_dist_maxConductance)(localGen);
		minConductance += (*gaussian_dist_minConductance)(localGen);
	if (minConductance >= maxConductance || maxConductance < 0 || minConductance < 0 ) {    // Conductance variation check
			puts("[Error] Conductance variation check not passed. The variation may be too large.");
			exit(-1);
		}
	}

	heightInFeatureSize = cmosAccess? 4 : 2;	// Cell height = 4F (1T1R) or 2F (cross-point)
	widthInFeatureSize = cmosAccess? 4 : 2;	// Cell width = 4F (1T1R) or 2F (cross-point)
}

double DigitalNVM::Read(double voltage) {	// Return read current (A)
	extern std::mt19937 gen;
	if (nonlinearIV) {
		// TODO: nonlinear read
		if (readNoise) {
			return voltage * conductance * (1 + (*gaussian_dist)(gen));
		} else {
			return voltage * conductance;
		}
	} else {
		if (readNoise) {
			return voltage * conductance * (1 + (*gaussian_dist)(gen));
		} else {
			return voltage * conductance;
		}
	}
}

void DigitalNVM::Write(int bitNew, double wireCapCol) {
	double conductanceNew;
	if (nonlinearIV) {  // Currently only for cross-point array
		if (bitNew == 1) {  // SET
			conductanceNew = maxConductance;
		} else {    // RESET
			conductanceNew = minConductance;
		}
		/* I-V nonlinearity */
		conductanceAtVwLTP = NonlinearConductance(conductance, NL, writeVoltageLTP, readVoltage, writeVoltageLTP);
		conductanceAtHalfVwLTP = NonlinearConductance(conductance, NL, writeVoltageLTP, readVoltage, writeVoltageLTP/2);
		conductanceAtVwLTD = NonlinearConductance(conductance, NL, writeVoltageLTD, readVoltage, writeVoltageLTD);
		conductanceAtHalfVwLTD = NonlinearConductance(conductance, NL, writeVoltageLTD, readVoltage, writeVoltageLTD/2);
		double conductanceNewAtVwLTP = NonlinearConductance(conductanceNew, NL, writeVoltageLTP, readVoltage, writeVoltageLTP);
		double conductanceNewAtHalfVwLTP = NonlinearConductance(conductanceNew, NL, writeVoltageLTP, readVoltage, writeVoltageLTP/2);
		double conductanceNewAtVwLTD = NonlinearConductance(conductanceNew, NL, writeVoltageLTD, readVoltage, writeVoltageLTD);
		double conductanceNewAtHalfVwLTD = NonlinearConductance(conductanceNew, NL, writeVoltageLTD, readVoltage, writeVoltageLTD/2);
		if (bitNew == 1 && bit == 0) {  // SET
			writeEnergy = writeVoltageLTP * writeVoltageLTP * (conductanceAtVwLTP + conductanceNewAtVwLTP)/2 * writePulseWidthLTP;    // Selected cell in SET phase
			writeEnergy += writeVoltageLTP * writeVoltageLTP * wireCapCol;  // Charging the cap of selected columns
			writeEnergy += writeVoltageLTD/2 * writeVoltageLTD/2 * conductanceNewAtHalfVwLTD * writePulseWidthLTD;    // Half-selected during RESET phase (use the new conductance value if SET phase is before RESET phase)
			writeEnergy += writeVoltageLTD/2 * writeVoltageLTD/2 * wireCapCol;
		} else if (bitNew == 0 && bit == 1) {    // RESET
			writeEnergy = writeVoltageLTP/2 * writeVoltageLTP/2 * conductanceAtHalfVwLTP * writePulseWidthLTP; // Half-selected during SET phase (use the old conductance value if SET phase is before RESET phase)
			writeEnergy += writeVoltageLTP/2 * writeVoltageLTP/2 * wireCapCol;
			writeEnergy += writeVoltageLTD * writeVoltageLTD * wireCapCol;  // Charging the cap of selected columns
			writeEnergy += writeVoltageLTD * writeVoltageLTD * (conductanceAtVwLTD + conductanceNewAtVwLTD)/2 * writePulseWidthLTD;    // Selected cell in RESET phase
		} else {	// Half-selected
			writeEnergy = writeVoltageLTP/2 * writeVoltageLTP/2 * conductanceAtHalfVwLTP * writePulseWidthLTP; // Half-selected during SET phase
			writeEnergy += writeVoltageLTP/2 * writeVoltageLTP/2 * wireCapCol;
			writeEnergy += writeVoltageLTD/2 * writeVoltageLTD/2 * conductanceAtHalfVwLTD * writePulseWidthLTD;	// Half-selected during RESET phase
			writeEnergy += writeVoltageLTD/2 * writeVoltageLTD/2 * wireCapCol;
		}
		/* Update the nonlinear conductances with new values */
		conductanceAtVwLTP = conductanceNewAtVwLTP;
		conductanceAtHalfVwLTP = conductanceNewAtHalfVwLTP;
		conductanceAtVwLTD = conductanceNewAtVwLTD;
		conductanceAtHalfVwLTD = conductanceNewAtHalfVwLTD;
	} else {    // If not cross-point array or not considering I-V nonlinearity
		if (bitNew == 1 && bit == 0) {	// SET
			/* Normal 1T1R */
			conductanceNew = maxConductance;
			writeEnergy = writeVoltageLTP * writeVoltageLTP * (conductance + conductanceNew)/2 * writePulseWidthLTP;	// Selected cell in SET phase
			writeEnergy += writeVoltageLTP * writeVoltageLTP * wireCapCol;	// Charging the cap of selected columns
			if (!cmosAccess) {	// Cross-point
				writeEnergy += writeVoltageLTD/2 * writeVoltageLTD/2 * conductanceNew * writePulseWidthLTD;    // Half-selected during RESET phase (use the new conductance value if SET phase is before RESET phase)
				writeEnergy += writeVoltageLTD/2 * writeVoltageLTD/2 * wireCapCol;
			}
		} else if (bitNew == 0 && bit == 1) {	// RESET
			/* Normal 1T1R */
			conductanceNew = minConductance;
			writeEnergy = writeVoltageLTD * writeVoltageLTD * (conductance + conductanceNew)/2 * writePulseWidthLTD;    // Selected cell in RESET phase
			writeEnergy += writeVoltageLTD * writeVoltageLTD * wireCapCol;  // Charging the cap of selected columns
			if (!cmosAccess) {  // Cross-point
				writeEnergy += writeVoltageLTP/2 * writeVoltageLTP/2 * conductance * writePulseWidthLTP;	// Half-selected during SET phase (use the old conductance value if SET phase is before RESET phase)
				writeEnergy += writeVoltageLTP/2 * writeVoltageLTP/2 * wireCapCol;
			}
		} else {	// No operation
			conductanceNew = (bitNew == 1)? maxConductance : minConductance;
		}
	}
	conductancePrev = conductance;
	conductance = conductanceNew;
	bitPrev = bit;
	bit = bitNew;
}