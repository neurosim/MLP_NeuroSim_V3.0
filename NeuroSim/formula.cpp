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
*
*	Rui Liu			Email: rliu51 at asu dot edu
********************************************************************************/

#include <cmath>
#include <cstdlib>
#include <iostream>
#include "constant.h"
#include "formula.h"
#include "../Param.h"

extern Param* param;

using namespace std;

/* Beyond 22 nm technology, the value capIdealGate is the sum of capIdealGate and capOverlap and capFringe */
double CalculateGateCap(double width, Technology tech) {
	return (tech.capIdealGate + tech.capOverlap + tech.capFringe) * width   // 3 * tech.capFringe
			+ tech.phyGateLength * tech.capPolywire;
}

double CalculateGateArea(	// Calculate layout area and width of logic gate given fixed layout height
		int gateType, int numInput,
		double widthNMOS, double widthPMOS,
		double heightTransistorRegion, Technology tech,
		double *height, double *width) {
	double	ratio = widthPMOS / (widthPMOS + widthNMOS);

	double maxWidthPMOS, maxWidthNMOS;
	int maxNumPFin, maxNumNFin;	/* Max number of fins for the specified cell height */
	double unitWidthRegionP, unitWidthRegionN;
	double widthRegionP, widthRegionN;
	double heightRegionP, heightRegionN;
	int numFoldedPMOS = 1, numFoldedNMOS = 1;

	if (param->processNode >= 22) {	// Bulk 
		if (ratio == 0) {	/* no PMOS */
			maxWidthPMOS = 0;
			maxWidthNMOS = heightTransistorRegion - (MIN_POLY_EXT_DIFF + MIN_GAP_BET_FIELD_POLY/2) * 2 * tech.featureSize;
		} else if (ratio == 1) {	/* no NMOS */
			maxWidthPMOS = heightTransistorRegion - (MIN_POLY_EXT_DIFF + MIN_GAP_BET_FIELD_POLY/2) * 2 * tech.featureSize;
			maxWidthNMOS = 0;
		} else {
			maxWidthPMOS = ratio * (heightTransistorRegion - MIN_GAP_BET_P_AND_N_DIFFS * tech.featureSize - (MIN_POLY_EXT_DIFF + MIN_GAP_BET_FIELD_POLY/2) * 2 * tech.featureSize);
			maxWidthNMOS = maxWidthPMOS / ratio * (1 - ratio);
		}

		if (widthPMOS > 0) {
			if (widthPMOS <= maxWidthPMOS) { /* No folding */
				unitWidthRegionP = 2 * (POLY_WIDTH + MIN_GAP_BET_GATE_POLY) * tech.featureSize;
				heightRegionP = widthPMOS;
			} else {	/* Folding */
				numFoldedPMOS = (int)(ceil(widthPMOS / maxWidthPMOS));
				unitWidthRegionP = (numFoldedPMOS + 1) * (POLY_WIDTH + MIN_GAP_BET_GATE_POLY) * tech.featureSize;
				heightRegionP = maxWidthPMOS;
			}
		} else {
			unitWidthRegionP = 0;
			heightRegionP = 0;
		}

		if (widthNMOS > 0) {
			if (widthNMOS <= maxWidthNMOS) { /* No folding */
				unitWidthRegionN = 2 * (POLY_WIDTH + MIN_GAP_BET_GATE_POLY) * tech.featureSize;
				heightRegionN = widthNMOS;
			} else {	/* Folding */
				numFoldedNMOS = (int)(ceil(widthNMOS / maxWidthNMOS));
				unitWidthRegionN = (numFoldedNMOS + 1) * (POLY_WIDTH + MIN_GAP_BET_GATE_POLY) * tech.featureSize;
				heightRegionN = maxWidthNMOS;
			}
		} else {
			unitWidthRegionN = 0;
			heightRegionN = 0;
		}		
	
	} else { //FinFET
		if (ratio == 0) {	/* no PFinFET */
			maxNumPFin = 0;
			maxNumNFin = (int)(floor((heightTransistorRegion - (MIN_POLY_EXT_DIFF + MIN_GAP_BET_FIELD_POLY/2) * 2 * tech.featureSize) / tech.PitchFin)) + 1;
		} else if (ratio == 1) {	/* no NFinFET */
			maxNumPFin = (int)(floor((heightTransistorRegion - (MIN_POLY_EXT_DIFF + MIN_GAP_BET_FIELD_POLY/2) * 2 * tech.featureSize) / tech.PitchFin)) + 1;
			maxNumNFin = 0;
		} else {
			maxNumPFin = (int)(floor(ratio * (heightTransistorRegion - MIN_GAP_BET_P_AND_N_DIFFS * tech.featureSize - (MIN_POLY_EXT_DIFF + MIN_GAP_BET_FIELD_POLY/2) * 2 * tech.featureSize) / tech.PitchFin)) + 1;
			maxNumNFin = (int)(floor( (1 - ratio) * (heightTransistorRegion - MIN_GAP_BET_P_AND_N_DIFFS * tech.featureSize - (MIN_POLY_EXT_DIFF + MIN_GAP_BET_FIELD_POLY/2) * 2 * tech.featureSize) / tech.PitchFin)) + 1;
		}
		
		int NumPFin = (int)(ceil(widthPMOS/(2 * tech.heightFin + tech.widthFin)));
		if (NumPFin > 0) {
			if (NumPFin <= maxNumPFin) { /* No folding */
				unitWidthRegionP = 2 * (POLY_WIDTH + MIN_GAP_BET_GATE_POLY) * tech.featureSize;
				heightRegionP = (NumPFin-1) * tech.PitchFin + 2 * tech.widthFin/2;
			} else {	/* Folding */
				numFoldedPMOS = (int)(ceil(NumPFin / maxNumPFin));
				unitWidthRegionP = (numFoldedPMOS + 1) * (POLY_WIDTH + MIN_GAP_BET_GATE_POLY) * tech.featureSize;
				heightRegionP = (maxNumPFin-1) * tech.PitchFin + 2 * tech.widthFin/2;
			}
		} else {
			unitWidthRegionP = 0;
			heightRegionP = 0;
		}
		
		int NumNFin = (int)(ceil(widthNMOS/(2 * tech.heightFin + tech.widthFin)));			
		if (NumNFin > 0) {
			if (NumNFin <= maxNumNFin) { /* No folding */
				unitWidthRegionN = 2 * (POLY_WIDTH + MIN_GAP_BET_GATE_POLY) * tech.featureSize;
				heightRegionN = (NumNFin-1) * tech.PitchFin + 2 * tech.widthFin/2;
			} else {	/* Folding */
				numFoldedNMOS = (int)(ceil(NumNFin / maxNumNFin));
				unitWidthRegionN = (numFoldedNMOS + 1) * (POLY_WIDTH + MIN_GAP_BET_GATE_POLY) * tech.featureSize;
				heightRegionN = (maxNumNFin-1) * tech.PitchFin + 2 * tech.widthFin/2;
			}
		} else {
			unitWidthRegionN = 0;
			heightRegionN = 0;
		}
	}

	switch (gateType) {
	case INV:
		widthRegionP = unitWidthRegionP;
		widthRegionN = unitWidthRegionN;
		break;
	case NOR:
		if (numFoldedPMOS == 1 && numFoldedNMOS == 1) {	// Need to subtract the source/drain sharing region
			widthRegionP = unitWidthRegionP * numInput
						- (numInput-1) * tech.featureSize * (POLY_WIDTH + MIN_GAP_BET_GATE_POLY);
			widthRegionN = unitWidthRegionN * numInput
						- (numInput-1) * tech.featureSize * (POLY_WIDTH + MIN_GAP_BET_GATE_POLY);
		} else {	// If either PMOS or NMOS has folding, there is no source/drain sharing among different PMOS and NMOS devices.
			widthRegionP = unitWidthRegionP * numInput;
			widthRegionN = unitWidthRegionN * numInput;
		}
		break;
	case NAND:
		if (numFoldedPMOS == 1 && numFoldedNMOS == 1) {	// Need to subtract the source/drain sharing region
			widthRegionP = unitWidthRegionP * numInput
						- (numInput-1) * tech.featureSize * (POLY_WIDTH + MIN_GAP_BET_GATE_POLY);
			widthRegionN = unitWidthRegionN * numInput
						- (numInput-1) * tech.featureSize * (POLY_WIDTH + MIN_GAP_BET_GATE_POLY);
		} else {	// If either PMOS or NMOS has folding, there is no source/drain sharing among different PMOS and NMOS devices.
			widthRegionP = unitWidthRegionP * numInput;
			widthRegionN = unitWidthRegionN * numInput;
		}
		break;
	default:
		widthRegionN = widthRegionP = 0;
	}
	*width = MAX(widthRegionN, widthRegionP);
	*height = heightTransistorRegion;	// Fixed standard cell height

	return (*width)*(*height);
}

void CalculateGateCapacitance(
		int gateType, int numInput,
		double widthNMOS, double widthPMOS,
		double heightTransistorRegion, Technology tech,
		double *capInput, double *capOutput) {

	double	ratio = widthPMOS / (widthPMOS + widthNMOS);
	double maxWidthPMOS = 0, maxWidthNMOS = 0;
	int maxNumPFin = 0, maxNumNFin = 0;	/* Max numbers of fin for the specified cell height */
	double unitWidthDrainP = 0, unitWidthDrainN = 0;
	double unitWidthSourceP = 0, unitWidthSourceN = 0;
	double widthDrainP = 0, widthDrainN = 0;
	double heightDrainP = 0, heightDrainN = 0;
	int numFoldedPMOS = 1, numFoldedNMOS = 1;
	double widthDrainSidewallP = 0, widthDrainSidewallN = 0;
	if (param->processNode >= 22) { // Bulk
		if (ratio == 0) {	/* no PMOS */
			maxWidthPMOS = 0;
			maxWidthNMOS = heightTransistorRegion - (MIN_POLY_EXT_DIFF + MIN_GAP_BET_FIELD_POLY/2) * 2 * tech.featureSize;
		} else if (ratio == 1) {	/* no NMOS */
			maxWidthPMOS = heightTransistorRegion - (MIN_POLY_EXT_DIFF + MIN_GAP_BET_FIELD_POLY/2) * 2 * tech.featureSize;
			maxWidthNMOS = 0;
		} else {
			maxWidthPMOS = ratio * (heightTransistorRegion - MIN_GAP_BET_P_AND_N_DIFFS * tech.featureSize - (MIN_POLY_EXT_DIFF + MIN_GAP_BET_FIELD_POLY/2) * 2 * tech.featureSize);
			maxWidthNMOS = maxWidthPMOS / ratio * (1 - ratio);
		}
		
		if (widthPMOS > 0) {
			if (widthPMOS <= maxWidthPMOS) { /* No folding */
				unitWidthDrainP = tech.featureSize * MIN_GAP_BET_GATE_POLY;
				unitWidthSourceP = unitWidthDrainP;
				heightDrainP = widthPMOS;
			} else {	/* Folding */
				numFoldedPMOS = (int)(ceil(widthPMOS / maxWidthPMOS));
				unitWidthDrainP = (int)ceil((double)(numFoldedPMOS+1)/2) * tech.featureSize * MIN_GAP_BET_GATE_POLY;	// Num of drain fingers >= num of source fingers
				unitWidthSourceP = (int)floor((double)(numFoldedPMOS+1)/2) * tech.featureSize * MIN_GAP_BET_GATE_POLY;
				heightDrainP = maxWidthPMOS;
			}
		} else {
			unitWidthDrainP = 0;
			unitWidthSourceP = 0;
			heightDrainP = 0;
		}	
		if (widthNMOS > 0) {
			if (widthNMOS <= maxWidthNMOS) { /* No folding */
				unitWidthDrainN = tech.featureSize * MIN_GAP_BET_GATE_POLY;
				unitWidthSourceN = unitWidthDrainN;
				heightDrainN = widthNMOS;
			} else {	/* Folding */
				numFoldedNMOS = (int)(ceil(widthNMOS / maxWidthNMOS));
				unitWidthDrainN = (int)ceil((double)(numFoldedNMOS+1)/2) * tech.featureSize * MIN_GAP_BET_GATE_POLY;
				unitWidthSourceN = (int)floor((double)(numFoldedNMOS+1)/2) * tech.featureSize * MIN_GAP_BET_GATE_POLY;
				heightDrainN = maxWidthNMOS;
			}
		} else {
			unitWidthDrainN = 0;
			unitWidthSourceN = 0;
			heightDrainN = 0;
		}
	
	} else { //FinFET
		if (ratio == 0) {	/* no PFinFET */
			maxNumPFin = 0;
			maxNumNFin = (int)(floor((heightTransistorRegion - (MIN_POLY_EXT_DIFF + MIN_GAP_BET_FIELD_POLY/2) * 2 * tech.featureSize) / tech.PitchFin)) + 1;
		} else if (ratio == 1) {	/* no NFinFET */
			maxNumPFin = (int)(floor((heightTransistorRegion - (MIN_POLY_EXT_DIFF + MIN_GAP_BET_FIELD_POLY/2) * 2 * tech.featureSize) / tech.PitchFin)) + 1;
			maxNumNFin = 0;
		} else {
			maxNumPFin = (int)(floor(ratio * (heightTransistorRegion - MIN_GAP_BET_P_AND_N_DIFFS * tech.featureSize - (MIN_POLY_EXT_DIFF + MIN_GAP_BET_FIELD_POLY/2) * 2 * tech.featureSize) / tech.PitchFin)) + 1;
			maxNumNFin = (int)(floor( (1 - ratio) * (heightTransistorRegion - MIN_GAP_BET_P_AND_N_DIFFS * tech.featureSize - (MIN_POLY_EXT_DIFF + MIN_GAP_BET_FIELD_POLY/2) * 2 * tech.featureSize) / tech.PitchFin)) + 1;
		}

		int NumPFin = (int)(ceil(widthPMOS/(2 * tech.heightFin + tech.widthFin)));	

		if (NumPFin > 0) {
			if (NumPFin <= maxNumPFin) { /* No folding */
				unitWidthDrainP = tech.featureSize * MIN_GAP_BET_GATE_POLY;
				unitWidthSourceP = unitWidthDrainP;
				heightDrainP = (NumPFin-1) * tech.PitchFin + 2 * tech.widthFin/2;
			} else {	/* Folding */
				numFoldedPMOS = (int)(ceil(NumPFin / maxNumPFin));
				unitWidthDrainP = (int)ceil((double)(numFoldedPMOS+1)/2) * tech.featureSize * MIN_GAP_BET_GATE_POLY;
				unitWidthSourceP = (int)floor((double)(numFoldedPMOS+1)/2) * tech.featureSize * MIN_GAP_BET_GATE_POLY;
				heightDrainP = (maxNumPFin-1) * tech.PitchFin + 2 * tech.widthFin/2;
			}
		} else {
			unitWidthDrainP = 0;
			unitWidthSourceP = 0;
			heightDrainP = 0;
		}

		int NumNFin = (int)(ceil(widthNMOS/(2 * tech.heightFin + tech.widthFin)));	

		if (NumNFin > 0) {
			if (NumNFin <= maxNumNFin) { /* No folding */
				unitWidthDrainN = tech.featureSize * MIN_GAP_BET_GATE_POLY;
				unitWidthSourceN = unitWidthDrainN;
				heightDrainN = (NumNFin-1) * tech.PitchFin + 2 * tech.widthFin/2;
			} else {	/* Folding */
				numFoldedNMOS = (int)(ceil(NumNFin / maxNumNFin));
				unitWidthDrainN = (int)ceil((double)(numFoldedNMOS+1)/2) * tech.featureSize * MIN_GAP_BET_GATE_POLY;
				unitWidthSourceN = (int)floor((double)(numFoldedNMOS+1)/2) * tech.featureSize * MIN_GAP_BET_GATE_POLY;
				heightDrainN = (maxNumNFin-1) * tech.PitchFin + 2 * tech.widthFin/2;
			}
		} else {
			unitWidthDrainN = 0;
			unitWidthSourceN = 0;
			heightDrainN = 0;
		}	
		
	}	

	switch (gateType) {
	case INV:
		if (widthPMOS > 0) {
			widthDrainP = unitWidthDrainP;
			// Folding=1: both drain and source has 1 side; folding=2: drain has 2 sides and source has 0 side... etc
			widthDrainSidewallP = widthDrainP * 2 + heightDrainP * (1+(numFoldedPMOS+1)%2);
		}
		if (widthNMOS > 0) {
			widthDrainN = unitWidthDrainN;
			widthDrainSidewallN = widthDrainN * 2 + heightDrainN * (1+(numFoldedPMOS+1)%2);
		}
		break;
	case NOR:
		// If either PMOS or NMOS has folding, there is no source/drain sharing among different PMOS and NMOS devices
		if (numFoldedPMOS == 1 && numFoldedNMOS == 1) { 
			if (widthPMOS > 0) {	// No need to consider the source capacitance in series PMOS because here source and drain shares
				widthDrainP = unitWidthDrainP * numInput;
				widthDrainSidewallP = widthDrainP * 2 + heightDrainP;
			}
			if (widthNMOS > 0) {	// The number of NMOS drains is not equal to the number of NMOS in parallel because drain can share
				widthDrainN = unitWidthDrainN * (int)floor((double)(numInput+1)/2);	// Use floor: assume num of source regions >= num of drain regions
				widthDrainSidewallN = widthDrainN * 2 + heightDrainN * (1-(numInput+1)%2);
			}
		} else {
			if (widthPMOS > 0) {	// Need to consider the source capacitance in series PMOS (excluding the top one)
				widthDrainP = unitWidthDrainP * numInput + (numInput-1) * unitWidthSourceP;
				widthDrainSidewallP = widthDrainP * 2
									+ heightDrainP * (1+(numFoldedPMOS+1)%2) * numInput			// Drain sidewalls
									+ heightDrainP * (1-(numFoldedPMOS+1)%2) * (numInput-1);	// Source sidewalls
			}
			if (widthNMOS > 0) {	// Drain cannot share between different NMOS
				widthDrainN = unitWidthDrainN * numInput;
				widthDrainSidewallN = widthDrainN * 2 + heightDrainN * (1+(numFoldedNMOS+1)%2) * numInput;
			}
		}
		break;
	case NAND:
		// If either PMOS or NMOS has folding, there is no source/drain sharing among different PMOS and NMOS devices
		if (numFoldedPMOS == 1 && numFoldedNMOS == 1) {
			if (widthPMOS > 0) {  // The number of PMOS drains is not equal to the number of PMOS in parallel because drain can share
				widthDrainP = unitWidthDrainP * (int)floor((double)(numInput+1)/2); // Use floor: assume num of source regions >= num of drain regions
				widthDrainSidewallP = widthDrainP * 2 + heightDrainP * (1-(numInput+1)%2);
			}
			if (widthNMOS > 0) {  // No need to consider the source capacitance in series NMOS because here source and drain shares
				widthDrainN = unitWidthDrainN * numInput;
				widthDrainSidewallN = widthDrainN * 2 + heightDrainN;
			}
		} else {
			if (widthPMOS > 0) {  // Drain cannot share between different PMOS
				widthDrainP = unitWidthDrainP * numInput;
				widthDrainSidewallP = widthDrainP * 2 + heightDrainP * (1+(numFoldedPMOS+1)%2) * numInput;
			}
			if (widthNMOS > 0) {  // Need to consider the source capacitance in series NMOS (excluding the bottom one)
				widthDrainN = unitWidthDrainN * numInput + (numInput-1) * unitWidthSourceN;
				widthDrainSidewallN = widthDrainN * 2
									+ heightDrainN * (1+(numFoldedNMOS+1)%2) * numInput         // Drain sidewalls
									+ heightDrainN * (1-(numFoldedNMOS+1)%2) * (numInput-1);    // Source sidewalls
			}
		}
		break;
	default:
		widthDrainN = widthDrainP = widthDrainSidewallP = widthDrainSidewallN = 0;
	}
	/* Junction capacitance */
	double capDrainBottomN = widthDrainN * heightDrainN * tech.capJunction;
	double capDrainBottomP = widthDrainP * heightDrainP * tech.capJunction;

	/* Sidewall capacitance */	// FIXME
	double capDrainSidewallN, capDrainSidewallP;
	capDrainSidewallP = widthDrainSidewallP * tech.capSidewall;
	capDrainSidewallN = widthDrainSidewallN * tech.capSidewall;

	/* Drain to channel capacitance */	// FIXME
	double capDrainToChannelN = numFoldedNMOS * heightDrainN * tech.capDrainToChannel;
	double capDrainToChannelP = numFoldedPMOS * heightDrainP * tech.capDrainToChannel;

	if (capOutput)
		*(capOutput) = capDrainBottomN + capDrainBottomP + capDrainSidewallN + capDrainSidewallP + capDrainToChannelN + capDrainToChannelP;
	if (capInput)
		*(capInput) = CalculateGateCap(widthNMOS, tech) + CalculateGateCap(widthPMOS, tech);

}


double CalculateDrainCap(
		double width, int type,
		double heightTransistorRegion, Technology tech) {
	double drainCap = 0;
	if (type == NMOS)
		CalculateGateCapacitance(INV, 1, width, 0, heightTransistorRegion, tech, NULL, &drainCap);
	else //PMOS
		CalculateGateCapacitance(INV, 1, 0, width, heightTransistorRegion, tech, NULL, &drainCap);
	return drainCap;
}

double CalculateGateLeakage(
		int gateType, int numInput,
		double widthNMOS, double widthPMOS,
		double temperature, Technology tech) {
	int tempIndex = (int)temperature - 300;
	if ((tempIndex > 100) || (tempIndex < 0)) {
		cout<<"Error: Temperature is out of range"<<endl;
		exit(-1);
	}
	double *leakN = tech.currentOffNmos;
	double *leakP = tech.currentOffPmos;
	double leakageN, leakageP;
	switch (gateType) {
	case INV:
		leakageN = widthNMOS * leakN[tempIndex];
		leakageP = widthPMOS * leakP[tempIndex];
		return (leakageN + leakageP)/2;
	case NOR:
		leakageN = widthNMOS * leakN[tempIndex] * numInput;
		if (numInput == 2) {
			return AVG_RATIO_LEAK_2INPUT_NOR * leakageN;
		}
		else {
			return AVG_RATIO_LEAK_3INPUT_NOR * leakageN;
		}
	case NAND:
		leakageP = widthPMOS * leakP[tempIndex] * numInput;
		if (numInput == 2) {
			return AVG_RATIO_LEAK_2INPUT_NAND * leakageP;
		}
		else {
			return AVG_RATIO_LEAK_3INPUT_NAND * leakageP;
		}
	default:
		return 0.0;
	}
}

double CalculateOnResistance(double width, int type, double temperature, Technology tech) {
	double r;
	int tempIndex = (int)temperature - 300;
	if ((tempIndex > 100) || (tempIndex < 0)) {
		cout<<"Error: Temperature is out of range"<<endl;
		exit(-1);
	}
	if (type == NMOS)
		r = tech.effectiveResistanceMultiplier * tech.vdd / (tech.currentOnNmos[tempIndex] * width);
	else
		r = tech.effectiveResistanceMultiplier * tech.vdd / (tech.currentOnPmos[tempIndex] * width);
	return r;
}

double CalculateTransconductance(double width, int type, Technology tech) {
	double gm;
	if (type == NMOS) {
		gm = (2*tech.current_gmNmos)*width/(0.7*tech.vdd-tech.vth);
	} else {//type==PMOS
		gm = (2*tech.current_gmPmos)*width/(0.7*tech.vdd-tech.vth); 
	}
	return gm;
}

double horowitz(double tr, double beta, double rampInput, double *rampOutput) {
	double alpha = 1 / rampInput / tr;
	double vs = 0.5;	/* Normalized switching voltage */
	beta = 0.5;	// Just use beta=0.5 as CACTI because we do not want to consider gm anymore
				// Need to delete this input argument in the future
	double result = tr * sqrt(log(vs) * log(vs) + 2 * alpha * beta * (1 - vs));
	if (rampOutput)
		*rampOutput = (1 - vs) / result;
	return result;
}

double CalculatePassGateArea(	// Calculate layout area, height and width of pass gate given the number of folding on the pass gate width
								// This function is for pass gate where the cell height can change. For normal standard cells, use CalculateGateArea() where the cell height is fixed
		double widthNMOS, double widthPMOS, Technology tech, int numFold, double *height, double *width) {
	
	*width = (numFold + 1) * (POLY_WIDTH + MIN_GAP_BET_GATE_POLY) * tech.featureSize;	// No folding means numFold=1

	if (param->processNode >= 22) {	// Bulk
		*height = widthPMOS/numFold + widthNMOS/numFold + MIN_GAP_BET_P_AND_N_DIFFS * tech.featureSize 
				+ (MIN_POLY_EXT_DIFF + MIN_GAP_BET_FIELD_POLY/2) * 2 * tech.featureSize;
	} else {	// FinFET
		int totalNumPFin = (int)(ceil(widthPMOS/(2 * tech.heightFin + tech.widthFin)));
		int totalNumNFin = (int)(ceil(widthNMOS/(2 * tech.heightFin + tech.widthFin)));
		int NumPFin = (int)(ceil((double)totalNumPFin/numFold));
		int NumNFin = (int)(ceil((double)totalNumNFin/numFold));
		double heightRegionP = (NumPFin-1) * tech.PitchFin + 2 * tech.widthFin/2;
		double heightRegionN = (NumNFin-1) * tech.PitchFin + 2 * tech.widthFin/2;
		*height = heightRegionP + heightRegionN + MIN_GAP_BET_P_AND_N_DIFFS * tech.featureSize
				+ (MIN_POLY_EXT_DIFF + MIN_GAP_BET_FIELD_POLY/2) * 2 * tech.featureSize;
	}
	return (*width)*(*height);
}

double NonlinearResistance(double R, double NL, double Vw, double Vr, double V) {	// Nonlinearity is the current ratio between Vw and V, and R means the resistance at Vr
	double R_NL = R * V/Vr * pow(NL, (Vr-V)/(Vw/2));
	return R_NL;
}

