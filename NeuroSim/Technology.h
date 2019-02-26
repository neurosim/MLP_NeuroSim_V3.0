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

#ifndef TECHNOLOGY_H_
#define TECHNOLOGY_H_

#include "typedef.h"

class Technology {
public:
	Technology();
	virtual ~Technology() {}

	/* Functions */
	void PrintProperty();
	void Initialize(int _featureSizeInNano, DeviceRoadmap _deviceRoadmap);
	
	/* Properties */
	bool initialized;	/* Initialization flag */
	int featureSizeInNano; /*Process feature size, Unit: nm */
	double featureSize;	/* Process feature size, Unit: m */
	double RRAMFeatureSize;	/* Process feature size of RRAM, Unit: m */
	DeviceRoadmap deviceRoadmap;	/* HP or LSTP */
	double vdd;			/* Supply voltage, Unit: V */
	double vth;				/* Threshold voltage, Unit: V */
	double heightFin;	/* Fin height, Unit: m */
	double widthFin;	/* Fin width, Unit: m */
	double PitchFin;	/* Fin pitch, Unit: m */
	double phyGateLength;	/* Physical gate length, Unit: m */
	double capIdealGate;	/* Ideal gate capacitance, Unit: F/m */
	double capFringe;		/* Fringe capacitance, Unit: F/m */
	double capJunction;		/* Junction bottom capacitance, Cj0, Unit: F/m^2 */
	double capOverlap;		/* Overlap capacitance, Cover in MASTAR, Unit: F/m */
	double capSidewall;		/* Junction sidewall capacitance, Cjsw, Unit: F/m */
	double capDrainToChannel;	/* Junction drain to channel capacitance, Cjswg, Unit: F/m */
	double buildInPotential;	/* Bottom junction built-in potential(PB in BSIM4 model), Unit: V */	
	double pnSizeRatio;		/* PMOS to NMOS size ratio */
	double effectiveResistanceMultiplier;	/* Extra resistance due to vdsat */
	double currentOnNmos[101];		/* NMOS saturation current, Unit: A/m */
	double currentOnPmos[101];		/* PMOS saturation current, Unit: A/m */
	double currentOffNmos[101];	/* NMOS off current (from 300K to 400K), Unit: A/m */
	double currentOffPmos[101]; /* PMOS off current (from 300K to 400K), Unit: A/m */
    double current_gmNmos;		/* NMOS current at 0.7*vdd for gm calculation, Unit: A/m/V*/ 
    double current_gmPmos;		/* PMOS current at 0.7*vdd for gm calculation, Unit: A/m/V*/ 
  
	double capPolywire;	/* Poly wire capacitance, Unit: F/m */
};

#endif /* TECHNOLOGY_H_ */
