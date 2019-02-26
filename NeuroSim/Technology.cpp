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

#include <iostream>
#include <cmath>
#include "Technology.h"

using namespace std;

Technology::Technology() {
	initialized = false;
}

void Technology::Initialize(int _featureSizeInNano, DeviceRoadmap _deviceRoadmap) {
	if (initialized)
		cout << "Warning: Already initialized!" << endl;

	featureSizeInNano = _featureSizeInNano;
	featureSize = _featureSizeInNano * 1e-9;
	deviceRoadmap = _deviceRoadmap;
	if (featureSizeInNano == 130) {	
		if (deviceRoadmap == HP) {
			/* PTM model: 130nm_HP.pm, from http://ptm.asu.edu/ */
			vdd = 1.3;
			vth = 128.4855e-3;
			phyGateLength = 1.3e-7;
			capIdealGate = 6.058401e-10;
			capFringe = 6.119807e-10;
			effectiveResistanceMultiplier = 1.54;	/* from CACTI */
			current_gmNmos=3.94E+02; 
			current_gmPmos=2.61E+02;
			currentOnNmos[0]   = 0.93e3; 
			currentOnNmos[10]  = 0.91e3; 
			currentOnNmos[20]  = 0.89e3; 
			currentOnNmos[30]  = 0.87e3; 
			currentOnNmos[40]  = 0.85e3; 
			currentOnNmos[50]  = 0.83e3; 
			currentOnNmos[60]  = 0.81e3; 
			currentOnNmos[70]  = 0.79e3; 
			currentOnNmos[80]  = 0.77e3; 
			currentOnNmos[90]  = 0.75e3; 
			currentOnNmos[100] = 0.74e3; 
			currentOnPmos[0]   = 0.43e3; 
			currentOnPmos[10]  = 0.41e3; 
			currentOnPmos[20]  = 0.38e3; 
			currentOnPmos[30]  = 0.36e3; 
			currentOnPmos[40]  = 0.34e3; 
			currentOnPmos[50]  = 0.32e3; 
			currentOnPmos[60]  = 0.30e3; 
			currentOnPmos[70]  = 0.28e3; 
			currentOnPmos[80]  = 0.26e3; 
			currentOnPmos[90]  = 0.25e3; 
			currentOnPmos[100] = 0.24e3; 
			currentOffNmos[0]  = 100.00e-3; 
			currentOffNmos[10] = 119.90e-3; 
			currentOffNmos[20] = 142.20e-3; 
			currentOffNmos[30] = 167.00e-3; 
			currentOffNmos[40] = 194.30e-3; 
			currentOffNmos[50] = 224.30e-3; 
			currentOffNmos[60] = 256.80e-3; 
			currentOffNmos[70] = 292.00e-3; 
			currentOffNmos[80] = 329.90e-3; 
			currentOffNmos[90] = 370.50e-3; 
			currentOffNmos[100]= 413.80e-3; 
			currentOffPmos[0]  = 100.20e-3; 
			currentOffPmos[10] = 113.60e-3; 
			currentOffPmos[20] = 127.90e-3; 
			currentOffPmos[30] = 143.10e-3; 
			currentOffPmos[40] = 159.10e-3; 
			currentOffPmos[50] = 175.80e-3; 
			currentOffPmos[60] = 193.40e-3; 
			currentOffPmos[70] = 211.70e-3; 
			currentOffPmos[80] = 230.80e-3; 
			currentOffPmos[90] = 250.70e-3; 
			currentOffPmos[100]= 271.20e-3;
			pnSizeRatio = currentOnNmos[0]/currentOnPmos[0];
		} else  {   //(deviceRoadmap == LSTP)
			/* PTM model: 130nm_LP.pm, from http://ptm.asu.edu/ */
			vdd = 1.3;
			vth = 466.0949e-3;
			phyGateLength = 1.3e-7;
			capIdealGate = 1.8574e-9;
			capFringe = 9.530642e-10;
			effectiveResistanceMultiplier = 1.54;	/* from CACTI */
			current_gmNmos=3.87E+01;
			current_gmPmos=5.67E+01;			
			currentOnNmos[0]   = 300.70; 
			currentOnNmos[10]  = 273.40; 
			currentOnNmos[20]  = 249.40; 
			currentOnNmos[30]  = 228.40; 
			currentOnNmos[40]  = 209.90; 
			currentOnNmos[50]  = 193.50; 
			currentOnNmos[60]  = 179.00; 
			currentOnNmos[70]  = 166.00; 
			currentOnNmos[80]  = 154.40; 
			currentOnNmos[90]  = 144.00; 
			currentOnNmos[100] = 134.60; 
			currentOnPmos[0]   = 150.70;
			currentOnPmos[10]  = 136.20;
			currentOnPmos[20]  = 123.60;
			currentOnPmos[30]  = 112.70;
			currentOnPmos[40]  = 103.20;
			currentOnPmos[50]  = 94.88 ;
			currentOnPmos[60]  = 87.54 ;
			currentOnPmos[70]  = 81.04 ;
			currentOnPmos[80]  = 75.25 ;
			currentOnPmos[90]  = 70.08 ;
			currentOnPmos[100] = 65.44 ;
			currentOffNmos[0]  = 100.20e-6;
			currentOffNmos[10] = 135.90e-6;
			currentOffNmos[20] = 181.20e-6;
			currentOffNmos[30] = 237.80e-6;
			currentOffNmos[40] = 307.30e-6;
			currentOffNmos[50] = 391.90e-6;
			currentOffNmos[60] = 493.30e-6;
			currentOffNmos[70] = 613.70e-6;
			currentOffNmos[80] = 755.30e-6;
			currentOffNmos[90] = 920.20e-6;
			currentOffNmos[100]= 1111.0e-6;
			currentOffPmos[0]  = 100.20e-6; 
			currentOffPmos[10] = 132.80e-6; 
			currentOffPmos[20] = 173.00e-6; 
			currentOffPmos[30] = 221.90e-6; 
			currentOffPmos[40] = 280.70e-6; 
			currentOffPmos[50] = 350.40e-6; 
			currentOffPmos[60] = 432.20e-6; 
			currentOffPmos[70] = 527.20e-6; 
			currentOffPmos[80] = 636.80e-6; 
			currentOffPmos[90] = 761.90e-6; 
			currentOffPmos[100]= 903.80e-6;
			pnSizeRatio = currentOnNmos[0]/currentOnPmos[0];
		} 
	} else if (featureSizeInNano == 90) {
		if (deviceRoadmap == HP) {
			/* PTM model: 90nm_HP.pm, from http://ptm.asu.edu/ */
			vdd = 1.2;
			vth = 146.0217e-3;
			phyGateLength = 9.0e-8;
			capIdealGate = 5.694423e-10;
			capFringe = 5.652302e-10;
			effectiveResistanceMultiplier = 1.54;	/* from CACTI */
			current_gmNmos=4.95E+02;
			current_gmPmos=3.16E+02;			
			currentOnNmos[0]   = 1.07e3;
			currentOnNmos[10]  = 1.05e3;
			currentOnNmos[20]  = 1.03e3;
			currentOnNmos[30]  = 1.01e3;
			currentOnNmos[40]  = 0.99e3;
			currentOnNmos[50]  = 0.97e3;
			currentOnNmos[60]  = 0.95e3;
			currentOnNmos[70]  = 0.93e3;
			currentOnNmos[80]  = 0.90e3;
			currentOnNmos[90]  = 0.88e3;
			currentOnNmos[100] = 0.86e3;
			currentOnPmos[0]   = 0.54e3; 
			currentOnPmos[10]  = 0.50e3; 
			currentOnPmos[20]  = 0.47e3; 
			currentOnPmos[30]  = 0.44e3; 
			currentOnPmos[40]  = 0.41e3; 
			currentOnPmos[50]  = 0.39e3; 
			currentOnPmos[60]  = 0.37e3; 
			currentOnPmos[70]  = 0.34e3; 
			currentOnPmos[80]  = 0.32e3; 
			currentOnPmos[90]  = 0.31e3; 
			currentOnPmos[100] = 0.29e3; 
			currentOffNmos[0]  = 100.8e-3;	
			currentOffNmos[10] = 120.8e-3;	
			currentOffNmos[20] = 143.4e-3;	
			currentOffNmos[30] = 168.6e-3;	
			currentOffNmos[40] = 196.6e-3;	
			currentOffNmos[50] = 227.4e-3;	
			currentOffNmos[60] = 261.1e-3;	
			currentOffNmos[70] = 297.7e-3;	
			currentOffNmos[80] = 337.3e-3;	
			currentOffNmos[90] = 379.8e-3;	
			currentOffNmos[100]= 425.4e-3;	
			currentOffPmos[0]  = 100.00e-3;
			currentOffPmos[10] = 114.00e-3;
			currentOffPmos[20] = 128.90e-3;
			currentOffPmos[30] = 144.80e-3;
			currentOffPmos[40] = 161.60e-3;
			currentOffPmos[50] = 179.30e-3;
			currentOffPmos[60] = 197.90e-3;
			currentOffPmos[70] = 217.40e-3;
			currentOffPmos[80] = 237.90e-3;
			currentOffPmos[90] = 259.10e-3;
			currentOffPmos[100]= 281.30e-3;
			pnSizeRatio = currentOnNmos[0]/currentOnPmos[0];
		} else {
			/* PTM model: 90nm_LP.pm, from http://ptm.asu.edu/ */
			vdd = 1.2;
			vth = 501.3229e-3;
			phyGateLength = 9.0e-8;
			capIdealGate = 1.5413e-10;
			capFringe = 9.601334e-10;
			effectiveResistanceMultiplier = 1.77;	/* from CACTI */
			current_gmNmos=4.38E+01;
			current_gmPmos=5.99E+01;			
			currentOnNmos[0]   = 346.30 ;
			currentOnNmos[10]  = 314.50 ;
			currentOnNmos[20]  = 286.80 ;
			currentOnNmos[30]  = 262.50 ;
			currentOnNmos[40]  = 241.20 ;
			currentOnNmos[50]  = 222.30 ;
			currentOnNmos[60]  = 205.60 ;
			currentOnNmos[70]  = 190.80 ;
			currentOnNmos[80]  = 177.50 ;
			currentOnNmos[90]  = 165.60 ;
			currentOnNmos[100] = 155.00 ;
			currentOnPmos[0]   = 200.30 ;
			currentOnPmos[10]  = 179.50 ;
			currentOnPmos[20]  = 161.90 ;
			currentOnPmos[30]  = 146.90 ;
			currentOnPmos[40]  = 133.90 ;
			currentOnPmos[50]  = 122.60 ;
			currentOnPmos[60]  = 112.80 ;
			currentOnPmos[70]  = 104.10 ;
			currentOnPmos[80]  = 96.47  ;
			currentOnPmos[90]  = 89.68  ;
			currentOnPmos[100] = 83.62  ;
			currentOffNmos[0]  = 100.00e-6;
			currentOffNmos[10] = 135.70e-6;
			currentOffNmos[20] = 181.10e-6;
			currentOffNmos[30] = 238.00e-6;
			currentOffNmos[40] = 308.50e-6;
			currentOffNmos[50] = 394.60e-6;
			currentOffNmos[60] = 498.50e-6;
			currentOffNmos[70] = 622.60e-6;
			currentOffNmos[80] = 769.30e-6;
			currentOffNmos[90] = 941.20e-6;
			currentOffNmos[100]= 1141.0e-6;
			currentOffPmos[0]  = 100.30e-6;
			currentOffPmos[10] = 133.20e-6;
			currentOffPmos[20] = 174.20e-6;
			currentOffPmos[30] = 224.40e-6;
			currentOffPmos[40] = 285.10e-6;
			currentOffPmos[50] = 357.60e-6;
			currentOffPmos[60] = 443.40e-6;
			currentOffPmos[70] = 543.70e-6;
			currentOffPmos[80] = 660.00e-6;
			currentOffPmos[90] = 793.80e-6;
			currentOffPmos[100]= 946.40e-6;
			pnSizeRatio = currentOnNmos[0]/currentOnPmos[0];
		}
	} else if (featureSizeInNano == 65) {
		if (deviceRoadmap == HP) {
			/* PTM model: 65nm_HP.pm, from http://ptm.asu.edu/ */
			vdd = 1.1;
			vth = 166.3941e-3;
			phyGateLength = 6.5e-8;
			capIdealGate = 4.868295e-10;
			capFringe = 5.270361e-10;
			effectiveResistanceMultiplier = 1.54;	/* from CACTI */
			current_gmNmos=5.72E+02;
			current_gmPmos=3.99E+02;			
			currentOnNmos[0]   = 1.12e3;
			currentOnNmos[10]  = 1.10e3;
			currentOnNmos[20]  = 1.08e3;
			currentOnNmos[30]  = 1.06e3;
			currentOnNmos[40]  = 1.04e3;
			currentOnNmos[50]  = 1.02e3;
			currentOnNmos[60]  = 1.00e3;
			currentOnNmos[70]  = 0.98e3;
			currentOnNmos[80]  = 0.95e3;
			currentOnNmos[90]  = 0.93e3;
			currentOnNmos[100] = 0.91e3;
			currentOnPmos[0]   = 0.70e3; 
			currentOnPmos[10]  = 0.66e3; 
			currentOnPmos[20]  = 0.62e3; 
			currentOnPmos[30]  = 0.58e3; 
			currentOnPmos[40]  = 0.55e3; 
			currentOnPmos[50]  = 0.52e3; 
			currentOnPmos[60]  = 0.49e3; 
			currentOnPmos[70]  = 0.46e3; 
			currentOnPmos[80]  = 0.44e3; 
			currentOnPmos[90]  = 0.41e3; 
			currentOnPmos[100] = 0.39e3; 
			currentOffNmos[0]  = 100.00e-3;	
			currentOffNmos[10] = 119.70e-3;	
			currentOffNmos[20] = 141.90e-3;	
			currentOffNmos[30] = 166.80e-3;	
			currentOffNmos[40] = 194.40e-3;	
			currentOffNmos[50] = 224.80e-3;	
			currentOffNmos[60] = 258.10e-3;	
			currentOffNmos[70] = 294.40e-3;	
			currentOffNmos[80] = 333.60e-3;	
			currentOffNmos[90] = 375.90e-3;	
			currentOffNmos[100]= 421.20e-3;	
			currentOffPmos[0]  = 100.10e-3;
			currentOffPmos[10] = 115.20e-3;
			currentOffPmos[20] = 131.50e-3;
			currentOffPmos[30] = 149.00e-3;
			currentOffPmos[40] = 167.60e-3;
			currentOffPmos[50] = 187.40e-3;
			currentOffPmos[60] = 208.40e-3;
			currentOffPmos[70] = 230.50e-3;
			currentOffPmos[80] = 253.70e-3;
			currentOffPmos[90] = 278.10e-3;
			currentOffPmos[100]= 303.60e-3;
			pnSizeRatio = currentOnNmos[0]/currentOnPmos[0];
		} else {
			/* PTM model: 65nm_LP.pm, from http://ptm.asu.edu/ */
			vdd = 1.1;
			vth = 501.6636e-3;
			phyGateLength = 6.5e-8;
			capIdealGate = 1.1926e-9;
			capFringe = 9.62148e-10;
			effectiveResistanceMultiplier = 1.77;	/* from CACTI */
			current_gmNmos=5.90E+01;
			current_gmPmos=6.75E+01;			
			currentOnNmos[0]   = 400.00 ;
			currentOnNmos[10]  = 363.90 ;
			currentOnNmos[20]  = 332.30 ;
			currentOnNmos[30]  = 304.70 ;
			currentOnNmos[40]  = 280.40 ;
			currentOnNmos[50]  = 258.90 ;
			currentOnNmos[60]  = 239.90 ;
			currentOnNmos[70]  = 223.00 ;
			currentOnNmos[80]  = 207.90 ;
			currentOnNmos[90]  = 194.30 ;
			currentOnNmos[100] = 182.10 ;
			currentOnPmos[0]   = 238.70 ;
			currentOnPmos[10]  = 216.10 ;
			currentOnPmos[20]  = 196.60 ;
			currentOnPmos[30]  = 179.70 ;
			currentOnPmos[40]  = 164.90 ;
			currentOnPmos[50]  = 152.00 ;
			currentOnPmos[60]  = 140.50 ;
			currentOnPmos[70]  = 130.40 ;
			currentOnPmos[80]  = 121.40 ;
			currentOnPmos[90]  = 113.30 ;
			currentOnPmos[100] = 106.10 ;
			currentOffNmos[0]  = 100.20e-6;
			currentOffNmos[10] = 137.50e-6;
			currentOffNmos[20] = 185.80e-6;
			currentOffNmos[30] = 247.20e-6;
			currentOffNmos[40] = 324.20e-6;
			currentOffNmos[50] = 419.30e-6;
			currentOffNmos[60] = 535.40e-6;
			currentOffNmos[70] = 675.70e-6;
			currentOffNmos[80] = 843.100e-6;
			currentOffNmos[90] = 1041.00e-6;
			currentOffNmos[100]= 1273.00e-6;
			currentOffPmos[0]  = 100.20e-6;
			currentOffPmos[10] = 135.40e-6;
			currentOffPmos[20] = 179.70e-6;
			currentOffPmos[30] = 234.90e-6;
			currentOffPmos[40] = 302.50e-6;
			currentOffPmos[50] = 384.30e-6;
			currentOffPmos[60] = 482.20e-6;
			currentOffPmos[70] = 598.00e-6;
			currentOffPmos[80] = 733.90e-6;
			currentOffPmos[90] = 891.60e-6;
			currentOffPmos[100]= 1073.00e-6;
			pnSizeRatio = currentOnNmos[0]/currentOnPmos[0];
		}
	} else if (featureSizeInNano == 45) {
		if (deviceRoadmap == HP) {
			/* PTM model: 45nm_HP.pm, from http://ptm.asu.edu/ */
			vdd = 1.0;
			vth = 171.0969e-3;
			phyGateLength = 4.5e-8;
			capIdealGate = 4.091305e-10;
			capFringe = 4.957928e-10;
			effectiveResistanceMultiplier = 1.54;	/* from CACTI */
			current_gmNmos=7.37E+02;
			current_gmPmos=6.30E+02;			
			currentOnNmos[0]   = 1.27e3;
			currentOnNmos[10]  = 1.24e3;
			currentOnNmos[20]  = 1.22e3;
			currentOnNmos[30]  = 1.19e3;
			currentOnNmos[40]  = 1.16e3;
			currentOnNmos[50]  = 1.13e3;
			currentOnNmos[60]  = 1.11e3;
			currentOnNmos[70]  = 1.08e3;
			currentOnNmos[80]  = 1.05e3;
			currentOnNmos[90]  = 1.02e3;
			currentOnNmos[100] = 1.00e3;
			currentOnPmos[0]   = 1.08e3; 
			currentOnPmos[10]  = 1.04e3; 
			currentOnPmos[20]  = 1.00e3; 
			currentOnPmos[30]  = 0.96e3; 
			currentOnPmos[40]  = 0.92e3; 
			currentOnPmos[50]  = 0.88e3; 
			currentOnPmos[60]  = 0.85e3; 
			currentOnPmos[70]  = 0.81e3; 
			currentOnPmos[80]  = 0.78e3; 
			currentOnPmos[90]  = 0.75e3; 
			currentOnPmos[100] = 0.72e3; 
			currentOffNmos[0]  = 100.00e-3;	
			currentOffNmos[10] = 120.70e-3;	
			currentOffNmos[20] = 144.10e-3;	
			currentOffNmos[30] = 170.50e-3;	
			currentOffNmos[40] = 199.80e-3;	
			currentOffNmos[50] = 232.30e-3;	
			currentOffNmos[60] = 268.00e-3;	
			currentOffNmos[70] = 307.10e-3;	
			currentOffNmos[80] = 349.50e-3;	
			currentOffNmos[90] = 395.40e-3;	
			currentOffNmos[100]= 444.80e-3;	
			currentOffPmos[0]  = 100.20e-3;
			currentOffPmos[10] = 118.70e-3;
			currentOffPmos[20] = 139.30e-3;
			currentOffPmos[30] = 162.00e-3;
			currentOffPmos[40] = 186.80e-3;
			currentOffPmos[50] = 213.90e-3;
			currentOffPmos[60] = 243.30e-3;
			currentOffPmos[70] = 274.90e-3;
			currentOffPmos[80] = 308.90e-3;
			currentOffPmos[90] = 345.20e-3;
			currentOffPmos[100]= 383.80e-3;
			pnSizeRatio = currentOnNmos[0]/currentOnPmos[0];
		} else {
			/* PTM model: 45nm_LP.pm, from http://ptm.asu.edu/ */
			vdd = 1;
			vth = 464.3718e-3;
			phyGateLength = 4.5e-8;
			capIdealGate = 8.930709e-10;
			capFringe = 8.849901e-10;
			effectiveResistanceMultiplier = 1.77;	/* from CACTI */
			current_gmNmos=1.32E+02;
			current_gmPmos=8.65E+01;			
			currentOnNmos[0]   = 500.20 ;
			currentOnNmos[10]  = 462.00 ;
			currentOnNmos[20]  = 427.80 ;
			currentOnNmos[30]  = 397.10 ;
			currentOnNmos[40]  = 369.40 ;
			currentOnNmos[50]  = 344.50 ;
			currentOnNmos[60]  = 322.10 ;
			currentOnNmos[70]  = 301.80 ;
			currentOnNmos[80]  = 283.40 ;
			currentOnNmos[90]  = 266.70 ;
			currentOnNmos[100] = 251.50 ;
			currentOnPmos[0]   = 300.00 ;
			currentOnPmos[10]  = 275.70 ;
			currentOnPmos[20]  = 254.20 ;
			currentOnPmos[30]  = 235.10 ;
			currentOnPmos[40]  = 218.10 ;
			currentOnPmos[50]  = 202.80 ;
			currentOnPmos[60]  = 189.20 ;
			currentOnPmos[70]  = 176.90 ;
			currentOnPmos[80]  = 165.80 ;
			currentOnPmos[90]  = 155.80 ;
			currentOnPmos[100] = 146.70 ;
			currentOffNmos[0]  = 100.00e-6;
			currentOffNmos[10] = 140.50e-6;
			currentOffNmos[20] = 193.90e-6;
			currentOffNmos[30] = 263.10e-6;
			currentOffNmos[40] = 351.40e-6;
			currentOffNmos[50] = 462.50e-6;
			currentOffNmos[60] = 600.30e-6;
			currentOffNmos[70] = 769.20e-6;
			currentOffNmos[80] = 973.900e-6;
			currentOffNmos[90] = 1219.00e-6;
			currentOffNmos[100]= 1511.00e-6;
			currentOffPmos[0]  = 100.20e-6;
			currentOffPmos[10] = 138.40e-6;
			currentOffPmos[20] = 187.60e-6;
			currentOffPmos[30] = 250.10e-6;
			currentOffPmos[40] = 328.10e-6;
			currentOffPmos[50] = 424.10e-6;
			currentOffPmos[60] = 540.90e-6;
			currentOffPmos[70] = 681.30e-6;
			currentOffPmos[80] = 848.30e-6;
			currentOffPmos[90] = 1045.00e-6;
			currentOffPmos[100]= 1275.00e-6;
			pnSizeRatio = currentOnNmos[0]/currentOnPmos[0];
		}
	} else if (featureSizeInNano == 32) {
		if (deviceRoadmap == HP) {
			/* PTM model: 32nm_HP.pm, from http://ptm.asu.edu/ */
			vdd = 0.9;
			vth = 194.4951e-3;
			phyGateLength = 3.4e-8;
			capIdealGate = 3.767721e-10;
			capFringe = 4.713762e-10;
			effectiveResistanceMultiplier = 1.54;	/* from CACTI */
			current_gmNmos=9.29E+02;
			current_gmPmos=6.73E+02;			
			currentOnNmos[0]   = 1.41e3;
			currentOnNmos[10]  = 1.38e3;
			currentOnNmos[20]  = 1.35e3;
			currentOnNmos[30]  = 1.31e3;
			currentOnNmos[40]  = 1.28e3;
			currentOnNmos[50]  = 1.25e3;
			currentOnNmos[60]  = 1.21e3;
			currentOnNmos[70]  = 1.18e3;
			currentOnNmos[80]  = 1.15e3;
			currentOnNmos[90]  = 1.12e3;
			currentOnNmos[100] = 1.08e3;
			currentOnPmos[0]   = 1.22e3; 
			currentOnPmos[10]  = 1.17e3; 
			currentOnPmos[20]  = 1.12e3; 
			currentOnPmos[30]  = 1.07e3; 
			currentOnPmos[40]  = 1.02e3; 
			currentOnPmos[50]  = 0.98e3; 
			currentOnPmos[60]  = 0.94e3; 
			currentOnPmos[70]  = 0.89e3; 
			currentOnPmos[80]  = 0.86e3; 
			currentOnPmos[90]  = 0.82e3; 
			currentOnPmos[100] = 0.78e3; 
			currentOffNmos[0]  = 100.30e-3;	
			currentOffNmos[10] = 120.40e-3;	
			currentOffNmos[20] = 143.10e-3;	
			currentOffNmos[30] = 168.60e-3;	
			currentOffNmos[40] = 197.00e-3;	
			currentOffNmos[50] = 228.40e-3;	
			currentOffNmos[60] = 262.90e-3;	
			currentOffNmos[70] = 300.60e-3;	
			currentOffNmos[80] = 341.70e-3;	
			currentOffNmos[90] = 386.10e-3;	
			currentOffNmos[100]= 433.90e-3;	
			currentOffPmos[0]  = 100.10e-3;
			currentOffPmos[10] = 119.00e-3;
			currentOffPmos[20] = 140.00e-3;
			currentOffPmos[30] = 163.30e-3;
			currentOffPmos[40] = 188.80e-3;
			currentOffPmos[50] = 216.70e-3;
			currentOffPmos[60] = 247.00e-3;
			currentOffPmos[70] = 279.70e-3;
			currentOffPmos[80] = 314.90e-3;
			currentOffPmos[90] = 352.60e-3;
			currentOffPmos[100]= 392.80e-3;
			pnSizeRatio = currentOnNmos[0]/currentOnPmos[0];
		} else {
			/* PTM model: 32nm_LP.pm, from http://ptm.asu.edu/ */
			vdd = 0.9;
			vth = 442.034e-3;
			phyGateLength = 3.4e-8;
			capIdealGate = 8.375279e-10;
			capFringe = 6.856677e-10;
			effectiveResistanceMultiplier = 1.77;	/* from CACTI */
			current_gmNmos=2.56E+02;
			current_gmPmos=1.19E+02;			
			currentOnNmos[0]   =600.20;
			currentOnNmos[10]  =562.80;
			currentOnNmos[20]  =528.20;
			currentOnNmos[30]  =496.20;
			currentOnNmos[40]  =466.80;
			currentOnNmos[50]  =439.70;
			currentOnNmos[60]  =414.80;
			currentOnNmos[70]  =391.90;
			currentOnNmos[80]  =370.70;
			currentOnNmos[90]  =351.30;
			currentOnNmos[100] =333.30;
			currentOnPmos[0]   = 400.00;
			currentOnPmos[10]  = 368.40;
			currentOnPmos[20]  = 340.30;
			currentOnPmos[30]  = 315.30;
			currentOnPmos[40]  = 292.90;
			currentOnPmos[50]  = 272.80;
			currentOnPmos[60]  = 254.80;
			currentOnPmos[70]  = 238.50;
			currentOnPmos[80]  = 223.80;
			currentOnPmos[90]  = 210.50;
			currentOnPmos[100] = 198.40;
			currentOffNmos[0]  = 100.10e-6;
			currentOffNmos[10] = 143.60e-6;
			currentOffNmos[20] = 202.10e-6;
			currentOffNmos[30] = 279.30e-6;
			currentOffNmos[40] = 379.50e-6;
			currentOffNmos[50] = 507.50e-6;
			currentOffNmos[60] = 668.80e-6;
			currentOffNmos[70] = 869.20e-6;
			currentOffNmos[80] = 1115.00e-6;
			currentOffNmos[90] = 1415.00e-6;
			currentOffNmos[100]= 1774.00e-6;
			currentOffPmos[0]  = 100.10e-6;
			currentOffPmos[10] = 140.70e-6;
			currentOffPmos[20] = 194.00e-6;
			currentOffPmos[30] = 262.50e-6;
			currentOffPmos[40] = 349.30e-6;
			currentOffPmos[50] = 457.70e-6;
			currentOffPmos[60] = 591.20e-6;
			currentOffPmos[70] = 753.70e-6;
			currentOffPmos[80] = 949.30e-6;
			currentOffPmos[90] = 1182.00e-6;
			currentOffPmos[100]= 1457.00e-6;
			pnSizeRatio = currentOnNmos[0]/currentOnPmos[0];
		}
	} else if (featureSizeInNano == 22) {
		if (deviceRoadmap == HP) {
			/* PTM model: 22nm.pm, from http://ptm.asu.edu/ */
			vdd = 0.85;
			vth = 208.9006e-3;
			phyGateLength = 2.0e-8;
			capIdealGate = 3.287e-10;
			capFringe = 4.532e-10;
			effectiveResistanceMultiplier = 1.54;	/* from CACTI */
			current_gmNmos=1.08E+03;
			current_gmPmos=6.98E+02;			
			currentOnNmos[0]   = 1.50e3;
			currentOnNmos[10]  = 1.47e3;
			currentOnNmos[20]  = 1.43e3;
			currentOnNmos[30]  = 1.39e3;
			currentOnNmos[40]  = 1.35e3;
			currentOnNmos[50]  = 1.31e3;
			currentOnNmos[60]  = 1.28e3;
			currentOnNmos[70]  = 1.24e3;
			currentOnNmos[80]  = 1.20e3;
			currentOnNmos[90]  = 1.17e3;
			currentOnNmos[100] = 1.13e3;
			currentOnPmos[0]   = 1.32e3;
			currentOnPmos[10]  = 1.25e3;
			currentOnPmos[20]  = 1.19e3;
			currentOnPmos[30]  = 1.13e3;
			currentOnPmos[40]  = 1.07e3;
			currentOnPmos[50]  = 1.02e3;
			currentOnPmos[60]  = 0.97e3;
			currentOnPmos[70]  = 0.92e3;
			currentOnPmos[80]  = 0.88e3;
			currentOnPmos[90]  = 0.84e3;
			currentOnPmos[100] = 0.80e3;
			currentOffNmos[0]  = 100.20e-3;	
			currentOffNmos[10] = 120.40e-3;	
			currentOffNmos[20] = 143.50e-3;	
			currentOffNmos[30] = 169.50e-3;	
			currentOffNmos[40] = 198.70e-3;	
			currentOffNmos[50] = 231.20e-3;	
			currentOffNmos[60] = 267.00e-3;	
			currentOffNmos[70] = 306.30e-3;	
			currentOffNmos[80] = 349.30e-3;	
			currentOffNmos[90] = 396.00e-3;	
			currentOffNmos[100]= 446.60e-3;	
			currentOffPmos[0]  = 100.20e-3;
			currentOffPmos[10] = 119.40e-3;
			currentOffPmos[20] = 140.80e-3;
			currentOffPmos[30] = 164.60e-3;
			currentOffPmos[40] = 190.90e-3;
			currentOffPmos[50] = 219.50e-3;
			currentOffPmos[60] = 250.70e-3;
			currentOffPmos[70] = 284.50e-3;
			currentOffPmos[80] = 320.90e-3;
			currentOffPmos[90] = 359.80e-3;
			currentOffPmos[100]= 401.50e-3;
			pnSizeRatio = currentOnNmos[0]/currentOnPmos[0];
		} else {
			/* PTM model: 22nm_LP.pm, from http://ptm.asu.edu/ */
			vdd = 0.85;
			vth = 419.915e-3;
			phyGateLength = 2.0e-8;
			capIdealGate = 5.245e-10;
			capFringe = 8.004e-10;
			effectiveResistanceMultiplier = 1.77;	/* from CACTI */
			current_gmNmos=4.56E+02;
			current_gmPmos=1.85E+02;			
			currentOnNmos[0]   = 791.90;
			currentOnNmos[10]  = 756.40;
			currentOnNmos[20]  = 722.20;
			currentOnNmos[30]  = 689.40;
			currentOnNmos[40]  = 658.10;
			currentOnNmos[50]  = 628.30;
			currentOnNmos[60]  = 600.00;
			currentOnNmos[70]  = 573.30;
			currentOnNmos[80]  = 548.00;
			currentOnNmos[90]  = 524.20;
			currentOnNmos[100] = 501.70;
			currentOnPmos[0]   = 600.20;
			currentOnPmos[10]  = 561.30;
			currentOnPmos[20]  = 525.50;
			currentOnPmos[30]  = 492.50;
			currentOnPmos[40]  = 462.20;
			currentOnPmos[50]  = 434.30;
			currentOnPmos[60]  = 408.70;
			currentOnPmos[70]  = 385.10;
			currentOnPmos[80]  = 363.40;
			currentOnPmos[90]  = 343.30;
			currentOnPmos[100] = 324.80;
			currentOffNmos[0]  = 100.00e-6;
			currentOffNmos[10] = 147.30e-6;
			currentOffNmos[20] = 212.10e-6;
			currentOffNmos[30] = 299.60e-6;
			currentOffNmos[40] = 415.30e-6;
			currentOffNmos[50] = 565.80e-6;
			currentOffNmos[60] = 758.90e-6;
			currentOffNmos[70] = 1003.00e-6;
			currentOffNmos[80] = 1307.00e-6;
			currentOffNmos[90] = 1682.00e-6;
			currentOffNmos[100]= 2139.00e-6;
			currentOffPmos[0]  = 100.00e-6;
			currentOffPmos[10] = 147.30e-6;
			currentOffPmos[20] = 212.10e-6;
			currentOffPmos[30] = 299.60e-6;
			currentOffPmos[40] = 415.30e-6;
			currentOffPmos[50] = 565.80e-6;
			currentOffPmos[60] = 758.90e-6;
			currentOffPmos[70] = 1003.00e-6;
			currentOffPmos[80] = 1307.00e-6;
			currentOffPmos[90] = 1682.00e-6;
			currentOffPmos[100]= 2139.00e-6;
			pnSizeRatio = currentOnNmos[0]/currentOnPmos[0];
		}
	} else if (featureSizeInNano == 14) {
		if (deviceRoadmap == HP) {
			/* PTM model: 14nfet_HP.pm,14pfet_HP.pm, from http://ptm.asu.edu/ */
			vdd = 0.8;
			vth = 172.9868e-3;
			heightFin = 2.3e-8;
			
			widthFin = 1.0e-8;
			PitchFin = 3.2e-8;			
			phyGateLength = 1.8e-8;
			capIdealGate = 1.2573e-9;
			capFringe = 0;
			effectiveResistanceMultiplier = 1.51;	/* from CACTI */
			current_gmNmos=1.33E+03;
			current_gmPmos=7.83E+02;			
			currentOnNmos[0]  = 1.6861e3;
			currentOnNmos[10] = 1.7108e3;
			currentOnNmos[20] = 1.7348e3;
			currentOnNmos[30] = 1.7583e3;
			currentOnNmos[40] = 1.7812e3;
			currentOnNmos[50] = 1.8035e3;
			currentOnNmos[60] = 1.8252e3;
			currentOnNmos[70] = 1.8464e3;
			currentOnNmos[80] = 1.8671e3;
			currentOnNmos[90] = 1.8872e3;
			currentOnNmos[100] =1.9068e3;
			currentOnPmos[0]  = 1.5504e3;
			currentOnPmos[10] = 1.5799e3;
			currentOnPmos[20] = 1.6091e3;
			currentOnPmos[30] = 1.638e3;
			currentOnPmos[40] = 1.6666e3;
			currentOnPmos[50] = 1.6948e3;
			currentOnPmos[60] = 1.7227e3;
			currentOnPmos[70] = 1.7504e3;
			currentOnPmos[80] = 1.7776e3;
			currentOnPmos[90] = 1.8045e3;
			currentOnPmos[100] =1.8311e3;
			currentOffNmos[0]  = 100.9152e-3;
			currentOffNmos[10] = 150.0786e-3;
			currentOffNmos[20] = 218.1063e-3;
			currentOffNmos[30] = 310.359e-3;
			currentOffNmos[40] = 433.1757e-3;
			currentOffNmos[50] = 593.9313e-3;
			currentOffNmos[60] = 801.0747e-3;
			currentOffNmos[70] = 1.0641;
			currentOffNmos[80] = 1.3938;
			currentOffNmos[90] = 1.8016;
			currentOffNmos[100] =2.30031;
			currentOffPmos[0]  = 98.6503e-3;
			currentOffPmos[10] = 157.5545e-3;
			currentOffPmos[20] = 245.1751e-3;
			currentOffPmos[30] = 372.5374e-3;
			currentOffPmos[40] = 553.7661e-3;
			currentOffPmos[50] = 806.5971e-3;
			currentOffPmos[60] = 1.1529;
			currentOffPmos[70] = 1.6190;
			currentOffPmos[80] = 2.2361;
			currentOffPmos[90] = 3.0405;
			currentOffPmos[100] =4.0733;
			pnSizeRatio = (int)(currentOnNmos[0]/currentOnPmos[0]);
		} else {
			/* PTM model: 14nfet_LP.pm,14pfet_LP.pm, from http://ptm.asu.edu/ */
			vdd = 0.8;
			vth = 382.1222e-3;
			heightFin = 2.3e-8;
			widthFin = 1.0e-8;
			PitchFin = 3.2e-8;
			phyGateLength = 1.8e-8;
			capIdealGate = 1.0572e-9;
			capFringe = 0;
			effectiveResistanceMultiplier = 1.76;	/* from CACTI */
			current_gmNmos=6.02E+02;
			current_gmPmos=2.36E+02;			
			currentOnNmos[0]  = 892.8173;
			currentOnNmos[10] = 911.7499;
			currentOnNmos[20] = 930.5372;
			currentOnNmos[30] = 949.1747;
			currentOnNmos[40] = 967.6577;
			currentOnNmos[50] = 985.9813;
			currentOnNmos[60] = 1.0041e3 ;
			currentOnNmos[70] = 1.0221e3;
			currentOnNmos[80] = 1.0399e3;
			currentOnNmos[90] = 1.0576e3 ;
			currentOnNmos[100] =1.0750e3;
			currentOnPmos[0]  = 819.8866;
			currentOnPmos[10] = 843.2728;
			currentOnPmos[20] = 866.7718;
			currentOnPmos[30] = 890.3688;
			currentOnPmos[40] = 914.0488;
			currentOnPmos[50] = 937.7967;
			currentOnPmos[60] = 961.5968;
			currentOnPmos[70] = 985.4335;
			currentOnPmos[80] = 1.0093e3;
			currentOnPmos[90] = 1.0332e3;
			currentOnPmos[100] =1.0570e3;
			currentOffNmos[0]  = 99.7866e-6;
			currentOffNmos[10] = 184.4553e-6;
			currentOffNmos[20] = 328.7707e-6;
			currentOffNmos[30] = 566.8658e-6;
			currentOffNmos[40] = 948.1816e-6;
			currentOffNmos[50] = 1.5425e-3;
			currentOffNmos[60] = 2.4460e-3;
			currentOffNmos[70] = 3.7885e-3;
			currentOffNmos[80] = 5.7416e-3;
			currentOffNmos[90] = 8.5281e-3;
			currentOffNmos[100] =1.24327e-2;;
			currentOffPmos[0]  = 102.3333e-6;
			currentOffPmos[10] = 203.4774e-6;
			currentOffPmos[20] = 389.0187e-6;
			currentOffPmos[30] = 717.5912e-6;
			currentOffPmos[40] = 1.2810e-3;
			currentOffPmos[50] = 2.2192e-3;
			currentOffPmos[60] = 3.7395e-3;
			currentOffPmos[70] = 6.1428e-3;
			currentOffPmos[80] = 9.8554e-3;
			currentOffPmos[90] = 1.54702e-2;
			currentOffPmos[100] =2.37959e-2;
			pnSizeRatio = (int)(currentOnNmos[0]/currentOnPmos[0]);
		}
	} else if (featureSizeInNano == 10) {
		if (deviceRoadmap == HP) {
			/* PTM model: 10nfet_HP.pm,10pfet_HP.pm, from http://ptm.asu.edu/ */
			vdd = 0.75;
			vth = 182.8605e-3;
			heightFin = 2.1e-8;
			widthFin = 9e-9;
			PitchFin = 2.8e-8;			
			phyGateLength = 1.4e-8;
			capIdealGate = 1.1418e-9;
			capFringe = 0;
			effectiveResistanceMultiplier = 1.49;	/* from CACTI */
			current_gmNmos=1.56E+03;
			current_gmPmos=8.02E+02;			
			currentOnNmos[0]  = 1.7691e3;
			currentOnNmos[10] = 1.7929e3;
			currentOnNmos[20] = 1.8162e3;
			currentOnNmos[30] = 1.8389e3;
			currentOnNmos[40] = 1.8609e3;
			currentOnNmos[50] = 1.8825e3;
			currentOnNmos[60] = 1.9035e3;
			currentOnNmos[70] = 1.9239e3;
			currentOnNmos[80] = 1.9438e3;
			currentOnNmos[90] = 1.9632e3;
			currentOnNmos[100] =1.9821e3;
			currentOnPmos[0]  = 1.6268e3;
			currentOnPmos[10] = 1.6561e3;
			currentOnPmos[20] = 1.6851e3;
			currentOnPmos[30] = 1.7138e3;
			currentOnPmos[40] = 1.7422e3;
			currentOnPmos[50] = 1.7703e3;
			currentOnPmos[60] = 1.798e3;
			currentOnPmos[70] = 1.8255e3;
			currentOnPmos[80] = 1.8525e3;
			currentOnPmos[90] = 1.8792e3;
			currentOnPmos[100] =1.9056e3;
			currentOffNmos[0]  = 100.1203e-3;
			currentOffNmos[10] = 148.6272e-3;
			currentOffNmos[20] = 215.6467e-3;
			currentOffNmos[30] = 306.4157e-3;
			currentOffNmos[40] = 427.1264e-3;
			currentOffNmos[50] = 584.9857e-3;
			currentOffNmos[60] = 788.2578e-3;
			currentOffNmos[70] = 1.0463; 
			currentOffNmos[80] = 1.3695; 
			currentOffNmos[90] = 1.7694; 
			currentOffNmos[100] =2.2584; 
			currentOffPmos[0]  = 1.6268e-3;
			currentOffPmos[10] = 157.8505e-3;
			currentOffPmos[20] = 245.2725e-3;
			currentOffPmos[30] = 372.2051e-3;
			currentOffPmos[40] = 552.6667e-3;
			currentOffPmos[50] = 804.2786e-3;
			currentOffPmos[60] = 1.1488; 
			currentOffPmos[70] = 1.6125; 
			currentOffPmos[80] = 2.2268; 
			currentOffPmos[90] = 3.0281; 
			currentOffPmos[100] =4.0584;
			pnSizeRatio = (int)(currentOnNmos[0]/currentOnPmos[0]);
		} else {
			/* PTM model: 10nfet_LP.pm,10pfet_LP.pm, from http://ptm.asu.edu/ */
			vdd = 0.75;
			vth = 390.5541e-3;
			heightFin = 2.1e-8;
			widthFin = 9e-9;
			PitchFin = 2.8e-8;			
			phyGateLength = 1.4e-8;
			capIdealGate = 9.418984e-10;
			capFringe = 0;
			effectiveResistanceMultiplier = 1.73;	/* from CACTI */
			current_gmNmos=8.22E+02;
			current_gmPmos=1.90E+02;			
			currentOnNmos[0]  = 862.4823;
			currentOnNmos[10] = 882.0505;
			currentOnNmos[20] = 901.514 ;
			currentOnNmos[30] = 920.8656;
			currentOnNmos[40] = 940.0977;
			currentOnNmos[50] = 959.2027;
			currentOnNmos[60] = 978.1731;
			currentOnNmos[70] = 997.0013;
			currentOnNmos[80] = 1.0157e3 ;
			currentOnNmos[90] = 1.0342e3 ;
			currentOnNmos[100] =1.0526e3 ;
			currentOnPmos[0]  = 774.9657;  
			currentOnPmos[10] = 799.5285;
			currentOnPmos[20] = 824.2549;
			currentOnPmos[30] = 849.1259;
			currentOnPmos[40] = 874.1225;
			currentOnPmos[50] = 899.2255;
			currentOnPmos[60] = 924.4159;
			currentOnPmos[70] = 949.6744;
			currentOnPmos[80] = 974.9818;
			currentOnPmos[90] = 1.0003e3;
			currentOnPmos[100] =1.0257e3;
			currentOffNmos[0]  = 99.6973e-6;
			currentOffNmos[10] = 184.4892e-6;
			currentOffNmos[20] = 329.1615e-6;
			currentOffNmos[30] = 568.0731e-6;
			currentOffNmos[40] = 951.0401e-6;
			currentOffNmos[50] = 1.5484e-3;
			currentOffNmos[60] = 2.4574e-3;
			currentOffNmos[70] = 3.8090e-3;
			currentOffNmos[80] = 5.7767e-3;
			currentOffNmos[90] = 8.5862e-3;
			currentOffNmos[100] =1.2525e-2;
			currentOffPmos[0]  = 100.5839e-6;
			currentOffPmos[10] = 200.2609e-6;
			currentOffPmos[20] = 383.3239e-6;
			currentOffPmos[30] = 707.8499e-6;
			currentOffPmos[40] = 1.2649e-3;
			currentOffPmos[50] = 2.1932e-3;
			currentOffPmos[60] = 3.6987e-3;
			currentOffPmos[70] = 6.0804e-3;
			currentOffPmos[80] = 9.7622e-3;
			currentOffPmos[90] = 1.53340e-2;
			currentOffPmos[100] =2.36007e-2;
			pnSizeRatio = (int)(currentOnNmos[0]/currentOnPmos[0]);
		}
	} else if (featureSizeInNano == 7) {
		if (deviceRoadmap == HP) {
			/* PTM model: 7nfet_HP.pm,7pfet_HP.pm, from http://ptm.asu.edu/ */
			vdd = 0.7;
			vth = 192.2339e-3;
			heightFin = 1.8e-8;
			widthFin = 7e-9;
			PitchFin = 2.2e-8;			
			phyGateLength = 1.1e-8;
			capIdealGate = 1.0487e-9;
			capFringe = 0;
			effectiveResistanceMultiplier = 1.45;	/* from CACTI */
			current_gmNmos=1.91E+03;
			current_gmPmos=8.02E+02;			
			currentOnNmos[0]   = 1912 ;
			currentOnNmos[10]  = 1937.6;
			currentOnNmos[20]  = 1962.6;
			currentOnNmos[30]  = 1987.1;
			currentOnNmos[40]  = 2011 ;
			currentOnNmos[50]  = 2034.4;
			currentOnNmos[60]  = 2057.2;
			currentOnNmos[70]  = 2079.5;
			currentOnNmos[80]  = 2101.3;
			currentOnNmos[90]  = 2122.6;
			currentOnNmos[100] = 2143.4;
			currentOnPmos[0]   = 1685.5; 
			currentOnPmos[10]  = 1716.4; 
			currentOnPmos[20]  = 1747.0; 
			currentOnPmos[30]  = 1777.4; 
			currentOnPmos[40]  = 1807.5; 
			currentOnPmos[50]  = 1837.4; 
			currentOnPmos[60]  = 1866.9; 
			currentOnPmos[70]  = 1896.1; 
			currentOnPmos[80]  = 1925.0; 
			currentOnPmos[90]  = 1953.6; 
			currentOnPmos[100] = 1981.8;
			currentOffNmos[0]  = 100.2258e-3;
			currentOffNmos[10] = 149.0252e-3;
			currentOffNmos[20] = 216.5654e-3;
			currentOffNmos[30] = 308.1967e-3;
			currentOffNmos[40] = 430.2635e-3;
			currentOffNmos[50] = 590.1731e-3;
			currentOffNmos[60] = 796.4489e-3;
			currentOffNmos[70] = 1.0588;
			currentOffNmos[80] = 1.3880;
			currentOffNmos[90] = 1.7960;
			currentOffNmos[100]= 2.2961;
			currentOffPmos[0]  = 97.9484e-3  ;
			currentOffPmos[10] = 156.3424e-3 ;
			currentOffPmos[20] = 243.1919e-3;
			currentOffPmos[30] = 369.4499e-3 ;
			currentOffPmos[40] = 549.1886e-3 ;
			currentOffPmos[50] = 800.1479e-3;
			currentOffPmos[60] = 1.1443      ;
			currentOffPmos[70] = 1.6083      ;
			currentOffPmos[80] = 2.2242      ;
			currentOffPmos[90] = 3.0295      ;
			currentOffPmos[100]= 4.0674      ;
			pnSizeRatio = (int)(currentOnNmos[0]/currentOnPmos[0]);
		} else {
			/* PTM model: 7nfet_LP.pm,7pfet_LP.pm, from http://ptm.asu.edu/ */
			vdd = 0.7;
			vth = 402.5252e-3;
			heightFin = 1.8e-8;
			widthFin = 7e-9;
			PitchFin = 2.2e-8;			
			phyGateLength = 1.1e-8;
			capIdealGate = 8.49489e-10;
			capFringe = 0;
			effectiveResistanceMultiplier = 1.73;	/* from CACTI */
			current_gmNmos=8.22E+02;
			current_gmPmos=1.45E+02;			
			currentOnNmos[0]  = 822.0573; 
			currentOnNmos[10] = 843.5584; 
			currentOnNmos[20] = 865.0229; 
			currentOnNmos[30] = 886.4385; 
			currentOnNmos[40] = 907.7931; 
			currentOnNmos[50] = 929.0751; 
			currentOnNmos[60] = 950.2729; 
			currentOnNmos[70] = 971.3751; 
			currentOnNmos[80] = 992.3706; 
			currentOnNmos[90] = 1.0132e3;
			currentOnNmos[100]= 1.0340e3; 
			currentOnPmos[0]  = 737.2425;
			currentOnPmos[10] = 763.7947;
			currentOnPmos[20] = 790.5774;
			currentOnPmos[30] = 817.5675;
			currentOnPmos[40] = 844.7417;
			currentOnPmos[50] = 872.0768;
			currentOnPmos[60] = 899.5498;
			currentOnPmos[70] = 927.1376;
			currentOnPmos[80] = 954.8176;
			currentOnPmos[90] = 982.567 ;
			currentOnPmos[100] =1.0104e3;
			currentOffNmos[0]  = 1.00E-04;
			currentOffNmos[10] = 1.85E-04;
			currentOffNmos[20] = 3.32E-04;
			currentOffNmos[30] = 5.74E-04;
			currentOffNmos[40] = 9.62E-04;
			currentOffNmos[50] = 1.5695e-3;
			currentOffNmos[60] = 2.4953e-3;
			currentOffNmos[70] = 3.8744e-3 ;
			currentOffNmos[80] = 5.8858e-3 ;
			currentOffNmos[90] = 8.7624e-3;
			currentOffNmos[100] =1.28025e-2;
			currentOffPmos[0]  = 100.9536e-3;
			currentOffPmos[10] = 201.3937e-3;
			currentOffPmos[20] = 386.2086e-3;
			currentOffPmos[30] = 714.4288e-3;
			currentOffPmos[40] = 1.2788e-3;
			currentOffPmos[50] = 2.2207e-3;
			currentOffPmos[60] = 3.7509e-3;
			currentOffPmos[70] = 6.1750e-3;
			currentOffPmos[80] = 9.9278e-3;
			currentOffPmos[90] = 1.56146e-2;
			currentOffPmos[100] =2.40633e-2;
			pnSizeRatio = (int)(currentOnNmos[0]/currentOnPmos[0]);
		}
	} else {
		cout<<"Error: CMOS Technology node "<< featureSizeInNano <<"nm is not supported"<<endl;
		exit(-1);
	}
	
	if (featureSizeInNano >= 22) {
		capOverlap = capIdealGate * 0.2;
	} else {
		capOverlap = 0;	// capOverlap and capFringe are included in capIdealGate in FinFET technology, so we let these two parameters 0
	}
	double cjd = 1e-3;			/* Bottom junction capacitance, Unit: F/m^2*/
	double cjswd = 2.5e-10;		/* Isolation-edge sidewall junction capacitance, Unit: F/m */
	double cjswgd = 0.5e-10;	/* Gate-edge sidewall junction capacitance, Unit: F/m */
	double mjd = 0.5;			/* Bottom junction capacitance grating coefficient */
	double mjswd = 0.33;		/* Isolation-edge sidewall junction capacitance grading coefficient */
	double mjswgd = 0.33;		/*   Gate-edge sidewall junction capacitance grading coefficient */
	buildInPotential = 0.9;		/* This value is from BSIM4 */
	capJunction = cjd / pow(1 + vdd / buildInPotential, mjd);
	capSidewall = cjswd / pow(1 + vdd / buildInPotential, mjswd);
	capDrainToChannel = cjswgd / pow(1 + vdd / buildInPotential, mjswgd);

	/* Properties not used so far */
	capPolywire = 0.0;	/* TO-DO: we need to find the values */

	/* Interpolate */
	for (int i = 1; i < 100; i++) {
		if (i % 10) {
			double a = currentOnNmos[i / 10 * 10];
			double b = currentOnNmos[i / 10 * 10 + 10];
			currentOnNmos[i] = a + (b-a) * (i % 10) / 10;
		}
	}

	for (int i = 1; i < 100; i++) {
		if (i % 10) {
			double a = currentOnPmos[i / 10 * 10];
			double b = currentOnPmos[i / 10 * 10 + 10];
			currentOnPmos[i] = a + (b-a) * (i % 10) / 10;
		}
	}

	for (int i = 1; i < 100; i++) {
		if (i % 10) {
			double a = currentOffNmos[i / 10 * 10];
			double b = currentOffNmos[i / 10 * 10 + 10];
			currentOffNmos[i] = a + (b-a) * (i % 10) / 10;
		}
	}

	for (int i = 1; i < 100; i++) {
		if (i % 10) {
			double a = currentOffPmos[i / 10 * 10];
			double b = currentOffPmos[i / 10 * 10 + 10];
			currentOffPmos[i] = a + (b-a) * (i % 10) / 10;
		}
	}

	initialized = true;
}
void Technology::PrintProperty() {
	// TODO
}
