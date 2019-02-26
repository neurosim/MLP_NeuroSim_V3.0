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

#ifndef _MEMCELL_H_
#define _MEMCELL_H_

#include "typedef.h"

class MemCell {
public:
	/* Properties */
	Type::MemCellType memCellType;	/* Memory cell type (like MRAM, PCRAM, etc.) */
	int processNode;        /* Cell original process technology node, Unit: nm*/
	double area;			/* Cell area, Unit: F^2 */
	double aspectRatio;		/* Cell aspect ratio, H/W */
	double widthInFeatureSize;	/* Cell width, Unit: F */
	double heightInFeatureSize;	/* Cell height, Unit: F */
    double widthAccessTransistor; /*width of the access transistor*/
	double resistanceOn;	/* Turn-on resistance */
	double resistanceOff;	/* Turn-off resistance */
	double minSenseVoltage; /* Minimum sense voltage */
	
	CellAccessType accessType;	/* Cell access type: CMOS, BJT, or diode */
	double featureSize;
	double accessVoltage;
	double readVoltage;
	double writeVoltage;
	double readPulseWidth;
	double writePulseWidth;
	bool nonlinearIV;	/* Consider I-V nonlinearity or not (Currently this option is for cross-point array. It is hard to have this option in pseudo-crossbar since it has an access transistor and the transistor's resistance can be comparable to RRAM's resistance after considering the nonlinearity. In this case, we have to iteratively find both the resistance and Vw across RRAM.) */
	double nonlinearity;	/* Current at write voltage / current at 1/2 write voltage */
	double resistanceAvg;
	double resCellAccess;
	double resMemCellOn;	// At on-chip Vr (different than the Vr in the reported measurement data)
	double resMemCellOff;	// At on-chip Vr (different than the Vr in the reported measurement data)
	double resMemCellAvg;	// At on-chip Vr (different than the Vr in the reported measurement data)
	double resMemCellOnAtHalfVw;
	double resMemCellOffAtHalfVw;
	double resMemCellAvgAtHalfVw;
	double resMemCellOnAtVw;
	double resMemCellOffAtVw;
	double resMemCellAvgAtVw;
	double capSRAMCell;
	int multipleCells;	/* Use multiple cells as one weight element to reduce the variation (only layout now) */

	/* Optional properties */
	double widthAccessCMOS;	/* The gate width of CMOS access transistor, Unit: F */
	double widthSRAMCellNMOS;	/* The gate width of NMOS in SRAM cells, Unit: F */
	double widthSRAMCellPMOS;	/* The gate width of PMOS in SRAM cells, Unit: F */
    
    double widthAccessNMOS;  // for the 2T1F,3T1C cell where access transistor is needed. 
    double widthAccessPMOS; 
    double widthStorageNode; // the width of the transistor as the storage node
};

#endif /* MEMCELL_H_ */
