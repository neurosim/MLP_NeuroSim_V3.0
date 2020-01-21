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

#ifndef ARRAY_H_
#define ARRAY_H_

#include <cstdlib>
#include "Cell.h"

class Array {
public:
	Cell ***cell;
	int arrayColSize, arrayRowSize, wireWidth;
	double unitLengthWireResistance;
	double wireResistanceRow, wireResistanceCol;
	double wireCapRow;	// Cap of the WL (cross-point) or BL (1T1R)
	double wireCapCol;	// Cap of the BL (cross-point) or SL (1T1R)
	double wireGateCapRow;	// Cap of 1T1R WL cap
	double wireCapBLCol;	// Cap of 1T1R BL cap in digital eNVM
	double readEnergy, writeEnergy;
    double transferReadEnergy, transferWriteEnergy, transferEnergy; // the transfer energy for the hybrid cell
	int numCellPerSynapse;	// For SRAM to use redundant cells to represent one synapse
	double writeEnergySRAMCell;	// Write energy per SRAM cell (will move this to SRAM cell level in the future)
	bool **weightChange;	// Specify if the weight value will change or not during weight update (for SRAM and digital eNVM)
    int refColumnNumber;
	/* Constructor */
    // code modified
	Array(int arrayColSize, int arrayRowSize, int wireWidth) {  
        this->arrayRowSize = arrayRowSize;
        this->arrayColSize = arrayColSize;
        this->wireWidth = wireWidth;
		readEnergy = 0;
		writeEnergy = 0;
        transferReadEnergy = transferWriteEnergy = 0;
        transferEnergy = 0;

		/* Initialize weightChange */
		weightChange = new bool*[arrayColSize];
		for (int col=0; col<arrayColSize; col++) {
			weightChange[col] = new bool[arrayRowSize];
		}
	}

	template <class memoryType>
	void Initialization(int numCellPerSynapse=1,bool refColumn = false) { // default value is 1
		/* Determine number of cells per synapse (SRAM and DigitalNVM) */
		this->numCellPerSynapse = numCellPerSynapse;

		/* Initialize memory cells */
        int cellsPerRow; // the number of columns 
        if(refColumn == true) 
            cellsPerRow = arrayColSize*numCellPerSynapse+2;
        else
            cellsPerRow = arrayColSize*numCellPerSynapse;
        cell = new Cell**[cellsPerRow];
		for (int col=0; col<cellsPerRow; col++) {
			cell[col] = new Cell*[arrayRowSize];
			for (int row=0; row<arrayRowSize; row++) {
				cell[col][row] = new memoryType(col, row);
			}
		}
        // initialize the conductance of the reference column
        if(refColumn == true)
        {
            refColumnNumber = arrayColSize*numCellPerSynapse; // the column number of the first reference column
            if(DigitalNVM *temp = dynamic_cast<DigitalNVM*>(this->cell[0][0]))
            {
                for(int row=0; row < this-> arrayRowSize; row++)
                {
                    static_cast<DigitalNVM*>(this->cell[refColumnNumber][row])->conductance = static_cast<DigitalNVM*>(this->cell[refColumnNumber][row])->minConductance;
                    static_cast<DigitalNVM*>(this->cell[refColumnNumber+1][row])->conductance = static_cast<DigitalNVM*>(this->cell[refColumnNumber+1][row])->maxConductance;
                }
            }    
        }
		
		/* Initialize interconnect wires */
		double AR;	// Aspect ratio of wire height to wire width
		double Rho;	// Resistivity
		switch(wireWidth) {// wireWidth is defined in the param.cpp
      case 800: AR = 2.10; Rho= 2.42e-8; break;
			case 200: AR = 2.10; Rho = 2.42e-8; break;
			case 100:	AR = 2.30; Rho = 2.73e-8; break;
			case 50:	AR = 2.34; Rho = 3.91e-8; break;
			case 40:	AR = 1.90; Rho = 4.03e-8; break;
			case 32:	AR = 1.90; Rho = 4.51e-8; break;
			case 22:	AR = 2.00; Rho = 5.41e-8; break;
			case 14:	AR = 2.10; Rho = 7.43e-8; break;
			case -1:	break;	// Ignore wire resistance or user define
			default:	exit(-1); puts("Wire width out of range"); 
		}
		double wireLength = wireWidth * 1e-9 * 2;	// 2F
		if (wireWidth == -1) {
			unitLengthWireResistance = 1.0;	// Use a small number to prevent numerical error for NeuroSim
			wireResistanceRow = 0;
			wireResistanceCol = 0;
		} else {
			unitLengthWireResistance =  Rho / ( wireWidth*1e-9 * wireWidth*1e-9 * AR );
			wireResistanceRow = unitLengthWireResistance * wireLength;
			wireResistanceCol = unitLengthWireResistance * wireLength;
		}
		wireCapRow = wireLength * 0.2e-15/1e-6;
		wireCapCol = wireLength * 0.2e-15/1e-6;
		wireGateCapRow = wireLength * 0.2e-15/1e-6;
		
	}

	double ReadCell(int x, int y,char*mode=NULL);	// x (column) and y (row) start from index 0
	void WriteCell(int x, int y, double deltaWeight, double weight, double maxWeight, double minWeight, bool regular);
	double GetMaxCellReadCurrent(int x, int y, char*mode=NULL);
	double GetMinCellReadCurrent(int x, int y, char*mode=NULL);
	double GetMediumCellReadCurrent(int x, int y);
	double ConductanceToWeight(int x, int y, double maxWeight, double minWeight,char* mode=NULL);
};

#endif
