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

#include "formula.h"
#include "Array.h"

int counter=0;
double Array::ReadCell(int x, int y, char* mode) {
    // mode is only for the 3T1C cell to select LSB or MSB
    // it should be "MSB_LTP","MSB_LTD" or "LSB" 
	if (AnalogNVM *temp = dynamic_cast<AnalogNVM*>(**cell)){ // Analog eNVM
		double readVoltage = static_cast<eNVM*>(cell[x][y])->readVoltage;
		double totalWireResistance;
		if (static_cast<eNVM*>(cell[x][y])->cmosAccess){  // 1T1R cell or 1T1C cell
			if (static_cast<AnalogNVM*>(cell[x][y])->FeFET) // FeFET
				totalWireResistance = (x + 1) * wireResistanceRow + (arrayRowSize - y) * wireResistanceCol; // do not need to consider the access resistance
            else // Normal
				totalWireResistance = (x + 1) * wireResistanceRow + (arrayRowSize - y) * wireResistanceCol + static_cast<eNVM*>(cell[x][y])->resistanceAccess;
		} 
        else 
			totalWireResistance = (x + 1) * wireResistanceRow + (arrayRowSize - y) * wireResistanceCol;
		double cellCurrent;
		if (static_cast<eNVM*>(cell[x][y])->nonlinearIV){
			// Bisection method to calculate read current with nonlinearity
			int maxIter = 30;
			double v1 = 0, v2 = readVoltage, v3;
			double wireCurrent;
			for (int iter=0; iter<maxIter; iter++){
				v3 = (v1 + v2)/2;
				wireCurrent = (readVoltage - v3)/totalWireResistance;
				cellCurrent = static_cast<AnalogNVM*>(cell[x][y])->Read(v3);
				if (wireCurrent > cellCurrent)
					v1 = v3;
				else
					v2 = v3;
			}
		} 
        else{	// No nonlinearity
			if (static_cast<eNVM*>(cell[x][y])->readNoise){
				extern std::mt19937 gen;
				cellCurrent = readVoltage / (1/static_cast<eNVM*>(cell[x][y])->conductance * (1 + (*static_cast<eNVM*>(cell[x][y])->gaussian_dist)(gen)) + totalWireResistance);
			} 
            else
				cellCurrent = readVoltage / (1/static_cast<eNVM*>(cell[x][y])->conductance + totalWireResistance);
		}
		return cellCurrent;
	} 
    else if (HybridCell *temp = dynamic_cast<HybridCell*>(**cell)){
        if(mode=="LSB"){
            extern std::mt19937 gen;
            double readVoltage_LSB =  static_cast<HybridCell*>(cell[x][y])->LSBcell.readVoltage;
            double totalWireResistance_LSB = (x + 1) * wireResistanceRow + (arrayRowSize - y) * wireResistanceCol;
            double cellCurrent_LSB;
            if (static_cast<HybridCell*>(cell[x][y])->LSBcell.readNoise)
				cellCurrent_LSB = readVoltage_LSB / (1/static_cast<HybridCell*>(cell[x][y])->LSBcell.conductance * (1 + (*static_cast<HybridCell*>(cell[x][y])->LSBcell.gaussian_dist)(gen)) + totalWireResistance_LSB);
            else
                cellCurrent_LSB = readVoltage_LSB / (1/static_cast<HybridCell*>(cell[x][y])->LSBcell.conductance + totalWireResistance_LSB);      
            return cellCurrent_LSB;
        }
        else if(mode=="MSB_LTP"){
            extern std::mt19937 gen;
            double readVoltage_MSB = static_cast<HybridCell*>(cell[x][y])->MSBcell_LTP.readVoltage;
            double totalWireResistance_MSB=  (x + 1) * wireResistanceRow + (arrayRowSize - y) * wireResistanceCol +static_cast<HybridCell*>(cell[x][y])->MSBcell_LTP.resistanceAccess;
            double cellCurrent_MSB_LTP;
            if (static_cast<HybridCell*>(cell[x][y])->MSBcell_LTP.readNoise) 
                cellCurrent_MSB_LTP = readVoltage_MSB / (1/static_cast<HybridCell*>(cell[x][y])->MSBcell_LTP.conductance * (1 + (*static_cast<HybridCell*>(cell[x][y])->MSBcell_LTP.gaussian_dist)(gen)) + totalWireResistance_MSB);
            else
                cellCurrent_MSB_LTP = readVoltage_MSB / (1/static_cast<HybridCell*>(cell[x][y])->MSBcell_LTP.conductance +  totalWireResistance_MSB);
           return cellCurrent_MSB_LTP;
        }
        else if(mode=="MSB_LTD"){
            extern std::mt19937 gen;
            double readVoltage_MSB = static_cast<HybridCell*>(cell[x][y])->MSBcell_LTD.readVoltage;  
            double totalWireResistance_MSB=  (x + 1) * wireResistanceRow + (arrayRowSize - y) * wireResistanceCol +(static_cast<HybridCell*>(cell[x][y])->MSBcell_LTP).resistanceAccess;
            double cellCurrent_MSB_LTD;          
            if (static_cast<HybridCell*>(cell[x][y])->MSBcell_LTD.readNoise) 
                cellCurrent_MSB_LTD = readVoltage_MSB / (1/static_cast<HybridCell*>(cell[x][y])->MSBcell_LTD.conductance * (1 + (*static_cast<HybridCell*>(cell[x][y])->MSBcell_LTD.gaussian_dist)(gen)) + totalWireResistance_MSB);
            else
                cellCurrent_MSB_LTD = readVoltage_MSB / (1/static_cast<HybridCell*>(cell[x][y])->MSBcell_LTD.conductance + totalWireResistance_MSB); 
            return cellCurrent_MSB_LTD;  
        }
        else printf("Please specify the correct reading mode\n");
    }
    else{ // SRAM or digital eNVM
		int weightDigits = 0;
		if (DigitalNVM *temp = dynamic_cast<DigitalNVM*>(**cell)) {	// Digital eNVM
			for (int n=0; n<numCellPerSynapse; n++){   // n=0 is LSB
				int colIndex = (x+1) * numCellPerSynapse - (n+1);
				double readVoltage = static_cast<eNVM*>(cell[colIndex][y])->readVoltage;
				double totalWireResistance;
				if (static_cast<eNVM*>(cell[colIndex][y])->cmosAccess) 
					totalWireResistance = (colIndex + 1) * wireResistanceRow + (arrayRowSize - y) * wireResistanceCol + static_cast<eNVM*>(cell[colIndex][y])->resistanceAccess; 
                else 
					totalWireResistance = (colIndex + 1) * wireResistanceRow + (arrayRowSize - y) * wireResistanceCol;
				double cellCurrent;
				if (static_cast<eNVM*>(cell[colIndex][y])->nonlinearIV) {
					/* Bisection method to calculate read current with nonlinearity */
					int maxIter = 30;
					double v1 = 0, v2 = readVoltage, v3;
					double wireCurrent;
					for (int iter=0; iter<maxIter; iter++){
						v3 = (v1 + v2)/2;
						wireCurrent = (readVoltage - v3)/totalWireResistance;
						cellCurrent = static_cast<DigitalNVM*>(cell[colIndex][y])->Read(v3);
						if (wireCurrent > cellCurrent)
							v1 = v3;
						else
							v2 = v3;
					}
				} 
                else{ // No nonlinearity 
					if (static_cast<eNVM*>(cell[colIndex][y])->readNoise){
						extern std::mt19937 gen;
						cellCurrent = readVoltage / (1/static_cast<eNVM*>(cell[colIndex][y])->conductance * (1 + (*static_cast<eNVM*>(cell[colIndex][y])->gaussian_dist)(gen)) + totalWireResistance);
					} 
                    else 
						cellCurrent = readVoltage / (1/static_cast<eNVM*>(cell[colIndex][y])->conductance + totalWireResistance);
				}
				// Current sensing
				int bit;
				if (cellCurrent >= static_cast<DigitalNVM*>(cell[colIndex][y])->refCurrent) 
					bit = 1;
                else 
					bit = 0;
				weightDigits += bit * pow(2, n);	// If the rightmost is LSB
			}
		} 
        else{	// SRAM
			for (int n=0; n<numCellPerSynapse; n++){ // n=0 is LSB
				weightDigits += static_cast<SRAM*>(cell[(x+1) * numCellPerSynapse - (n+1)][y])->bit * pow(2, n);    // If the rightmost is LSB
			}
		}
		return weightDigits;
	}
}

void Array::WriteCell(int x, int y, double deltaWeight, double weight, double maxWeight, double minWeight, 
						bool regular /* False: ideal write, True: regular write considering device properties */){
	// TODO: include wire resistance
	if (AnalogNVM *temp = dynamic_cast<AnalogNVM*>(**cell)){ // Analog eNVM
        if (regular)	// Regular write
			static_cast<AnalogNVM*>(cell[x][y])->Write(deltaWeight, weight, minWeight, maxWeight);
        else{	
			double conductance = 0;
			double maxConductance = static_cast<eNVM*>(cell[x][y])->maxConductance;
			double minConductance = static_cast<eNVM*>(cell[x][y])->minConductance;
			conductance = (weight-minWeight)/(maxWeight-minWeight) * (maxConductance - minConductance);
			if (conductance > maxConductance) 
				conductance = maxConductance; 
            else if (conductance < minConductance) 
				conductance = minConductance;
			static_cast<eNVM*>(cell[x][y])->conductance = conductance;
		}
	}
    else if(HybridCell*temp = dynamic_cast<HybridCell*>(**cell)){
        double weightLSB = this->ConductanceToWeight(x,y, maxWeight, minWeight, "LSB");
		
        if (regular) // Regular write
            static_cast<HybridCell*>(cell[x][y])->Write((static_cast<HybridCell*>(cell[x][y])->significance+1)*deltaWeight, weightLSB, minWeight, maxWeight);
            //static_cast<HybridCell*>(cell[x][y])->Write(deltaWeight, weightLSB, minWeight, maxWeight);
        else{
			double conductance = 0;
			double maxConductance = static_cast<HybridCell*>(cell[x][y])->LSBcell.maxConductance;
			double minConductance = static_cast<HybridCell*>(cell[x][y])->LSBcell.minConductance;
			//weight = weight/static_cast<HybridCell*>(cell[x][y])->significance;
			conductance = (weight-minWeight)/(maxWeight-minWeight) * (maxConductance - minConductance)+minConductance;
			if (conductance > maxConductance) 
				conductance = maxConductance;
            else if (conductance < minConductance) 
				conductance = minConductance;
			static_cast<HybridCell*>(cell[x][y])->LSBcell.conductance = conductance;            
        }    
    }
    else{    // SRAM or digital eNVM 
		// firstly need to truncate weight(-1, +1) to weight(0, 1), then truncate to weight(0, numLevel)
		int numLevel = pow(2, numCellPerSynapse);
		weightChange[x][y] = (deltaWeight != 0)? true : false; // only do update for the cells with Delta weight !=0
		int targetWeightDigits = (int)((weight + 1)/2 * (numLevel-1));// mapping (-1,+1) to (0,1), the number of conductance levels that need to update
		
        int maxWeightDigits = pow(2, numCellPerSynapse) - 1;
		if (targetWeightDigits > maxWeightDigits)
			targetWeightDigits = maxWeightDigits;
		else if (targetWeightDigits < 0)
			targetWeightDigits = 0;		
		/* Write new weight and calculate write energy */
		if (DigitalNVM *temp = dynamic_cast<DigitalNVM*>(**cell)){ // Digital eNVM
			for (int n=0; n<numCellPerSynapse; n++){ // n=0 is LSB
				int bitNew = ((targetWeightDigits >> n) & 1); //get the nth bit to write to
				/* Write new weight */
				if (static_cast<eNVM*>(cell[x][y])->cmosAccess) // 1T1R
					static_cast<DigitalNVM*>(cell[(x+1) * numCellPerSynapse - (n+1)][y])->Write(bitNew, wireCapBLCol);
                else // Cross-point
					static_cast<DigitalNVM*>(cell[(x+1) * numCellPerSynapse - (n+1)][y])->Write(bitNew, wireCapCol);
			}
		} 
        else {
			static_cast<SRAM*>(cell[x * numCellPerSynapse][y])->writeEnergy = 0;    // Use the MSB cell to store the info of the write energy of the synapse
			for (int n=0; n<numCellPerSynapse; n++){   // n=0 is LSB
				int bit = static_cast<SRAM*>(cell[(x+1) * numCellPerSynapse - (n+1)][y])->bit;
				int bitNew = ((targetWeightDigits >> n) & 1);
				if (bit != bitNew) // Consume write energy if the new bit is different than the current bit
					static_cast<SRAM*>(cell[x * numCellPerSynapse][y])->writeEnergy += writeEnergySRAMCell; // Currently this writeEnergySRAMCell is the array level parameter
				/* Write new weight */
				static_cast<SRAM*>(cell[(x+1) * numCellPerSynapse - (n+1)][y])->bitPrev = bit;	// If the rightmost is LSB
				static_cast<SRAM*>(cell[(x+1) * numCellPerSynapse - (n+1)][y])->bit = bitNew;	// If the rightmost is LSB
			}
		}
	}
}

double Array::GetMaxCellReadCurrent(int x, int y, char* mode) { 
    // two mode: "LSB", "MSB". For hybrid cell only
    if(AnalogNVM*temp = dynamic_cast<AnalogNVM*>(**cell)) 
	    return static_cast<AnalogNVM*>(cell[x][y])->GetMaxReadCurrent();
    else if (HybridCell* temp = dynamic_cast<HybridCell*>(**cell)){   
        if(mode=="LSB")
            return static_cast<HybridCell*>(cell[x][y])->LSBcell.GetMaxReadCurrent();
        else if(mode =="MSB")
            return static_cast<HybridCell*>(cell[x][y])->MSBcell_LTP.GetMaxReadCurrent();
        else 
            printf("Please specify the correct reading mode\n");            
    }
}

double Array::GetMinCellReadCurrent(int x, int y, char*mode) {
    // two mode: "LSB", "MSB". For hybrid cell only
    if(AnalogNVM*temp = dynamic_cast<AnalogNVM*>(**cell)) 
	    return static_cast<AnalogNVM*>(cell[x][y])->GetMinReadCurrent();
    else if (HybridCell* temp = dynamic_cast<HybridCell*>(**cell)){   
        if(mode=="LSB")
            return static_cast<HybridCell*>(cell[x][y])->LSBcell.GetMinReadCurrent();
        else if(mode =="MSB")
            return static_cast<HybridCell*>(cell[x][y])->MSBcell_LTP.GetMinReadCurrent();
        else 
            printf("Please specify the correct reading mode\n");            
    }
}

double Array::GetMediumCellReadCurrent(int x, int y) {  
    double Imax, Imin;
    if(AnalogNVM*temp = dynamic_cast<AnalogNVM*>(**cell)){
	     Imax = static_cast<AnalogNVM*>(cell[x][y])->GetMaxReadCurrent();
         Imin = static_cast<AnalogNVM*>(cell[x][y])->GetMinReadCurrent();
    }
    else if(HybridCell*temp = dynamic_cast<HybridCell*>(**cell)){
             Imax = static_cast<HybridCell*>(cell[x][y])->LSBcell.GetMaxReadCurrent();
	         Imin = static_cast<HybridCell*>(cell[x][y])->LSBcell.GetMinReadCurrent();
    }
	return (Imax+Imin)/2;
}

// convert the conductance to -1~1 
double Array::ConductanceToWeight(int x, int y, double maxWeight, double minWeight, char* mode) {
	if (AnalogNVM *temp = dynamic_cast<AnalogNVM*>(**cell)){	// Analog eNVM
		/* Measure current */
		double I = this->ReadCell(x, y); // for AnalogNVM, read the current and convert it into conductance
		/* Convert current to weight */
		double Imax = static_cast<AnalogNVM*>(cell[x][y])->GetMaxReadCurrent(); // the current when Conductance is the minimum
		double Imin = static_cast<AnalogNVM*>(cell[x][y])->GetMinReadCurrent(); // the current when Conductance is the maximum
		if (I<Imin)
			I = Imin;
		else if (I>Imax)
			I = Imax;
		return (I-Imin) / (Imax-Imin) * (maxWeight-minWeight) + minWeight;
	}
    else if (HybridCell *temp = dynamic_cast<HybridCell*>(**cell)){
		double I = this->ReadCell(x, y,"LSB"); // for 3T1C cell, read the current and convert it into conductance
		double Imax = static_cast<HybridCell*>(cell[x][y])->LSBcell.GetMaxReadCurrent(); // the current when Conductance is the minimum
		double Imin = static_cast<HybridCell*>(cell[x][y])->LSBcell.GetMinReadCurrent(); // the current when Conductance is the maximum
		if (I<Imin)
			I = Imin;
		else if (I>Imax)
			I = Imax;
		double weightLSB= (I-Imin) / (Imax-Imin) * (maxWeight-minWeight) + minWeight; 
        if(mode == "LSB")
            return weightLSB;
        else{
            double I_LTP = this->ReadCell(x,y,"MSB_LTP");
            double I_LTD = this->ReadCell(x,y,"MSB_LTD");
            Imax = static_cast<HybridCell*>(cell[x][y])->MSBcell_LTP.GetMaxReadCurrent();
            Imin = static_cast<HybridCell*>(cell[x][y])->MSBcell_LTP.GetMinReadCurrent();
            if(I_LTP>Imax)
                I_LTP=Imax;
            else if(I_LTP<Imin)
                I_LTP <Imin;
                
            if(I_LTD>Imax)
                I_LTD=Imax;
            else if(I_LTD<Imin)
                I_LTD <Imin;

            double weightMSB =     (I_LTP-I_LTD)/(Imax-Imin)*(maxWeight-0);
            double weightMSB_LTP = (I_LTP-Imin)/(Imax-Imin)*(maxWeight-0)+0;
            double weightMSB_LTD = (I_LTD-Imin)/(Imax-Imin)*(maxWeight-0)+0;            
            double weight = (static_cast<HybridCell*>(cell[x][y])->significance*weightMSB+weightLSB)/(static_cast<HybridCell*>(cell[x][y])->significance+1);
            if(weight > maxWeight)
                weight = maxWeight;
            else if (weight < minWeight)
                weight = minWeight;
            if(mode=="MSB")
                return weightMSB;
            else if(mode=="MSB_LTP")
                return weightMSB_LTP;
            else if(mode=="MSB_LTD")
                return weightMSB_LTD;
            else // return the weight of the cell is not specified;
                return weight;
        }
    }
    else{	// SRAM or digital eNVM
		double weightDigits = this->ReadCell(x, y);
		int weightDigitsMax = pow(2, numCellPerSynapse) - 1;
		return (weightDigits / weightDigitsMax) * (maxWeight - minWeight) + minWeight;
	}
}

