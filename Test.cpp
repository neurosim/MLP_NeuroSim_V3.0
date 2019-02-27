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

#include <cstdio>
#include <iostream>
#include <vector>
#include <random>
#include "formula.h"
#include "Param.h"
#include "Array.h"
#include "Mapping.h"
#include "NeuroSim.h"
#include "Cell.h"

extern Param *param;

extern std::vector< std::vector<double> > testInput;
extern std::vector< std::vector<int> > dTestInput;
extern std::vector< std::vector<double> > testOutput;

extern std::vector< std::vector<double> > weight1;
extern std::vector< std::vector<double> > weight2;

extern Technology techIH;
extern Technology techHO;
extern Array *arrayIH;
extern Array *arrayHO;
extern SubArray *subArrayIH;
extern SubArray *subArrayHO;
extern Adder adderIH;
extern Mux muxIH;
extern RowDecoder muxDecoderIH;
extern DFF dffIH;
extern Subtractor subtractorIH;
extern Adder adderHO;
extern Mux muxHO;
extern RowDecoder muxDecoderHO;
extern DFF dffHO;
extern Subtractor subtractorHO;

extern int correct;		// # of correct prediction

/* Validation */
void Validate() {
	int numBatchReadSynapse;    // # of read synapses in a batch read operation (decide later)
	double outN1[param->nHide]; // Net input to the hidden layer [param->nHide]
	double a1[param->nHide];    // Net output of hidden layer [param->nHide] also the input of hidden layer to output layer
	int da1[param->nHide];  // Digitized net output of hidden layer [param->nHide] also the input of hidden layer to output layer
	double outN2[param->nOutput];   // Net input to the output layer [param->nOutput]
	double a2[param->nOutput];  // Net output of output layer [param->nOutput]
	double tempMax;
	int countNum;
	correct = 0;

	double sumArrayReadEnergyIH = 0;   // Use a temporary variable here since OpenMP does not support reduction on class member
	double sumNeuroSimReadEnergyIH = 0;   // Use a temporary variable here since OpenMP does not support reduction on class member
	double sumReadLatencyIH = 0;    // Use a temporary variable here since OpenMP does not support reduction on class member
    double readVoltageIH,readVoltageHO,readVoltageMSB;
    double readPulseWidthIH,readPulseWidthHO,readPulseWidthMSB;
	double sumArrayReadEnergyHO = 0;    // Use a temporary variable here since OpenMP does not support reduction on class member
	double sumNeuroSimReadEnergyHO = 0; // Use a temporary variable here since OpenMP does not support reduction on class member
	double sumReadLatencyHO = 0;    // Use a temporary variable here since OpenMP does not support reduction on class member
    if(eNVM* temp = dynamic_cast<eNVM*>(arrayIH->cell[0][0]))
    {
        readVoltageIH = static_cast<eNVM*>(arrayIH->cell[0][0])->readVoltage;
        readVoltageHO = static_cast<eNVM*>(arrayHO->cell[0][0])->readVoltage;
        readPulseWidthIH = static_cast<eNVM*>(arrayIH->cell[0][0])->readPulseWidth;
	    readPulseWidthHO = static_cast<eNVM*>(arrayHO->cell[0][0])->readPulseWidth;
    }
    
    #pragma omp parallel for private(outN1, a1, da1, outN2, a2, tempMax, countNum, numBatchReadSynapse) reduction(+: correct, sumArrayReadEnergyIH, sumNeuroSimReadEnergyIH, sumArrayReadEnergyHO, sumNeuroSimReadEnergyHO, sumReadLatencyIH, sumReadLatencyHO)
	for (int i = 0; i < param->numMnistTestImages; i++)
	{
		// Forward propagation
		/* First layer from input layer to the hidden layer */
		std::fill_n(outN1, param->nHide, 0);
		std::fill_n(a1, param->nHide, 0);
		if (param->useHardwareInTestingFF) {    // Hardware
			for (int j=0; j<param->nHide; j++) {
				if (AnalogNVM *temp = dynamic_cast<AnalogNVM*>(arrayIH->cell[0][0])) {  // Analog eNVM
					if (static_cast<eNVM*>(arrayIH->cell[0][0])->cmosAccess) {  // 1T1R
						sumArrayReadEnergyIH += arrayIH->wireGateCapRow * techIH.vdd * techIH.vdd * param->nInput; // All WLs open
					}
				} else if (DigitalNVM *temp = dynamic_cast<DigitalNVM*>(arrayIH->cell[0][0])) { // Digital eNVM
					if (static_cast<eNVM*>(arrayIH->cell[0][0])->cmosAccess) {  // 1T1R
						sumArrayReadEnergyIH += arrayIH->wireGateCapRow * techIH.vdd * techIH.vdd;  // Selected WL
					} else {    // Cross-point
						sumArrayReadEnergyIH += arrayIH->wireCapRow * techIH.vdd * techIH.vdd * (param->nInput - 1);    // Unselected WLs
					}
				}
				
                for (int n=0; n<param->numBitInput; n++) {
					double pSumMaxAlgorithm = pow(2, n) / (param->numInputLevel - 1) * arrayIH->arrayRowSize;   // Max algorithm partial weighted sum for the nth vector bit (if both max input value and max weight are 1)
					if (AnalogNVM *temp = dynamic_cast<AnalogNVM*>(arrayIH->cell[0][0])) {  // Analog eNVM
						double Isum = 0;    // weighted sum current
						double IsumMax = 0; // Max weighted sum current
						double IsumMin = 0; // Max weighted sum current
						double inputSum = 0;    // Weighted sum current of input vector * weight=1 column
						for (int k=0; k<param->nInput; k++) {
							if ((dTestInput[i][k]>>n) & 1) {    // if the nth bit of dTestInput[i][k] is 1
								Isum += arrayIH->ReadCell(j,k);
								inputSum += arrayIH->GetMediumCellReadCurrent(j,k);
								sumArrayReadEnergyIH += arrayIH->wireCapRow * readVoltageIH * readVoltageIH;   // Selected BLs (1T1R) or Selected WLs (cross-point)
							}
							IsumMax += arrayIH->GetMaxCellReadCurrent(j,k);
							IsumMin += arrayIH->GetMinCellReadCurrent(j,k);
						}
						sumArrayReadEnergyIH += Isum * readVoltageIH * readPulseWidthIH;
						int outputDigits = 2*(CurrentToDigits(Isum, IsumMax-IsumMin)-CurrentToDigits(inputSum, IsumMax-IsumMin));
						outN1[j] += DigitsToAlgorithm(outputDigits, pSumMaxAlgorithm);
					} 
                    else {
                            bool digitalNVM = false; 
                            bool parallelRead = false;
                            if(DigitalNVM*temp = dynamic_cast<DigitalNVM*>(arrayIH->cell[0][0]))
                            {    digitalNVM = true;
                                if(static_cast<DigitalNVM*>(arrayIH->cell[0][0])->parallelRead == true) 
								{
                                    parallelRead = true;
                                }
                            }
                            if(digitalNVM && parallelRead) // parallel read-out for DigitalNVM
                            {
                                    double Imax = static_cast<DigitalNVM*>(arrayIH->cell[0][0])->avgMaxConductance*static_cast<DigitalNVM*>(arrayIH->cell[0][0])->readVoltage;
                                    double Imin = static_cast<DigitalNVM*>(arrayIH->cell[0][0])->avgMinConductance*static_cast<DigitalNVM*>(arrayIH->cell[0][0])->readVoltage;
                                    double Isum = 0;    // weighted sum current
							        double IsumMax = 0; // Max weighted sum current
							        double inputSum = 0;    // Weighted sum current of input vector * weight=1 column
                                    int Dsum=0;
                                    int DsumMax = 0;
                                    int Dref = 0;
                                    for (int w=0;w<param->numWeightBit;w++){
                                        int colIndex = (j+1) * param->numWeightBit - (w+1);  // w=0 is the LSB
									    for (int k=0; k<param->nInput; k++) 
                                        {
										    if((dTestInput[i][k]>>n) & 1){ // accumulate the current along a column
											    Isum += static_cast<DigitalNVM*>(arrayIH->cell[colIndex ][k])->conductance*static_cast<DigitalNVM*>(arrayIH->cell[colIndex ][k])->readVoltage;
                                                inputSum += static_cast<DigitalNVM*>(arrayIH->cell[arrayIH->refColumnNumber][k])->conductance*static_cast<DigitalNVM*>(arrayIH->cell[arrayIH->refColumnNumber][k])->readVoltage;
										    }
									    }
                                        int outputDigits = (int) (Isum /(Imax-Imin)); // the output at the ADC of this column
                                        int outputDigitsRef = (int) (inputSum/(Imax-Imin));
                                        outputDigits = outputDigits-outputDigitsRef;
                                            
                                        Dref = (int)(inputSum/Imin);
                                        Isum=0;
                                        inputSum=0;
                                        Dsum += outputDigits*(int) pow(2,w);  // get the weight represented by the column
                                        DsumMax += param->nInput*(int) pow(2,w); // the maximum weight that can be represented by this column
        
                                    }
                                    outN1[j] += (double)(Dsum - Dref*(pow(2,param->numWeightBit-1)-1)) / DsumMax * pSumMaxAlgorithm;
                                    sumArrayReadEnergyIH  += static_cast<DigitalNVM*>(arrayIH->cell[0][0])->readEnergy * arrayIH->numCellPerSynapse * arrayIH->arrayRowSize;
                            }
                            else
                            {	 // Digital NVM or SRAM row-by-row readout				
							    int Dsum = 0;
							    int DsumMax = 0;
							    int inputSum = 0;
							    for (int k=0; k<param->nInput; k++) {
								    if ((dTestInput[i][k]>>n) & 1) {    // if the nth bit of dInput[i][k] is 1
									    Dsum += (int)(arrayIH->ReadCell(j,k));
									    inputSum += pow(2, arrayIH->numCellPerSynapse-1) - 1;   // get the digital weights of the dummy column as reference
								    }
								    DsumMax += pow(2, arrayIH->numCellPerSynapse) - 1;
							    }
							    if (DigitalNVM *temp = dynamic_cast<DigitalNVM*>(arrayIH->cell[0][0])) {    // Digital eNVM
								    sumArrayReadEnergyIH  += static_cast<DigitalNVM*>(arrayIH->cell[0][0])->readEnergy * arrayIH->numCellPerSynapse * arrayIH->arrayRowSize;
							    } 
                                else {    // SRAM
								    sumArrayReadEnergyIH  += static_cast<SRAM*>(arrayIH->cell[0][0])->readEnergy * arrayIH->numCellPerSynapse * arrayIH->arrayRowSize;
							    }
							    outN1[j] += (double)(Dsum - inputSum) / DsumMax * pSumMaxAlgorithm;
							}
                    }
                }
				a1[j] = sigmoid(outN1[j]);
				da1[j] = round_th(a1[j]*(param->numInputLevel-1), param->Hthreshold);
			}

			numBatchReadSynapse = (int)ceil((double)param->nHide/param->numColMuxed);
			#pragma omp critical    // Use critical here since NeuroSim class functions may update its member variables
			for (int j=0; j<param->nHide; j+=numBatchReadSynapse) {
				int numActiveRows = 0;  // Number of selected rows for NeuroSim
				for (int n=0; n<param->numBitInput; n++) {
					for (int k=0; k<param->nInput; k++) {
						if ((dTestInput[i][k]>>n) & 1) {    // if the nth bit of dTestInput[i][k] is 1
							numActiveRows++;
						}
					}
				}
				subArrayIH->activityRowRead = (double)numActiveRows/param->nInput/param->numBitInput;
				sumNeuroSimReadEnergyIH += NeuroSimSubArrayReadEnergy(subArrayIH);
				sumNeuroSimReadEnergyIH += NeuroSimNeuronReadEnergy(subArrayIH, adderIH, muxIH, muxDecoderIH, dffIH, subtractorIH);
				sumReadLatencyIH += NeuroSimSubArrayReadLatency(subArrayIH);
				sumReadLatencyIH += NeuroSimNeuronReadLatency(subArrayIH, adderIH, muxIH, muxDecoderIH, dffIH, subtractorIH);
			}
		} else {    // Algorithm
			for (int j=0; j<param->nHide; j++){
				for (int k=0; k<param->nInput; k++){
					outN1[j] += testInput[i][k] * weight1[j][k];
				}
				a1[j] = sigmoid(outN1[j]);
			}
		}

		/* Second layer from hidden layer to the output layer */
		tempMax = 0;
		countNum = 0;
		std::fill_n(outN2, param->nOutput, 0);
		std::fill_n(a2, param->nOutput, 0);
		if (param->useHardwareInTestingFF) {  // Hardware
			for (int j=0; j<param->nOutput; j++) {
				if (AnalogNVM *temp = dynamic_cast<AnalogNVM*>(arrayHO->cell[0][0])) {  // Analog eNVM
					if (static_cast<eNVM*>(arrayHO->cell[0][0])->cmosAccess) {  // 1T1R
						sumArrayReadEnergyHO += arrayHO->wireGateCapRow * techHO.vdd * techHO.vdd * param->nHide; // All WLs open
					}
				} else if (DigitalNVM *temp = dynamic_cast<DigitalNVM*>(arrayHO->cell[0][0])) {
					if (static_cast<eNVM*>(arrayHO->cell[0][0])->cmosAccess) {  // 1T1R
						sumArrayReadEnergyHO += arrayHO->wireGateCapRow * techHO.vdd * techHO.vdd;  // Selected WL
					} else {    // Cross-point
						sumArrayReadEnergyHO += arrayHO->wireCapRow * techHO.vdd * techHO.vdd * (param->nHide - 1); // Unselected WLs
					}
				}else if (AnalogNVM *temp = dynamic_cast<AnalogNVM*>(arrayHO->cell[0][0]))  // Analog eNVM
						sumArrayReadEnergyHO += arrayHO->wireGateCapRow * techHO.vdd * techHO.vdd * param->nHide; // All WLs open

				for (int n=0; n<param->numBitInput; n++) {
					double pSumMaxAlgorithm = pow(2, n) / (param->numInputLevel - 1) * arrayHO->arrayRowSize;    // Max algorithm partial weighted sum for the nth vector bit (if both max input value and max weight are 1)
					if (AnalogNVM *temp = dynamic_cast<AnalogNVM*>(arrayHO->cell[0][0])) {  // Analog NVM
						double Isum = 0;    // weighted sum current
						double IsumMax = 0; // Max weighted sum current
            double IsumMin = 0;
						double a1Sum = 0;   // Weighted sum current of a1 vector * weight=1 column
						for (int k=0; k<param->nHide; k++) {
							if ((da1[k]>>n) & 1) {    // if the nth bit of da1[k] is 1
								Isum += arrayHO->ReadCell(j,k);
								a1Sum += arrayHO->GetMediumCellReadCurrent(j,k);
								sumArrayReadEnergyHO += arrayHO->wireCapRow * readVoltageHO * readVoltageHO;  
							}
							IsumMax += arrayHO->GetMaxCellReadCurrent(j,k);
              IsumMin += arrayHO->GetMinCellReadCurrent(j,k);
						}
						sumArrayReadEnergyHO += Isum * readVoltageHO * readPulseWidthHO;
						int outputDigits = 2*(CurrentToDigits(Isum, IsumMax-IsumMin)-CurrentToDigits(a1Sum, IsumMax-IsumMin));
						outN2[j] += DigitsToAlgorithm(outputDigits, pSumMaxAlgorithm);
                        
					} 
                    else 
                        {// SRAM or digital eNVM
                            bool digitalNVM = false; 
                            bool parallelRead = false;
                            if(DigitalNVM*temp = dynamic_cast<DigitalNVM*>(arrayHO->cell[0][0]))
                            {    digitalNVM = true;
                                if(static_cast<DigitalNVM*>(arrayHO->cell[0][0])->parallelRead == true) 
								{
                                    parallelRead = true;
                                }
                            }
                            if(digitalNVM && parallelRead)
                            {
                                double Imin = static_cast<DigitalNVM*>(arrayHO->cell[0][0])->avgMinConductance*static_cast<DigitalNVM*>(arrayHO->cell[0][0])->readVoltage;
                                double Imax = static_cast<DigitalNVM*>(arrayHO->cell[0][0])->avgMaxConductance*static_cast<DigitalNVM*>(arrayHO->cell[0][0])->readVoltage;
                                double Isum = 0;    // weighted sum current
                                double IsumMax = 0; // Max weighted sum current
                                double inputSum = 0;    // Weighted sum current of input vector * weight=1 column
                                int Dsum=0;
                                int DsumMax = 0;
                                int Dref = 0;
                                for (int w=0;w<param->numWeightBit;w++){
                                    int colIndex = (j+1) * param->numWeightBit - (w+1);  // w=0 is the LSB
                                    for (int k=0; k<param->nHide; k++) {
                                        if ((da1[k]>>n) & 1) { // accumulate the current along a column
                                            Isum += static_cast<DigitalNVM*>(arrayHO->cell[colIndex][k])->conductance*static_cast<DigitalNVM*>(arrayHO->cell[colIndex][k])->readVoltage;
                                            //inputSum += Imin;
                                            inputSum += static_cast<DigitalNVM*>(arrayHO->cell[arrayHO->refColumnNumber][k])->conductance*static_cast<DigitalNVM*>(arrayHO->cell[arrayHO->refColumnNumber][k])->readVoltage;                                            
                                        }
                                    }
                                    int outputDigits = (int) (Isum /(Imax-Imin)); 
                                    int outputDigitsRef = (int) (inputSum/(Imax-Imin));
                                    outputDigits = outputDigits-outputDigitsRef;
 
                                    Dref = (int)(inputSum/Imin);
                                    Isum=0;
                                    inputSum=0;
                                    Dsum += outputDigits*(int) pow(2,w);  // get the weight represented by the column
                                    DsumMax += param->nHide*(int) pow(2,w); // the maximum weight that can be represented by this column                                        
                                }
                                sumArrayReadEnergyHO += static_cast<DigitalNVM*>(arrayHO->cell[0][0])->readEnergy * arrayHO->numCellPerSynapse * arrayHO->arrayRowSize;
                                outN2[j] += (double)(Dsum - Dref*(pow(2,param->numWeightBit-1)-1)) / DsumMax * pSumMaxAlgorithm;
                            }
                            else
                            {                            
							    int Dsum = 0;
							    int DsumMax = 0;
							    int a1Sum = 0;
							    for (int k=0; k<param->nHide; k++) {
								    if ((da1[k]>>n) & 1) {    // if the nth bit of da1[k] is 1
									    Dsum += (int)(arrayHO->ReadCell(j,k));
									    a1Sum += pow(2, arrayHO->numCellPerSynapse-1) - 1;    // get current of Dummy Column as reference
								    }
								    DsumMax += pow(2, arrayHO->numCellPerSynapse) - 1;
							    } 
							    if (DigitalNVM *temp = dynamic_cast<DigitalNVM*>(arrayHO->cell[0][0])) {    // Digital eNVM
								    sumArrayReadEnergyHO += static_cast<DigitalNVM*>(arrayHO->cell[0][0])->readEnergy * arrayHO->numCellPerSynapse * arrayHO->arrayRowSize;
							    } 
                                else {
								    sumArrayReadEnergyHO += static_cast<SRAM*>(arrayHO->cell[0][0])->readEnergy * arrayHO->numCellPerSynapse * arrayHO->arrayRowSize;
							    }
							    outN2[j] += (double)(Dsum - a1Sum) / DsumMax * pSumMaxAlgorithm;
                            }
						} 
				}
				a2[j] = sigmoid(outN2[j]);
				if (a2[j] > tempMax) {
					tempMax = a2[j];
					countNum = j;
				}
			}

			numBatchReadSynapse = (int)ceil((double)param->nOutput/param->numColMuxed);
			#pragma omp critical    // Use critical here since NeuroSim class functions may update its member variables
			for (int j=0; j<param->nOutput; j+=numBatchReadSynapse) {
				int numActiveRows = 0;  // Number of selected rows for NeuroSim
				for (int n=0; n<param->numBitInput; n++) {
					for (int k=0; k<param->nHide; k++) {
						if ((da1[k]>>n) & 1) {    // if the nth bit of da1[k] is 1
							numActiveRows++;
						}
					}
				}
				subArrayHO->activityRowRead = (double)numActiveRows/param->nHide/param->numBitInput;
				sumNeuroSimReadEnergyHO += NeuroSimSubArrayReadEnergy(subArrayHO);
				sumNeuroSimReadEnergyHO += NeuroSimNeuronReadEnergy(subArrayHO, adderHO, muxHO, muxDecoderHO, dffHO, subtractorHO);
				sumReadLatencyHO += NeuroSimSubArrayReadLatency(subArrayHO);
				sumReadLatencyHO += NeuroSimNeuronReadLatency(subArrayHO, adderHO, muxHO, muxDecoderHO, dffHO, subtractorHO);
			}
		} else {    // Algorithm
			for (int j=0; j<param->nOutput; j++) {
				for (int k=0; k<param->nHide; k++) {
					outN2[j] += a1[k] * weight2[j][k];
				}
				a2[j] = sigmoid(outN2[j]);
				if (a2[j] > tempMax) {
					tempMax = a2[j];
					countNum = j;
				}
			}
		}
		if (testOutput[i][countNum] == 1) {
			correct++;
		}
	}
	if (!param->useHardwareInTraining) {    // Calculate the classification latency and energy only for offline classification
		arrayIH->readEnergy += sumArrayReadEnergyIH;
		subArrayIH->readDynamicEnergy += sumNeuroSimReadEnergyIH;
		arrayHO->readEnergy += sumArrayReadEnergyHO;
		subArrayHO->readDynamicEnergy += sumNeuroSimReadEnergyHO;
		subArrayIH->readLatency += sumReadLatencyIH;
		subArrayHO->readLatency += sumReadLatencyHO;
	}
}

