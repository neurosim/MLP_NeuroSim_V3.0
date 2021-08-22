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

#ifndef CONSTANT_H_
#define CONSTANT_H_

#define INV		0
#define NOR		1
#define NAND	2

#define NMOS	0
#define PMOS	1

#define MAX_NMOS_SIZE	100
#define MIN_NMOS_SIZE	2	//1.5

#define MAX_TRANSISTOR_HEIGHT 28

#define MIN_GAP_BET_P_AND_N_DIFFS	3.5 //2
#define MIN_GAP_BET_SAME_TYPE_DIFFS	1.6	//1.5
#define MIN_GAP_BET_GATE_POLY		2.8	//1.5
#define MIN_GAP_BET_CONTACT_POLY	0.7	//0.75
#define CONTACT_SIZE				1.3	//1
#define MIN_WIDTH_POWER_RAIL		3.4	//2
#define MIN_POLY_EXT_DIFF			1.0	// Minimum poly extension beyond diffusion region
#define MIN_GAP_BET_FIELD_POLY		1.6	// Field poly means the poly above the field oxide (outside the active region)
#define POLY_WIDTH					1.0
#define M2_PITCH					3.2
#define M3_PITCH					2.8

#define AVG_RATIO_LEAK_2INPUT_NAND 0.48
#define AVG_RATIO_LEAK_3INPUT_NAND 0.31
#define AVG_RATIO_LEAK_2INPUT_NOR  0.95
#define AVG_RATIO_LEAK_3INPUT_NOR  0.62

#define W_SENSE_P		7.5
#define W_SENSE_N		3.75
#define W_SENSE_ISO		12.5
#define W_SENSE_EN		5.0
#define W_SENSE_MUX		9.0

#define IR_DROP_TOLERANCE 			0.1

#define HEIGHT_WIDTH_RATIO_LIMIT	5

#define RATIO_READ_THRESHOLD_VS_VOLTAGE	0.2

#define	ROW_MODE	0	// Connect to rows
#define	COL_MODE	1	// Connect to columns

#endif /* CONSTANT_H_ */

