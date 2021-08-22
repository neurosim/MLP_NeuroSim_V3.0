#ifndef WLNewDecoderDriver_H_
#define WLNewDecoderDriver_H_

#include "FunctionUnit.h"
#include "constant.h"
#include "typedef.h"
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"

class WLNewDecoderDriver : public FunctionUnit {
public:
	WLNewDecoderDriver(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell);
	virtual ~WLNewDecoderDriver();
	const InputParameter& inputParameter;
	const Technology& tech;
	const MemCell& cell;

	/* Functions */
	void PrintProperty(const char* str);
	void Initialize(int _numWLRow);
	void CalculateArea(double _newHeight, double _newWidth, AreaModify _option);
	void CalculateLatency(double _rampInput, double _capLoad, double _resLoad, double numRead, double numWrite);
	void CalculatePower(double numRead, double numWrite);
	

	/* Properties */
	bool initialized;	/* Initialization flag */
	bool invalid;      /*Invalidatio flag */
	//double inputCap;	/* Input capacitance, unit: F */
	double capLoad;	/* Output capacitance, unit: F */
	double resLoad;	/* Output resistance, unit: ohm */
	int numWLRow;
	double widthNandN, widthNandP, widthInvN, widthInvP, widthTgN, widthTgP;
	double capNandInput, capNandOutput, capInvInput, capInvOutput, capTgGateN, capTgGateP, capTgDrain;
	double resTg;
	double rampInput, rampOutput;
	bool multifunctional;
	bool neuro;
};

#endif /* WLNewDecoderDriver_H_ */