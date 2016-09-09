
#ifndef KALMAN_H
#define KALMAN_H

#include "exports.h"

class DLL_EXPORT CKalManFilter
{
public:
	CKalManFilter();
	//CKalManFilter(int nDim);
	~CKalManFilter();

	double Predict(double in);

private:
	CvKalman* kalman;
	CvMat* state;
	CvMat* process_noise;
	CvMat* measurement;	
};








#endif
