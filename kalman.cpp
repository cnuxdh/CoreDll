

#include "cv.h"
#include "highgui.h"
#include "cvtypes.h"

#include "kalman.h"


CKalManFilter::CKalManFilter()
{
	const float A[] = { 1, 1, 0, 1 };

	//IplImage* img = cvCreateImage( cvSize(500,500), 8, 3 );
	kalman = cvCreateKalman( 2, 1, 0 );
	state  = cvCreateMat( 2, 1, CV_32FC1 );      
	process_noise = cvCreateMat( 2, 1, CV_32FC1 );
	measurement = cvCreateMat( 1, 1, CV_32FC1 );

	memcpy( kalman->transition_matrix->data.fl, A, sizeof(A));
	cvSetIdentity( kalman->measurement_matrix, cvRealScalar(1) );
	cvSetIdentity( kalman->process_noise_cov, cvRealScalar(1e-5) );
	cvSetIdentity( kalman->measurement_noise_cov, cvRealScalar(1e-1) );
	cvSetIdentity( kalman->error_cov_post, cvRealScalar(1));
}


CKalManFilter::~CKalManFilter()
{
	cvReleaseKalman(&kalman);
	cvReleaseMat(&state);
	cvReleaseMat(&process_noise);
	cvReleaseMat(&measurement);
}

double CKalManFilter::Predict(double in)
{
	CvRNG rng = cvRNG(-1);

	//predict
	const CvMat* prediction = cvKalmanPredict( kalman, 0 );
	float predict_value = prediction->data.fl[0];

	//update measurement
	measurement->data.fl[0] = in;
	//cvMatMulAdd( kalman->measurement_matrix, state, measurement, measurement );
	//float measurement_value = measurement->data.fl[0];

	//correct
	cvKalmanCorrect( kalman, measurement );
	//cvRandArr( &rng, process_noise, CV_RAND_NORMAL, cvRealScalar(0),
	//	         cvRealScalar(sqrt(kalman->process_noise_cov->data.fl[0])));
	//cvMatMulAdd( kalman->transition_matrix, state, process_noise, state );

	//out = predict_value;
	return predict_value;
}