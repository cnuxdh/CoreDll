#include "ransac.h"





#define T_DIST 4                // thres. for distance in RANSAC algorithm



//****************************************************************
// Compute the homography matrix H
// i.e., solve the optimization problem min ||Ah||=0 s.t. ||h||=1
// where A is 2n*9, h is 9*1
// input: n (number of pts pairs)
// p1, p2 (coresponded pts pairs x and x')
// output: 3*3 matrix H
//****************************************************************
void ComputeH(int n, CvPoint2D64f *p1, CvPoint2D64f *p2, CvMat *H)
{
	int i;
	CvMat *A = cvCreateMat(2*n, 9, CV_64FC1);
	CvMat *U = cvCreateMat(2*n, 2*n, CV_64FC1);
	CvMat *D = cvCreateMat(2*n, 9, CV_64FC1);
	CvMat *V = cvCreateMat(9, 9, CV_64FC1);
	
	cvZero(A);
	for(i=0; i<n; i++)
	{
		// 2*i row
		
		cvmSet(A,2*i,3,-p1[i].x);
		cvmSet(A,2*i,4,-p1[i].y);
		cvmSet(A,2*i,5,-1);
		cvmSet(A,2*i,6,p2[i].y*p1[i].x);
		cvmSet(A,2*i,7,p2[i].y*p1[i].y);
		cvmSet(A,2*i,8,p2[i].y);
		
		// 2*i+1 row
		cvmSet(A,2*i+1,0,p1[i].x);
		cvmSet(A,2*i+1,1,p1[i].y);
		cvmSet(A,2*i+1,2,1);
		cvmSet(A,2*i+1,6,-p2[i].x*p1[i].x);
		cvmSet(A,2*i+1,7,-p2[i].x*p1[i].y);
		cvmSet(A,2*i+1,8,-p2[i].x);
	}
	
	// SVD
	// The flags cause U and V to be returned transposed
	// Therefore, in OpenCV, A = U^T D V
	cvSVD(A, D, U, V, CV_SVD_U_T|CV_SVD_V_T);
	
	// take the last column of V^T, i.e., last row of V
	for(i=0; i<9; i++)
		cvmSet(H, i/3, i%3, cvmGet(V, 8, i));
	cvReleaseMat(&A);
	cvReleaseMat(&U);
	cvReleaseMat(&D);
	cvReleaseMat(&V);
}


//**********************************************************************
// Compute number of inliers by computing distance under a perticular H
// distance = d(Hx, x') + d(invH x', x)
// input: num (number of pts pairs)
// p1, p2 (coresponded pts pairs x and x')
// H (the homography matrix)
// output: inlier_mask (masks to indicate pts of inliers in p1, p2)
// dist_std (std of the distance among all the inliers)
// return: number of inliers
//**********************************************************************
int ComputeNumberOfInliers(int num, CvPoint2D64f *p1, CvPoint2D64f *p2, CvMat *H,
						   CvMat *inlier_mask, double *dist_std)
{
	int i, num_inlier;
	double curr_dist, sum_dist, mean_dist;
	CvPoint2D64f tmp_pt;
	CvMat *dist = cvCreateMat(num, 1, CV_64FC1);
	CvMat *x = cvCreateMat(3,1,CV_64FC1);
	CvMat *xp = cvCreateMat(3,1,CV_64FC1);
	CvMat *pt = cvCreateMat(3,1,CV_64FC1);
	CvMat *invH = cvCreateMat(3,3,CV_64FC1);
	
	cvInvert(H, invH);
	// check each correspondence
	sum_dist = 0;
	num_inlier = 0;
	cvZero(inlier_mask);
	for(i=0; i<num; i++)
	{
		// initial point x
		cvmSet(x,0,0,p1[i].x);
		cvmSet(x,1,0,p1[i].y);
		cvmSet(x,2,0,1);
		
		// initial point x'
		cvmSet(xp,0,0,p2[i].x);
		cvmSet(xp,1,0,p2[i].y);
		cvmSet(xp,2,0,1);
		
		// d(Hx, x')
		cvMatMul(H, x, pt);
		tmp_pt.x = (int)(cvmGet(pt,0,0)/cvmGet(pt,2,0));
		tmp_pt.y = (int)(cvmGet(pt,1,0)/cvmGet(pt,2,0));
		curr_dist = pow(tmp_pt.x-p2[i].x, 2.0) + pow(tmp_pt.y-p2[i].y, 2.0);
		
		// d(x, invH x')
		cvMatMul(invH, xp, pt);
		tmp_pt.x = (int)(cvmGet(pt,0,0)/cvmGet(pt,2,0));
		tmp_pt.y = (int)(cvmGet(pt,1,0)/cvmGet(pt,2,0));
		curr_dist += pow(tmp_pt.x-p1[i].x, 2.0) + pow(tmp_pt.y-p1[i].y, 2.0);
		
		if(curr_dist < T_DIST)
		{
			// an inlier
			num_inlier++;
			cvmSet(inlier_mask,i,0,1);
			cvmSet(dist,i,0,curr_dist);
			sum_dist += curr_dist;
		}
	}
	
	// Compute the standard deviation of the distance
	mean_dist = sum_dist/(double)num_inlier;
	*dist_std = 0;
	for(i=0; i<num; i++){
		if(cvmGet(inlier_mask,i,0) == 1)
			*dist_std += pow(cvmGet(dist,i,0)-mean_dist,2.0);
	}
	*dist_std /= (double) (num_inlier -1);
	cvReleaseMat(&dist);
	cvReleaseMat(&x);
	cvReleaseMat(&xp);
	cvReleaseMat(&pt);
	cvReleaseMat(&invH);
	return num_inlier;
}


//*****************************************
// Check colinearity of a set of pts
// input: p (pts to be checked)
// num (ttl number of pts)
// return true if some pts are coliner
// false if not
//*****************************************
bool isColinear(int num, CvPoint2D64f *p)
{
	int i,j,k;
	bool iscolinear;
	double value;
	CvMat *pt1 = cvCreateMat(3,1,CV_64FC1);
	CvMat *pt2 = cvCreateMat(3,1,CV_64FC1);
	CvMat *pt3 = cvCreateMat(3,1,CV_64FC1);
	CvMat *line = cvCreateMat(3,1,CV_64FC1);
	
	iscolinear = false;
	// check for each 3 points combination
	for(i=0; i<num-2; i++){
		cvmSet(pt1,0,0,p[i].x);
		cvmSet(pt1,1,0,p[i].y);
		cvmSet(pt1,2,0,1);
		for(j=i+1; j<num-1; j++){
			cvmSet(pt2,0,0,p[j].x);
			cvmSet(pt2,1,0,p[j].y);
			cvmSet(pt2,2,0,1);
			// compute the line connecting pt1 & pt2
			cvCrossProduct(pt1, pt2, line);
			for(k=j+1; k<num; k++){
				cvmSet(pt3,0,0,p[k].x);
				cvmSet(pt3,1,0,p[k].y);
				cvmSet(pt3,2,0,1);
				// check whether pt3 on the line
				value = cvDotProduct(pt3, line);
				if(abs(value) < 10e-2){
					iscolinear = true;
					break;
				}
			}
			if(iscolinear == true) break;
		}
		if(iscolinear == true) break;
	}
	cvReleaseMat(&pt1);
	cvReleaseMat(&pt2);
	cvReleaseMat(&pt3);
	cvReleaseMat(&line);
	return iscolinear;
}

//**********************************************************************
// finding the normalization matrix x' = T*x, where T={s,0,tx, 0,s,ty, 0,0,1}
// compute T such that the centroid of x' is the coordinate origin (0,0)T
// and the average distance of x' to the origin is sqrt(2)
// we can derive that tx = -scale*mean(x), ty = -scale*mean(y),
// scale = sqrt(2)/(sum(sqrt((xi-mean(x)^2)+(yi-mean(y))^2))/n)
// where n is the total number of points
// input: num (ttl number of pts)
// p (pts to be normalized)
// output: T (normalization matrix)
// p (normalized pts)
// NOTE: because of the normalization process, the pts coordinates should
// has accurcy as "float" or "double" instead of "int"
//**********************************************************************
void Normalization(int num, CvPoint2D64f *p, CvMat *T)
{
	double scale, tx, ty;
	double meanx, meany;
	double value;
	int i;
	CvMat *x = cvCreateMat(3,1,CV_64FC1);
	CvMat *xp = cvCreateMat(3,1,CV_64FC1);
	
	meanx = 0;
	meany = 0;
	for(i=0; i<num; i++){
		meanx += p[i].x;
		meany += p[i].y;
	}
	meanx /= (double)num;
	meany /= (double)num;
	
	value = 0;
	for(i=0; i<num; i++)
		value += sqrt(pow(p[i].x-meanx, 2.0) + pow(p[i].y-meany, 2.0));
	value /= (double)num;
	
	scale = sqrt(2.0)/value;
	tx = -scale * meanx;
	ty = -scale * meany;
	
	cvZero(T);
	cvmSet(T,0,0,scale);
	cvmSet(T,0,2,tx);
	cvmSet(T,1,1,scale);
	cvmSet(T,1,2,ty);
	cvmSet(T,2,2,1.0);
	
	//Transform x' = T*x
	for(i=0; i<num; i++){
		cvmSet(x,0,0,p[i].x);
		cvmSet(x,1,0,p[i].y);
		cvmSet(x,2,0,1.0);
		cvMatMul(T,x,xp);
		p[i].x = cvmGet(xp,0,0)/cvmGet(xp,2,0);
		p[i].y = cvmGet(xp,1,0)/cvmGet(xp,2,0);
	}
	cvReleaseMat(&x);
	cvReleaseMat(&xp);
}


//*****************************************************************************
// RANSAC algorithm
// input: num (ttl number of pts)
// m1, m2 (pts pairs)
// output: inlier_mask (indicate inlier pts pairs in (m1, m2) as 1; outlier: 0)
// H (the best homography matrix) m2=H.m1
//*****************************************************************************
void RANSAC_homography(int num, CvPoint2D64f *m1, CvPoint2D64f *m2, CvMat *H,CvMat *inlier_mask)
{
	int i,j;
	int N = 1000, s = 4, sample_cnt = 0;
	double e, p = 0.99;
	int numinlier, MAX_num;
	double curr_dist_std, dist_std;
	bool iscolinear;
	CvPoint2D64f *curr_m1 = new CvPoint2D64f[s];
	CvPoint2D64f *curr_m2 = new CvPoint2D64f[s];
	int *curr_idx = new int[s];
	
	CvMat *curr_inlier_mask = cvCreateMat(num,1,CV_64FC1);
	CvMat *curr_H = cvCreateMat(3,3,CV_64FC1);
	CvMat *T1 = cvCreateMat(3,3,CV_64FC1);
	CvMat *T2 = cvCreateMat(3,3,CV_64FC1);
	CvMat *invT2 = cvCreateMat(3,3,CV_64FC1);
	CvMat *tmp_pt = cvCreateMat(3,1,CV_64FC1);
	
	// RANSAC algorithm (reject outliers and obtain the best H)
	srand(134);
	MAX_num = -1;
	while(N > sample_cnt)
	{
		// for a randomly chosen non-colinear correspondances
		iscolinear = true;
		while(iscolinear == true)
		{
			iscolinear = false;
			for(i=0; i<s; i++)
			{
				// randomly select an index
				curr_idx[i] = rand()%num;
				for(j=0; j<i; j++)
				{
					if(curr_idx[i] == curr_idx[j])
					{
						iscolinear = true;
						break;
					}
				}
				if(iscolinear == true) break;
				curr_m1[i].x = m1[curr_idx[i]].x;
				curr_m1[i].y = m1[curr_idx[i]].y;
				curr_m2[i].x = m2[curr_idx[i]].x;
				curr_m2[i].y = m2[curr_idx[i]].y;
			}
			
			// Check whether these points are colinear
			if(iscolinear == false)
				iscolinear = isColinear(s, curr_m1);
		}
		
		// Nomalized DLT
		Normalization(s, curr_m1, T1); //curr_m1 <- T1 * curr_m1
		Normalization(s, curr_m2, T2); //curr_m2 <- T2 * curr_m2
		
		// Compute the homography matrix H = invT2 * curr_H * T1
		ComputeH(s, curr_m1, curr_m2, curr_H);
		cvInvert(T2, invT2);
		cvMatMul(invT2, curr_H, curr_H); // curr_H <- invT2 * curr_H
		cvMatMul(curr_H, T1, curr_H); // curr_H <- curr_H * T1
		
		// Calculate the distance for each putative correspondence
		// and compute the number of inliers
		numinlier =ComputeNumberOfInliers(num, m1, m2, curr_H, curr_inlier_mask, &curr_dist_std);
		
		// Update a better H
		if(numinlier > MAX_num || (numinlier == MAX_num && curr_dist_std <dist_std))
		{
			MAX_num = numinlier;
			cvCopy(curr_H, H);
			cvCopy(curr_inlier_mask, inlier_mask);
			dist_std = curr_dist_std;
		}
		// update number N by Algorithm 4.5
		e = 1 - (double)numinlier / (double)num;
		N = (int)(log(1-p)/log(1-pow(1-e,s)));
		sample_cnt++;
	}
	// Optimal estimation using all the inliers
	delete curr_m1, curr_m2, curr_idx;
	cvReleaseMat(&curr_H);
	cvReleaseMat(&T1);
	cvReleaseMat(&T2);
	cvReleaseMat(&invT2);
	cvReleaseMat(&tmp_pt);
	cvReleaseMat(&curr_inlier_mask);
}
