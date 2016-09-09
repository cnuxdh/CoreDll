
#ifndef DELAUNAY_H
#define DELAUNAY_H

#include "exports.h"

#include "../../Corelib/Triangles.h"


#include "cv.h"
#include "highgui.h"


void DLL_EXPORT GenerateDelaunayTriangle(CTINClass* pTin, double* px, double* py, int np);
void DLL_EXPORT DrawDelaunayTriangle(IplImage* pImage, CTINClass* pTin);
int  DLL_EXPORT SelectTriangle(double x, double y, CTINClass* pTin, int* pIndex);


DLL_EXPORT CvSubdiv2D* init_delaunay( CvMemStorage* storage, CvRect rect );
DLL_EXPORT void draw_subdiv_point( IplImage* img, CvPoint2D32f fp, CvScalar color );
DLL_EXPORT void draw_subdiv_edge( IplImage* img, CvSubdiv2DEdge edge, CvScalar color );
DLL_EXPORT void draw_subdiv( IplImage* img, CvSubdiv2D* subdiv,
				 CvScalar delaunay_color, CvScalar voronoi_color );

DLL_EXPORT void locate_point( CvSubdiv2D* subdiv, CvPoint2D32f fp, IplImage* img,
				              CvScalar active_color );

DLL_EXPORT void draw_subdiv_facet( IplImage* img, CvSubdiv2DEdge edge );

DLL_EXPORT void paint_voronoi( CvSubdiv2D* subdiv, IplImage* img );





#endif