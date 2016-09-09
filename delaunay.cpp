
#include "stdio.h"

#include "delaunay.h"
#include "../../Corelib/CommonFuncs.h"


void GenerateDelaunayTriangle(CTINClass* pTin, double* px, double* py, int np)
{
	int i,j;
	
	//CTINClass* pTin = NULL;		
	//pTin = new CTINClass("aaa");		

	pTin->BeginAddPoints();
	for(i=0; i<np; i++)
	{
		pTin->AddPoint(px[i], py[i], i);
	}
	pTin->EndAddPoints();
	pTin->FastConstruct();
	pTin->EnableIntersection( false );
}



void DrawDelaunayTriangle(IplImage* pImage, CTINClass* pTin)
{	
	//draw triangle
	LONG nTriangleNum = 0 ;		
	TRIANGLE **tris = pTin->SelectTriangles(&nTriangleNum, 0, 0, 0, 0 );	
	CvPoint pts[3];
	int i,k;
	int index;

	for(i=0; i<nTriangleNum; i++)
	{		
		for( k=0; k<3; k++)
		{		
			pts[k].x = (*tris)->vertex[k]->x;
			pts[k].y = (*tris)->vertex[k]->y;
			index    = (*tris)->vertex[k]->attr;
		}
		cvLine(pImage, pts[0], pts[1], CV_RGB(255,255,255), 1, 8, 0);
		cvLine(pImage, pts[1], pts[2], CV_RGB(255,255,255), 1, 8, 0);
		cvLine(pImage, pts[0], pts[2], CV_RGB(255,255,255), 1, 8, 0);
		*tris++;	
	}	
}

//given the point, find the triangle containing this point
int  SelectTriangle(double x, double y, CTINClass* pTin, int* pIndex)
{
	LONG nTriangleNum = 0 ;		
	TRIANGLE **tris = pTin->SelectTriangles(&nTriangleNum, 0, 0, 0, 0 );	
	double px[3];
	double py[3];
	int i,k;
	int index[3];
	bool bIsFind = false;

	for(i=0; i<nTriangleNum; i++)
	{		
		for( k=0; k<3; k++)
		{		
			px[k] = (*tris)->vertex[k]->x;
			py[k] = (*tris)->vertex[k]->y;
			index[k]  = (*tris)->vertex[k]->attr;
		}
		
		//judge if the point is in the triangle
		if( pnpoly(3, px, py, x, y) )
		{
			for(k=0; k<3; k++)
				pIndex[k] = index[k];
			
			return 1;
		}	

		*tris++;	
	}	

	return 0;
}



/* the script demostrates iterative construction of
delaunay triangulation and voronoi tesselation */

CvSubdiv2D* init_delaunay( CvMemStorage* storage, CvRect rect )
{
	CvSubdiv2D* subdiv;

	subdiv = cvCreateSubdiv2D( CV_SEQ_KIND_SUBDIV2D, sizeof(*subdiv),
		sizeof(CvSubdiv2DPoint),
		sizeof(CvQuadEdge2D),
		storage );
	cvInitSubdivDelaunay2D( subdiv, rect );

	return subdiv;
}


void draw_subdiv_point( IplImage* img, CvPoint2D32f fp, CvScalar color )
{
	cvCircle( img, cvPoint(cvRound(fp.x), cvRound(fp.y)), 3, color, CV_FILLED, 8, 0 );
}


void draw_subdiv_edge( IplImage* img, CvSubdiv2DEdge edge, CvScalar color )
{
	CvSubdiv2DPoint* org_pt;
	CvSubdiv2DPoint* dst_pt;
	CvPoint2D32f org;
	CvPoint2D32f dst;
	CvPoint iorg, idst;

	org_pt = cvSubdiv2DEdgeOrg(edge);
	dst_pt = cvSubdiv2DEdgeDst(edge);

	if( org_pt && dst_pt )
	{
		org = org_pt->pt;
		dst = dst_pt->pt;

		iorg = cvPoint( cvRound( org.x ), cvRound( org.y ));
		idst = cvPoint( cvRound( dst.x ), cvRound( dst.y ));

		cvLine( img, iorg, idst, color, 1, CV_AA, 0 );
	}
}


void draw_subdiv( IplImage* img, CvSubdiv2D* subdiv,
				 CvScalar delaunay_color, CvScalar voronoi_color )
{
	CvSeqReader  reader;
	int i, total = subdiv->edges->total;
	int elem_size = subdiv->edges->elem_size;

	cvStartReadSeq( (CvSeq*)(subdiv->edges), &reader, 0 );

	for( i = 0; i < total; i++ )
	{
		CvQuadEdge2D* edge = (CvQuadEdge2D*)(reader.ptr);

		if( CV_IS_SET_ELEM( edge ))
		{
			draw_subdiv_edge( img, (CvSubdiv2DEdge)edge + 1, voronoi_color );
			draw_subdiv_edge( img, (CvSubdiv2DEdge)edge, delaunay_color );
		}

		CV_NEXT_SEQ_ELEM( elem_size, reader );
	}
}


void locate_point( CvSubdiv2D* subdiv, CvPoint2D32f fp, IplImage* img,
				  CvScalar active_color )
{
	CvSubdiv2DEdge e;
	CvSubdiv2DEdge e0 = 0;
	CvSubdiv2DPoint* p = 0;

	cvSubdiv2DLocate( subdiv, fp, &e0, &p );

	if( e0 )
	{
		e = e0;
		do
		{
			draw_subdiv_edge( img, e, active_color );
			e = cvSubdiv2DGetEdge(e,CV_NEXT_AROUND_LEFT);
		}
		while( e != e0 );
	}

	draw_subdiv_point( img, fp, active_color );
}


void draw_subdiv_facet( IplImage* img, CvSubdiv2DEdge edge )
{
	CvSubdiv2DEdge t = edge;
	int i, count = 0;
	CvPoint* buf = 0;

	// count number of edges in facet
	do
	{
		count++;
		t = cvSubdiv2DGetEdge( t, CV_NEXT_AROUND_LEFT );
	} while (t != edge );

	buf = (CvPoint*)malloc( count * sizeof(buf[0]));

	// gather points
	t = edge;
	for( i = 0; i < count; i++ )
	{
		CvSubdiv2DPoint* pt = cvSubdiv2DEdgeOrg( t );
		if( !pt ) break;
		buf[i] = cvPoint( cvRound(pt->pt.x), cvRound(pt->pt.y));
		t = cvSubdiv2DGetEdge( t, CV_NEXT_AROUND_LEFT );
	}

	if( i == count )
	{
		CvSubdiv2DPoint* pt = cvSubdiv2DEdgeDst( cvSubdiv2DRotateEdge( edge, 1 ));
		for(int k=0; k<count; k++)
		{
			printf("%d %d  ", buf[k].x, buf[k].y );
		}
		printf("\n");

		cvFillConvexPoly( img, buf, count, CV_RGB(rand()&255,rand()&255,rand()&255), CV_AA, 0 );
		cvPolyLine( img, &buf, &count, 1, 1, CV_RGB(0,0,0), 1, CV_AA, 0);
		draw_subdiv_point( img, pt->pt, CV_RGB(0,0,0));
         
		printf("%lf %lf \n", pt->pt.x, pt->pt.y);
		printf("count: %d \n\n", count);

	}
	free( buf );
}

void paint_voronoi( CvSubdiv2D* subdiv, IplImage* img )
{
	CvSeqReader  reader;
	int i, total = subdiv->edges->total;
	int elem_size = subdiv->edges->elem_size;

	cvCalcSubdivVoronoi2D( subdiv );

	cvStartReadSeq( (CvSeq*)(subdiv->edges), &reader, 0 );

	printf("total: %d \n", total);
	for( i=0; i<total; i++ )
	{
		CvQuadEdge2D* edge = (CvQuadEdge2D*)(reader.ptr);
		if( CV_IS_SET_ELEM( edge ))
		{
			CvSubdiv2DEdge e = (CvSubdiv2DEdge)edge;

			// left
			draw_subdiv_facet( img, cvSubdiv2DRotateEdge( e, 1 ));

			// right
			draw_subdiv_facet( img, cvSubdiv2DRotateEdge( e, 3 ));
		}
		CV_NEXT_SEQ_ELEM( elem_size, reader );
	}
}

