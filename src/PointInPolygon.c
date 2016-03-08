#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <complex.h>
#include <limits.h>
#define min(a, b) (((a) < (b)) ? (a) : (b));

struct pointInSpace{
	double x;
	double y;
	double z;
};
struct cube3dcoord{
  int x;
  int y;
  int z;
};
struct planeEquation{
	double a;
	double b;
	double c;
	double d;
};
struct DICOM_OrientationMatrix{
	double a11;		// 1st element of ImageOrientationPatient
	double a21;		// 2nd element of ImageOrientationPatient
	double a31;		// 3rd element of ImageOrientationPatient	
	double a12;		// 4th element of ImageOrientationPatient
	double a22;		// 5th element of ImageOrientationPatient
	double a32;		// 6th element of ImageOrientationPatient
	double Sx;		// 1st element of ImagePositionPatient
	double Sy;		// 2nd element of ImagePositionPatient
	double Sz;		// 3nd element of ImagePositionPatient
	double XpixelSpacing;	// 1st element of Pixel Spacing
	double YpixelSpacing;	// 2nd element of Pixel Spacing
	struct DICOM_OrientationMatrix *next;
	struct DICOM_OrientationMatrix *prev;
};

/*
 * Function for calculating the <x,y,z> coords of a voxel into the space
 * starting from rows (Ny) and columns (Nx) of the 2D Matrix of doses
 * Nx	:		columns
 * Ny 	:		rows
 * DOM  :		Dicom Orientation Matrix
 */ 
struct pointInSpace get3DPosFromNxNy(int Nx, int Ny, struct DICOM_OrientationMatrix DOM) {
       struct pointInSpace pis;
        //printf("\n a11=%lf, Nx=%d Ny=%d, a12=%lf, Sx=%lf  \n",DOM.a11,Nx, Ny,  DOM.a12, DOM.Sx);
        pis.x = DOM.a11 * Nx * DOM.XpixelSpacing + DOM.a12 * Ny * DOM.YpixelSpacing + DOM.Sx;
        pis.y = DOM.a21 * Nx * DOM.XpixelSpacing + DOM.a22 * Ny * DOM.YpixelSpacing + DOM.Sy;
        pis.z = DOM.a31 * Nx * DOM.XpixelSpacing + DOM.a32 * Ny * DOM.YpixelSpacing + DOM.Sz;
        return pis;
}
/*
 * Function for calculating the distance between a point and a plane
 * pS   :		point in space 
 * pE	:		plane equation
 */ 
double getPointPlaneDistance(struct pointInSpace pS, struct planeEquation pE) {
        //printf("1=%lf, 2=%lf  ",(pE.a * pS.x  +  pE.b * pS.y  +  pE.c * pS.z  +  pE.d),sqrt( pE.a * pE.a + pE.b * pE.b + pE.c * pE.c));
        return fabs(pE.a * pS.x  +  pE.b * pS.y  +  pE.c * pS.z  +  pE.d) / sqrt( pE.a * pE.a + pE.b * pE.b + pE.c * pE.c);
}


/*
 * Function for calculating the equation of a plane starting from 3 points in the space
 * Pa   :		1st point
 * Pb	:		2nd point
 * Pc   :		3rd point
 */ 
struct planeEquation getPlaneEquationBetween3Points(struct pointInSpace Pa, struct pointInSpace Pb, struct pointInSpace Pc) {
        struct planeEquation plane;
        plane.a = (Pb.y - Pa.y) * (Pc.z - Pa.z) - (Pc.y - Pa.y) * (Pb.z - Pa.z);
        plane.b = (Pb.z - Pa.z) * (Pc.x - Pa.x) - (Pc.z - Pa.z) * (Pb.x - Pa.x);
        plane.c = (Pb.x - Pa.x) * (Pc.y - Pa.y) - (Pc.x - Pa.x) * (Pb.y - Pa.y);
        plane.d = - ((plane.a * Pa.x) + (plane.b * Pa.y) + (plane.c * Pa.z));
        return plane;
}
/*
 * Function for calculating the distance between two points
 * a{Xa, Ya} and b{Xb, Yb}
 */ 
double distance(double Xa, double Ya, double Xb, double Yb) {
	return (sqrt(pow(Xa - Xb, 2) + pow(Ya - Yb, 2)));	
}

/*
 * Function for calculating the quadratic norm of
 * a complex number {a + bI}
 */

double CNorm(double complex z) {
	return(cabs(z));
}
/*
 * Function for mapping a Polygon in Complex plane
 */
void MapPolygon(double *totalX, double *totalY, int *Offset, double complex *MappedPolygon) {
	int n;

	for (n=1; n< Offset[1]; n++) {  // start counting from the value after the 1st Offset value
		MappedPolygon[n - 1] = totalX[n] + totalY[n] * I;
	}
}

/*
 *  Algorithm for Point In Polygon
 *  calculated over a multiple series of contours each one on a different slice.
 *	Each slice is ordered according NumSlices value.
 *	Variables to be addressed to the C code are:
 *	totalY: 	vector of all Y coordinates of contours (one appended to the other)
 *	totalX: 	vector of all X coordinates of contours (one appended to the other)
 *	offset:		vector of all shifts for moving from one contour to the other
 *	X, Y:		coordinates position of x & y points matrix of single plan of image 
 *	PIPvector:	vector of Points in Polygon to be arranged in an array in R
 *	nX, nY:		number of X & Y coordinates in X & Y
 *	NumSlices:	number of axial slices to be scanned for PIP
 *	FullZ:		vector of slices containing contours {1} or empty {0}
 */
void MultiPIP (double *X, double *Y, double *totalX, double *totalY, int *nX, int *nY,  
		int *NumSlices, int *Offset, int *PIPvector, int *FullZ) {
	int n, i, j, x, y, c;

	for (n=0; n< *NumSlices; n++) {		// loop among slices
		if (FullZ[n]==1) {				// loop for analyzing only full slices
			// loop for filling full slices by PIP values 0: point outside, 1: point inside
			for (y=0; y< *nY; y++) {	// loop through Y axis
				for (x=0; x< *nX; x++) {// loop through X axis
					c=0;				// default value is 0: outside
					for (i = Offset[n], j = Offset[n + 1] - 1; i < Offset[n + 1] ; j = i++) {
						if ( ((totalY[i]>Y[y]) != (totalY[j]>Y[y])) &&
								(X[x] < (totalX[j]-totalX[i]) * (Y[y]-totalY[i]) / (totalY[j]-totalY[i]) + totalX[i]) )
							c = !c;		// each cross changes the value of output by opposite
					}
					PIPvector[n * (*nX) * (*nY) + x + y * (*nX)]= c;  // fill the output values
				}
			}
		}
	}
} 

/*
 *  Algorithm for Point In Polygon on Oblique planes
 *  calculated over a multiple series of contours each one on a different slice.
 *	Each slice is ordered according NumSlices value.
 *	Variables to be addressed to the C code are:
 *	totalY: 	vector of all Y coordinates of contours (one appended to the other)
 *	totalX: 	vector of all X coordinates of contours (one appended to the other)
 *	offset:		vector of all shifts for moving from one contour to the other
 *	X, Y:		coordinates position of x & y points matrix of single plan of image
 *	PIPvector:	vector of Points in Polygon to be arranged in an array in R
 *	nX, nY:		number of X & Y coordinates in X & Y
 *	NumSlices:	number of axial slices to be scanned for PIP
 *	FullZ:		vector of slices containing contours {1} or empty {0}
 *	DICOMv:		DICOM orientation vector, vectorized DICOM orientation matrix multiplied by pixel spacing
 */
void MultiPIPObl (double *totalX, double *totalY, int *nX, int *nY,
		int *NumSlices, int *Offset, int *PIPvector, int *FullZ, double *DICOMv) {
	int n, i, j, x, y, c;
//	double X, Y;
	struct pointInSpace point;
	struct DICOM_OrientationMatrix DOM;


	//XY = (double*)calloc(1, sizeof(double)*(*nrow * *ncol));  // matrix of point to be checked
	for (n=0; n< *NumSlices; n++) {		// loop among slices
		if (FullZ[n]==1) {				// loop for analyzing only full slices
			// loop for filling full slices by PIP values 0: point outside, 1: point inside
			for (y=0; y< *nY; y++) {	// loop through Y axis
				for (x=0; x< *nX; x++) {// loop through X axis
					DOM.a11= DICOMv[n * 9 + 0];	// position 0 in array
					DOM.a21= DICOMv[n * 9 + 1];
					DOM.a31= DICOMv[n * 9 + 2];
					DOM.a12= DICOMv[n * 9 + 3];
					DOM.a22= DICOMv[n * 9 + 4];
					DOM.a32= DICOMv[n * 9 + 5];
					DOM.Sx=  DICOMv[n * 9 + 6];
					DOM.Sy=  DICOMv[n * 9 + 7];
					DOM.Sz=  DICOMv[n * 9 + 8];
					DOM.XpixelSpacing=1;
					DOM.YpixelSpacing=1;
					point=get3DPosFromNxNy(x, y, DOM);
					//printf("\npoint.x %lf point.y %lf, x %d, y %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", point.x, point.y, x, y,DOM.a11,DOM.a21,DOM.a31,DOM.a12,DOM.a22,DOM.a32,DOM.Sx,DOM.Sy,DOM.Sz,DOM.XpixelSpacing,DOM.YpixelSpacing);
					c=0;				// default value is 0: outside
					for (i = Offset[n], j = Offset[n + 1] - 1; i < Offset[n + 1] ; j = i++) {
						if ( ((totalY[i]>point.y) != (totalY[j]>point.y)) &&
								(point.x < (totalX[j]-totalX[i]) * (point.y-totalY[i]) / (totalY[j]-totalY[i]) + totalX[i]) )
							c = !c;		// each cross changes the value of output by opposite
					}
					PIPvector[n * (*nX) * (*nY) + x + y * (*nX)]= c;  // fill the output values
				}
			}
		}
	}
}
/*
 nvert :   number of vertex
 vertx :   x coords of vertex points
 verty :   y coords of vertex points
 
*/

int isThePointInsideThePoly(int nvert, double *vertx, double *verty, double testx, double testy, int fromPosition, int toPosition)
{
    int i, j, c = 0;
    

  for (i = fromPosition, j = fromPosition+nvert-1; i < fromPosition+nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
   (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
  }  
  return c;
}
/*
 *  Algorithm for Point In Polygon on Oblique planes
 *  calculated over a multiple series of contours each one on a different slice.
 *	Each slice is ordered according NumSlices value.
 *	Variables to be addressed to the C code are:
 *      PIPvector       vector of Points in Polygon to be arranged in an array in R
 *	totalY: 	vector of all Y coordinates of contours (one appended to the other)
 *	totalX: 	vector of all X coordinates of contours (one appended to the other)
 *      numberOfPoints: number of points declared in totalX/Y
 *      nX, nY, nZ:	x,y,z dimensions of the voxel matrix
 *      arrayAssociationROIandSlice: associates each point in totalX/Y with a slice. terminator (-100000) means the ROI is closed
 *	arrayDOM:		DICOM orientation vector, vectorized DICOM orientation matrix multiplied by pixel spacing
 * 
 *      minX,maxX,minY,maxY: min, max values of coords in order to preserve time. In future they will be calculated directly into this C code...
 */
void NewMultiPIPObl (int *PIPvector, double *totalX, double *totalY, int *numberOfPoints,
                int *nX, int *nY, int *nZ,
		            int *arrayAssociationROIandSlice, 
                double *arrayDOM,
                double *minX, double *maxX, double *minY, double *maxY) {
	int x, y, z, c, fromPosition, toPosition;
	struct pointInSpace point; 
	struct DICOM_OrientationMatrix DOM;

        int baseCursor = 0;
        for (baseCursor = 0; baseCursor < *numberOfPoints; baseCursor++ ) {
            
           // check if a new ROI is beginning (AND you are not looking the last point!) )
            if ( totalX[ baseCursor ] == -10000 && baseCursor < (*numberOfPoints-1) ) {
                
                // get the corresponding z position
                z = arrayAssociationROIandSlice[ baseCursor + 1 ];
                
                // get the DOM
                DOM.a11= arrayDOM[ z * 9 + 0 ]; DOM.a21= arrayDOM[ z * 9 + 1 ];
                DOM.a31= arrayDOM[ z * 9 + 2 ]; DOM.a12= arrayDOM[ z * 9 + 3 ];
                DOM.a22= arrayDOM[ z * 9 + 4 ]; DOM.a32= arrayDOM[ z * 9 + 5 ];
                DOM.Sx=  arrayDOM[ z * 9 + 6 ]; DOM.Sy=  arrayDOM[ z * 9 + 7 ];
                DOM.Sz=  arrayDOM[ z * 9 + 8 ];
                DOM.XpixelSpacing = 1;  DOM.YpixelSpacing = 1;

                // now Run to find the beginning (clear) and the end of the 
                // points of the ROI, in the given array
                fromPosition = baseCursor;
                for ( toPosition = fromPosition+1; totalX[toPosition]!=-10000; toPosition++) ;;
                toPosition--;                

                // check all the voxel on such z-slice
                // loop through Y axis
                for ( y = 0; y < *nY; y++ ) {	
                        // loop through X axis
                        for ( x = 0; x < *nX; x++ ) {
                                point = get3DPosFromNxNy(x, y, DOM);
                                
                                // check the box, just to avoid redundant computation
                                if(point.x < *minX || point.x > *maxX || point.y < *minY || point.y > *maxY) continue;  

                                // now check if the given point is in the poly
                                c = isThePointInsideThePoly(toPosition-fromPosition-1, totalX, totalY, 
                                        point.x, point.y,  fromPosition+1,  toPosition+1);

                                if( c == 1 )   {                                             
                                    PIPvector[z * (*nX) * (*nY) + x + y * (*nX)] = !PIPvector[z * (*nX) * (*nY) + x + y * (*nX)];
                                }                                
                        }
                }
            }
        }
} 

/*
 * function that calculates the distance between two points "a" and "b"
 * having coordinates (aX, aY) and (bX, bY)
 */
void dist (double *aX, double *aY, double *bX, double *bY, double *dista) {
	//*dist = sqrt(pow(*aX - *bX, 2) + pow(*aY - *bY, 2));
	*dista = distance(*aX, *aY, *bX, *bY);
}

/*
 * Function that computes the quadratic norm
 * exported in R
 */
void QNorm (double *R, double *Img, double *result) {
	*result = CNorm(*R + *Img * I);
} 

/*
 * Function for calculating the minimum distance of a Point matrix from verteces of a Polygon
 * totalX:	vector of X coordinates of Polygon
 * totalY:	vector of Y coordinates of Polygon
 * X:		X coordinates of Points to be tested
 * Y:		Y coordinates of Points to be tested
 * PointV:	vector of result vector of distances
 * num:		length of X vector
 * lP:		length of Polygon
 */
void DistPolygon(double *totalX, double *totalY, double *X, double *Y, double *PointV, int *num, int *lP) {
	int m, n, x, y;					// counters
	double tempD, tempDA, tempDB;	// temp distances, overall, from point A, from point B
	double defD, result;			// other temp values for iteration
	double complex Pz, Az, Bz, z;	// complex numbers for calculation	
	
	m=0;
	for (x=0; x< *num; x++) {
		for (y=0; y< *num; y++) {
			tempD=DBL_MAX;			  // setting the upper bound of distance
			for (n=0; n<*lP-2; n++) { // start mapping all other verteces in the Polygon
				Pz = X[x] + Y[y] * I;
				Az = totalX[n] + totalY[n] * I;
				Bz = totalX[n+1] + totalY[n+1] * I;
				z = (Pz - Az)/(Bz - Az);
				if ((creal(z)<0)||(creal(z)>1)) {  // projection out of the segment
					tempDA=distance(X[x], Y[y], totalX[n], totalY[n]);
					tempDB=distance(X[x], Y[y], totalX[n+1], totalY[n+1]);
					defD=min(tempDA, tempDB);					
				} else {						   // projection inside the segment
					defD= distance(totalX[n+1]*cimag(z), totalY[n+1]*cimag(z), totalX[n]*cimag(z), totalY[n]*cimag(z));				
				}
				result=min(defD, tempD);	
				tempD=result;				
			}			
			PointV[m]=tempD;
			m++;
		}
	}
	//printf("Distance=%lf", tempD);
}

/*
 * test for min function
 */
void Minimum(double *a, double *b, double *c) {
	*c=min(*a, *b);
}

/*
 * test for complex division
 */
void ComplexDivision(double *ra, double *ia, double *rb, double *ib, double *resR, double *resI) {
	double complex za, zb, res;
	za = *ra + *ia * I;
	zb = *rb + *ib * I;
	res = za/zb;
	*resR = creal(res);
	*resI = cimag(res);
	printf("%lf %lfi", *resR, *resI);
}


/*
 * Algorithm for Point In Polygon and Distance from Polygon together
 *  calculated over a multiple series of contours each one on a different slice.
 *	Each slice is ordered according NumSlices value.
 *	Variables to be addressed to the C code are:
 *	totalY: 	vector of all Y coordinates of contours (one appended to the other)
 *	totalX: 	vector of all X coordinates of contours (one appended to the other)
 *	Offset:		vector of all shifts for moving from one contour to the other
 *	X, Y:		coordinates position of x & y points matrix of single plan of image
 *	PIPvector:	vector of Points in Polygon to be arranged in an array in R
 *	nX, nY:		number of X & Y coordinates in X & Y
 *	NumSlices:	number of axial slices to be scanned for PIP
 *	FullZ:		vector of slices containing contours {1} or empty {0}
 *	DISTvector: vector of distances of points from the polygon
 *	Origin:		Origin coordinate of X
 */ 

void MultiPIPDist (double *X, double *Y, double *totalX, double *totalY, int *nX, int *nY,  
  	int *NumSlices, int *Offset, int *PIPvector, int *FullZ, double *DISTvector, double *Origin) {
	int n, i, j, x, y, c;
  double tempD, tempDA, tempDB;	// temp distances, overall, from point A, from point B
	double defD, result;			// other temp values for iteration
	double complex Pz, Az, Bz, z;	// complex numbers for calculation
  
  printf("\nCalculating Point in Polygons...");
	for (n=0; n< *NumSlices; n++) {		// loop among slices
		if (FullZ[n]==1) {				// loop for analyzing only full slices
			// loop for filling full slices by PIP values 0: point outside, 1: point inside
			for (y=0; y< *nY; y++) {	// loop through Y axis
				for (x=0; x< *nX; x++) {// loop through X axis
					c=0;				// default value is 0: outside
					for (i = Offset[n], j = Offset[n + 1] - 1; i < Offset[n + 1] ; j = i++) {
						if ( ((totalY[i]>Y[y]) != (totalY[j]>Y[y])) &&
								(X[x] < (totalX[j]-totalX[i]) * (Y[y]-totalY[i]) / (totalY[j]-totalY[i]) + totalX[i]) )
							c = !c;		// each cross changes the value of output by opposite
					}
					PIPvector[n * (*nX) * (*nY) + x + y * (*nX)]= c;  // fill the output values
				}
			}
		}
	}
  printf("\nCalculating Distances...");
  for (n=0; n< *NumSlices; n++) {
    if (FullZ[n]==1) {
      for (y=0; y< *nY; y++) {
        for (x=0; x< *nX; x++) {
          tempD=DBL_MAX;  // setting the upper bound of distance
          for (i = Offset[n] + 1, j = Offset[n + 1] - 1; i < Offset[n + 1] ; j = i++) {
            if ((totalX[j]!= *Origin) && (totalX[i]!= *Origin)) {  // start the calculation of distances
              Pz = X[x] + Y[y] * I;
              Az = totalX[j] + totalY[j] * I;
              Bz = totalX[i] + totalY[i] * I;
				      z = (Pz - Az)/(Bz - Az);
				      if ((creal(z)<0)||(creal(z)>1)) {  // projection out of the segment
					      tempDA=distance(X[x], Y[y], totalX[j], totalY[j]);
					      tempDB=distance(X[x], Y[y], totalX[i], totalY[i]);
					      defD=min(tempDA, tempDB);					
				        } else {						   // projection inside the segment
					        defD= distance(totalX[i]*cimag(z), totalY[i]*cimag(z), totalX[j]*cimag(z), totalY[j]*cimag(z));				
				        }
				      result=min(defD, tempD);	
				      tempD=result;
            }
            if (PIPvector[n * (*nX) * (*nY) + x + y * (*nX)]==1) {
              DISTvector[n * (*nX) * (*nY) + x + y * (*nX)]=tempD;
              } else {
                DISTvector[n * (*nX) * (*nY) + x + y * (*nX)]=-tempD;
              }
          }
        }
      }
    }
  }
} 


/*
void MultiPIPDist (double *X, double *Y, double *totalX, double *totalY, int *nX, int *nY,
		int *NumSlices, int *Offset, int *PIPvector, int *FullZ, double *DISTvector, double *Origin) {
	int n, m, i, j, x, y, c;
	double tempD, tempDA, tempDB;	// temp distances, overall, from point A, from point B
	double defD, result;			// other temp values for iteration
	double complex Pz, Az, Bz, z;	// complex numbers for calculation

	for (n=0; n< *NumSlices; n++) {		// loop among slices
		if (FullZ[n]==1) {				// loop for analyzing only full slices
			// loop for filling full slices by PIP values 0: point outside, 1: point inside
			for (y=0; y< *nY; y++) {	// loop through Y axis
				for (x=0; x< *nX; x++) {// loop through X axis
					c=0;				// default value is 0: outside
					tempD=DBL_MAX;		// setting the upper bound of distance
					for (i = Offset[n], j = Offset[n + 1] - 1; i < Offset[n + 1] ; j = i++) {
						if ( ((totalY[i]>Y[y]) != (totalY[j]>Y[y])) &&
								(X[x] < (totalX[j]-totalX[i]) * (Y[y]-totalY[i]) / (totalY[j]-totalY[i]) + totalX[i]) )
							c = !c;		// each cross changes the value of output by opposite
						if ((totalX[i]!=*Origin) && (totalX[j]!=*Origin)) {   // start the check of point distance from Polygon
							Pz = X[x] + Y[y] * I;  // mapping points in Imaginary space
							Az = totalX[i] + totalY[i] * I;
							Bz = totalX[j] + totalY[j] * I;
							z = (Pz - Az)/(Bz - Az);
							if ((creal(z)<0)||(creal(z)>1)) {  // projection out of the segment
								tempDA=distance(X[x], Y[y], totalX[i], totalY[i]);
								tempDB=distance(X[x], Y[y], totalX[j], totalY[j]);
								defD=min(tempDA, tempDB);
							} else {					   // projection inside the segment
								defD= cimag(z)*distance(totalX[j], totalY[j], totalX[i], totalY[i]); 
							}
							tempD=min(defD, tempD);
						} else tempD=min(defD, tempD);
					}
					PIPvector[n * (*nX) * (*nY) + x + y * (*nX)]= c;  		// fill the output values
					if (c==0) {
						DISTvector[n * (*nX) * (*nY) + x + y * (*nX)]=-tempD;	// fill the distance values
					} else {
						DISTvector[n * (*nX) * (*nY) + x + y * (*nX)]=tempD;
					}	// fill the distance values
				}
			}
		}
	}
}  */

/*
 * Function for calculating SIX TIMES the signed volume
 * of a tetrahedron in a triangular mesh
 */
double SignedVolumeOfTriangle(double p1X, double p1Y, double p1Z, 
		double p2X, double p2Y, double p2Z, double p3X, double p3Y, double p3Z) {
	double v321 = p3X*p2Y*p1Z;
	double v231 = p2X*p3Y*p1Z;
	double v312 = p3X*p1Y*p2Z;
	double v132 = p1X*p3Y*p2Z;
	double v213 = p2X*p1Y*p3Z;
	double v123 = p1X*p2Y*p3Z;
	//printf("\ndVolume=%lf", (double)(1.0/6.0));
	return (-v321 + v231 + v312 - v132 - v213 + v123);	
}

/*
 * Function for calculating the DOUBLE of area of a facet in a triangular mesh
 */
double FacetSurface(double p1X, double p1Y, double p1Z, 
		double p2X, double p2Y, double p2Z, double p3X, double p3Y, double p3Z) {
	double ax = p2X - p1X;
	double ay = p2Y - p1Y;
	double az = p2Z - p1Z;
	double bx = p3X - p1X;
	double by = p3Y - p1Y;
	double bz = p3Z - p1Z;
	double cx = ay*bz - az*by;
	double cy = az*bx - ax*bz;
	double cz = ax*by - ay*bx;
	//printf("\np2X*p3Y-p3X*p2Y=%lf",p2X * p3Y - p3X * p2Y );
	return sqrt(cx*cx + cy*cy + cz*cz);
}

/* 
 * Function for calulating the volume of a mesh
 */
void MeshVolume(double *X, double *Y, double *Z, int *numT, int *V1, int *V2, int *V3, double *Volume) {
	int n;			// counter
	*Volume=0; 		// initial volume
	for (n=0; n<*numT; n++) {
		*Volume = *Volume + SignedVolumeOfTriangle(X[V1[n]], Y[V1[n]], Z[V1[n]], X[V2[n]], Y[V2[n]], Z[V2[n]], X[V3[n]], Y[V3[n]], Z[V3[n]]);
		//printf("\nX = %lf Y = %lf Z = %lf n = %d", X[V1[n]], Y[V1[n]], Z[V1[n]], n);
	}
	*Volume = fabs(*Volume * (double)(1.0/6.0));  // absolute value of a double
}

/* 
 * Function for calculating the mesh surface
 */
void MeshSurface(double *X, double *Y, double *Z, int *numT, int *V1, int *V2, int *V3, double *Surface) {
	int n;			// counter
	*Surface = 0;	// initial surface
	for (n=0; n<*numT; n++) {
		*Surface = *Surface + FacetSurface(X[V1[n]], Y[V1[n]], Z[V1[n]], X[V2[n]], Y[V2[n]], Z[V2[n]], X[V3[n]], Y[V3[n]], Z[V3[n]]);
	}
	*Surface=0.5 * *Surface;
}

/*
 * Function to calculate the index position of an array which refers to a 3D matrix
 * x,y,z are the coors of the matrix you are interested in 
 * nx,ny,nz  are the dimensions of the matrix along the 3 axes
 *
 */
int posDecod(int x, int y, int z, int nx, int ny, int nz) {
  return( z*ny*nx+ y*nx + x   );
}
/*
 * Function to calculate the raw surface from a 3D voxel matrix
 + arr the input matrix, serialized
 * nX,nY,nZ  are the dimensions of the matrix along the 3 axes
 * pSX, pSY,pSZ are the pixelSpacing along the axec
 * surface is the output
 */

void rawSurface(double *arr, int *nX, int *nY, int *nZ, double *pSX, double *pSY, double *pSZ, double *surface) {
  
  int z,y,x;
  // reset surface value
  *surface = 0;
  // loop the 3D-matrix
  for( z = 0 ; z < *nZ ; z++ ) {
    for( y = 0 ; y < *nY ; y++ ) {
      for( x = 0 ; x < *nX ; x++ ) {
        
        if( arr[  posDecod(x,y,z,*nX,*nY,*nZ) ] !=0 ) {
          
          // if a non-zero voxel is on the border, surface cannot be calculated
          // for now skipped but you now... in the future....
  //        if( x==0 || x==*nX || y==0 || y==*nY || z==0 || z==*nZ) {*surface = -1; return; }
          
          // is it a border-voxel in each possible direction?
          if(z+1!=*nZ) {if( arr[  posDecod(x,y,z+1,*nX,*nY, *nZ ) ] == 0 ) *surface+= (*pSX) * (*pSY);}
          if(z-1>0) {if( arr[  posDecod(x,y,z-1,*nX,*nY, *nZ ) ] == 0 ) *surface+= (*pSX) * (*pSY);}
          if(y+1!=*nY) {if( arr[  posDecod(x,y+1,z,*nX, *nY, *nZ ) ] == 0 ) *surface+= (*pSX) * (*pSZ);}
          if(y-1>0) {if( arr[  posDecod(x,y-1,z,*nX,*nY, *nZ ) ] == 0 ) *surface+= (*pSX) * (*pSZ);}
          if(x+1!=*nX) {if( arr[  posDecod(x+1,y,z, *nX, *nY, *nZ ) ] == 0 ) *surface+= (*pSY) * (*pSZ);}
          if(x-1>0) {if( arr[  posDecod(x-1,y,z, *nX, *nY, *nZ ) ] == 0 ) *surface+= (*pSY) * (*pSZ);        }
          
          if(x-1==0) *surface+= (*pSY) * (*pSZ);
          if(x+1==*nX) *surface+= (*pSY) * (*pSZ); 
          if(y-1==0) *surface+= (*pSX) * (*pSZ);
          if(y+1==*nY) *surface+= (*pSX) * (*pSZ);          
          if(z-1==0) *surface+= (*pSX) * (*pSY);
          if(z+1==*nZ) *surface+= (*pSX) * (*pSY);           
            
        }
        
      }      
    }
  }
  return;
}

void erosion( double *cube, int *nX, int *nY, int *nZ, int *mx, int *my, int *mz, int *iterator) {
  int x,y,z,center,ct;
  if( *iterator >= 10) return;  // just to avoid infinite loops
  if(*mx == 0 && *my ==0 && *mz ==0 ) return;
  

  // loop per ogni elemento del cubo
  for( z=0; z<*nZ; z++ ) {
    for( y=0; y<*nY; y++ ) {
      for( x=0; x<*nX; x++) {
        // prendi l'offset relativo al punto in esame
        center = posDecod(x,y,z,*nX,*nY,*nZ);
        // se != 0 vediamo l'intorno
        if(cube[center]>0) {
          if( *mx>0 ){
            //if(x==0 || x==*nX) cube[center]=-1;
            if(x==0 || x>=(*nX-1)) cube[center]=-1;
            else if( cube[posDecod(x-1,y,z,*nX,*nY,*nZ)]==0 || cube[posDecod(x+1,y,z,*nX,*nY,*nZ)]==0  ) cube[center]=-1;
          }
          if( *my>0 ){
            // if(y==0 || y==*nY) cube[center]=-1;
            if(y==0 || y==(*nY-1)) cube[center]=-1;
            else if( cube[posDecod(x,y-1,z,*nX,*nY,*nZ)]==0 || cube[posDecod(x,y+1,z,*nX,*nY,*nZ)]==0  ) cube[center]=-1;
          }        
          if( *mz>0 ){
            // if(z==0 || z==*nZ) cube[center]=-1;
            if(z==0 || z==(*nZ-1)) cube[center]=-1;
            else if( cube[posDecod(x,y,z-1,*nX,*nY,*nZ)]==0 || cube[posDecod(x,y,z+1,*nX,*nY,*nZ)]==0  )  cube[center]=-1;
          }
        }
      }
    }
}
  // decrementa i vincoli sui margini
  if( *mx>0 ) *mx = *mx - 1;
  if( *my>0 ) *my = *my - 1;  
  if( *mz>0 ) *mz = *mz - 1;
  // trasforma tutti i '-1' in '0' per l'iterazione successiva
  //for(ct=0; ct<= posDecod((*nX-1),(*nY-1),(*nZ-1),*nX,*nY,*nZ); ct++) {
  for(ct=0; ct<= ((*nX)*(*nY)*(*nZ)-1); ct++) {
    if(cube[ct]==-1) cube[ct]=0;
  }
  // rilancia ricorsivamente
  *iterator = *iterator + 1;
  erosion( cube, nX, nY, nZ, mx, my, mz, iterator );
}
/*
 * give back row, cloumn, slice from cIndec (index of aa linearized matrix)
 */
void RCSfromC( int cIndex, struct pointInSpace *point, int NX, int NY, int NZ ) {
  int x,y,z;
  z = (int)(cIndex / (NY*NX));
  y = (int)((cIndex - z*(NY*NX) )/NX);
  x = cIndex - NX*y - z*NX*NY;
  point->x = x;
  point->y = y;
  point->z = z;
//  printf("\n punto calcolato in cIndex = %d RCS: %d %d %d, pxyz=%d,%d,%d (nxnynz=%d,%d,%d)",cIndex,x,y,z,point->x,point->y,point->z,(int)NX,(int)NY,(int)NZ);
  return;
}
void getNxNyFrom3D(double Px, double Py, 
                   double a11, double a21, double a31, double a12, 
                   double a22, double a32, double a14, double a24, struct cube3dcoord *coord) {
  coord->y = (a21*(Px-a14)+a11*(a24-Py)) / (a12*a21-a22*a11 ) ;
  coord->x = (a22*(Px-a14)+a12*(a24-Py)) / (a11*a22-a21*a12 ) ;
  return;
}
double trilineareGeneralizzata(
    double x0y0z0, double x0y0z1, double x0y1z0, 
    double x0y1z1, double x1y0z0, double x1y0z1, 
    double x1y1z0, double x1y1z1, 
    double x0,double y0, double z0, 
    double dx1x0, double dy1y0, double dz1z0, 
    double x, double y,double z) {
  
  double xd,yd,zd,c00,c01,c10,c11,c0,c1,c;
  xd = (x-x0)/dx1x0;
  yd = (y-y0)/dy1y0;
  zd = (z-z0)/dz1z0;
  c00 = x0y0z0*(1-xd)+x1y0z0*xd;
  c10 = x0y1z0*(1-xd)+x1y1z0*xd;
  c01 = x0y0z1*(1-xd)+x1y0z1*xd;
  c11 = x0y1z1*(1-xd)+x1y1z1*xd;
  c0 = c00*(1-yd)+c10*yd;
  c1 = c01*(1-yd)+c11*yd;
  c = c0*(1-zd)+c1*zd;
  return c;
}
/*
 * NAME: interpolaPerPuntiSparsi
 * DESCRIZIONE: calcola le coordinate di ogni voxel di una matrice di voxel 
 *              partendo da i DOM delle varie slice.
 *              
 *  NX,NY,NZ: the number of row, clumns, slices
 *  directionCosinesArray: the array with the direction cosines, uno per ogni slice (9 per slice)
 *  xNew,yNew,zNew: le coordinate dei punti x,y,z di cui voglio calcolare l'interpolata
 *  numberOfInputPoints : numero dei punti da interpolare
 *  interpolatedValues :  array contenente i valori interpolati (output)
 *  originalValues : array contenente la matrice originale (input)
 */
void interpolaPerPuntiSparsi(
    int *NX, int *NY,int *NZ,
    double *directionCosinesArray,  
    double *xNew, double *yNew, double *zNew, 
    int *numberOfInputPoints,
    double *interpolatedValues, double *originalValues)  {
  
  int x,y,z,ct,punto;
  struct pointInSpace point;
  struct cube3dcoord cubeCoord;
  struct DICOM_OrientationMatrix DOM;  
  double *xOld, *yOld, *zOld;
  int *xOldN, *yOldN, *zOldN;  // to remove
  double distanza;
  int numeroDiVoxel,x0,x1,y0,y1,z0,z1;
  double valoreCalcolato;
  struct DICOM_OrientationMatrix *DOMStruct;
  
  numeroDiVoxel = (*NX) * (*NY) * (*NZ) ;
  
  xOld = (double*)calloc( numeroDiVoxel , sizeof(double));
  yOld = (double*)calloc( numeroDiVoxel , sizeof(double));
  zOld = (double*)calloc( numeroDiVoxel , sizeof(double));
  
  xOldN = (int*)calloc( numeroDiVoxel , sizeof(int)); // to remove
  yOldN = (int*)calloc( numeroDiVoxel , sizeof(int)); // to remove
  zOldN = (int*)calloc( numeroDiVoxel , sizeof(int)); // to remove
  
  DOMStruct = (struct DICOM_OrientationMatrix*)calloc((*NZ), sizeof(struct DICOM_OrientationMatrix));
  
  if( xOld == NULL || yOld==NULL || zOld==NULL || xOldN==NULL || yOldN==NULL || zOldN==NULL) {
    printf("\n Error in allocating memory");
    return;
  }
  
  // ------------------------------------
  // calcola gli array con le posizioni per ogni punto della matrice di input
  // (della matrice interpolante)
  // ------------------------------------
  ct=0;
  for( z = 0; z < (*NZ); z++ ) {
    for( y = 0; y < (*NY); y++ ) {
      for( x = 0; x < (*NX); x++ ) {
        DOM.a11= directionCosinesArray[z * 9 + 0];	// position 0 in array
        DOM.a21= directionCosinesArray[z * 9 + 1];
        DOM.a31= directionCosinesArray[z * 9 + 2];
        DOM.a12= directionCosinesArray[z * 9 + 3];
        DOM.a22= directionCosinesArray[z * 9 + 4];
        DOM.a32= directionCosinesArray[z * 9 + 5];
        DOM.Sx=  directionCosinesArray[z * 9 + 6];
        DOM.Sy=  directionCosinesArray[z * 9 + 7];
        DOM.Sz=  directionCosinesArray[z * 9 + 8];
        DOM.XpixelSpacing=1;    DOM.YpixelSpacing=1;  // lo immagino già nella DOM
        point = get3DPosFromNxNy(x, y, DOM);
        xOld[ct] = point.x;         yOld[ct] = point.y;         zOld[ct] = point.z;
        xOldN[ct] = x;         yOldN[ct] = y;         zOldN[ct] = z;
        DOMStruct[z] = DOM;
        ct++;
      }
    }
  }

  // ------------------------------------
  // ora scorri l'array dei punti da interpolare e cerca i vertici adiacenti
  // ------------------------------------
  // per ogni punto da interpolare
  for(punto=0; punto < *numberOfInputPoints; punto++ ) {
    double deltaMinore = 1000000000;
    double deltaDOMMinore = 1000000000;
    int posizioneMinore=0;    
    int posizioneDOMMinore = 0;
    // prendi le coordinate riga,colonna,slide del punto da interpolare
    // scorri i DOM, prendi quello con la Z più vicina
    for(z=0; z<(*NZ); z++) {
      if( fabs(DOMStruct[z].Sz - zNew[punto]) < deltaDOMMinore  ) {
        deltaDOMMinore = fabs(DOMStruct[z].Sz - zNew[punto]) ; 
        posizioneDOMMinore = z;
      }
    }
    // prendi le coordinate riga,colonna,slide del punto da interpolare
    getNxNyFrom3D(xNew[punto], yNew[punto], 
                  DOMStruct[posizioneDOMMinore].a11, DOMStruct[posizioneDOMMinore].a21, 
                  DOMStruct[posizioneDOMMinore].a31, DOMStruct[posizioneDOMMinore].a12, 
                  DOMStruct[posizioneDOMMinore].a22, DOMStruct[posizioneDOMMinore].a32, 
                  DOMStruct[posizioneDOMMinore].Sx, DOMStruct[posizioneDOMMinore].Sy, 
                  &cubeCoord);
    // guarda in un intorno del punto identificato per cercare il punto più vicino.
    // assumo che la coordinata z sia corretta (l'ho calcolata prima)
    int arrIndici[9];
    int counter;
    
    arrIndici[0] = posDecod(cubeCoord.x-1, cubeCoord.y-1, posizioneDOMMinore, *NX, *NY, *NZ);
    arrIndici[1] = posDecod(cubeCoord.x,   cubeCoord.y-1, posizioneDOMMinore, *NX, *NY, *NZ);
    arrIndici[2] = posDecod(cubeCoord.x+1, cubeCoord.y-1, posizioneDOMMinore, *NX, *NY, *NZ);
    arrIndici[3] = posDecod(cubeCoord.x-1, cubeCoord.y,   posizioneDOMMinore, *NX, *NY, *NZ);
    arrIndici[4] = posDecod(cubeCoord.x,   cubeCoord.y,   posizioneDOMMinore, *NX, *NY, *NZ);
    arrIndici[5] = posDecod(cubeCoord.x+1, cubeCoord.y,   posizioneDOMMinore, *NX, *NY, *NZ);
    arrIndici[6] = posDecod(cubeCoord.x-1, cubeCoord.y+1, posizioneDOMMinore, *NX, *NY, *NZ);
    arrIndici[7] = posDecod(cubeCoord.x,   cubeCoord.y+1, posizioneDOMMinore, *NX, *NY, *NZ);
    arrIndici[8] = posDecod(cubeCoord.x+1, cubeCoord.y+1, posizioneDOMMinore, *NX, *NY, *NZ);

    for(counter=0; counter<9; counter++) {
      ct = arrIndici[counter];
      distanza = sqrt((xOld[ct] - xNew[punto])*(xOld[ct] - xNew[punto]) +
        (yOld[ct] - yNew[punto])*(yOld[ct] - yNew[punto]) +
        (zOld[ct] - zNew[punto])*(zOld[ct] - zNew[punto]));
      if( distanza <  deltaMinore ) {  
        deltaMinore = distanza;  
        posizioneMinore = ct;   
      }
    }
    /*
    // cerchiamo il punto, nella matrice interpolante,
    // più vicino in politica di distanza euclidea
    for(  ct=0 ; ct < numeroDiVoxel ; ct++  ) {
      distanza = sqrt((xOld[ct] - xNew[punto])*(xOld[ct] - xNew[punto]) +
              (yOld[ct] - yNew[punto])*(yOld[ct] - yNew[punto]) +
              (zOld[ct] - zNew[punto])*(zOld[ct] - zNew[punto]));
      if( distanza <  deltaMinore ) {  
        deltaMinore = distanza;  
        posizioneMinore = ct;   
      }
    }
    */
//    printf("\n\n  nx=%d ny=%d nz=%d",cubeCoord.x,cubeCoord.y,posizioneDOMMinore);
//    printf("\n\n  nxOld=%d nyOld=%d nzOld=%d",xOldN[posizioneMinore],yOldN[posizioneMinore],zOldN[posizioneMinore]);
//    return;
    // ora in 'posizionMinore' c'è l'indice corrispondente al punto della matrice interpolante più
    // prossimo al punto 'punto' da analizzare
    
    // prendi la posizione riga, colonna, slide del punto più vicino corrispondente
    RCSfromC( posizioneMinore, &point, *NX, *NY, *NZ );
    
    // ora cerca i vertici adiacenti
    int cx0,cx1,cy0,cy1,cz0,cz1;
    cx0 = posDecod(point.x-1, point.y, point.z, *NX, *NY, *NZ);
    cx1 = posDecod(point.x+1, point.y, point.z, *NX, *NY, *NZ);
    cy0 = posDecod(point.x, point.y-1, point.z, *NX, *NY, *NZ);
    cy1 = posDecod(point.x, point.y+1, point.z, *NX, *NY, *NZ);
    cz0 = posDecod(point.x, point.y, point.z-1, *NX, *NY, *NZ);
    cz1 = posDecod(point.x, point.y, point.z+1, *NX, *NY, *NZ);

    double x00, x01, y00, y01, z00, z01;
    int x00Adder,x01Adder,y00Adder,y01Adder,z00Adder,z01Adder;
    
    if( fabs(xNew[punto]-xOld[cx0])-fabs(xNew[punto]-xOld[cx1]) <0 ) {
      x00 = xOld[cx0]; x01 = xOld[posizioneMinore]; 
      x00Adder = -1;  x01Adder = 0;
      }  
    else { 
      x00 = xOld[posizioneMinore];  x01=xOld[cx1]; 
      x00Adder = 0;  x01Adder = +1;
      }
    if( fabs(yNew[punto]-yOld[cy0])-fabs(yNew[punto]-yOld[cy1]) <0) {
      y00 = yOld[cy0]; y01 = yOld[posizioneMinore];
      y00Adder = -1; y01Adder = 0;
      }  
    else { 
      y00 = yOld[posizioneMinore];  y01=yOld[cy1]; 
      y00Adder = 0; y01Adder = +1;
      }
    if( fabs(zNew[punto]-zOld[cz0])-fabs(zNew[punto]-zOld[cz1]) <0 ) {
      z00 = zOld[cz0]; z01 = zOld[posizioneMinore];
      z00Adder = -1; z01Adder = +1;
      }  
    else { 
      z00 = zOld[posizioneMinore];  z01=zOld[cz1]; 
      z00Adder = 0; z01Adder = +1;
      }    
    
    int v000,v001,v010,v011,v100,v101,v110,v111;
    v000 = posDecod(point.x + x00Adder, point.y + y00Adder, point.z + z00Adder, *NX, *NY, *NZ);
    v001 = posDecod(point.x + x00Adder, point.y + y00Adder, point.z + z01Adder, *NX, *NY, *NZ);
    v010 = posDecod(point.x + x00Adder, point.y + y01Adder, point.z + z00Adder, *NX, *NY, *NZ);
    v011 = posDecod(point.x + x00Adder, point.y + y01Adder, point.z + z01Adder, *NX, *NY, *NZ);
    v100 = posDecod(point.x + x01Adder, point.y + y00Adder, point.z + z00Adder, *NX, *NY, *NZ);
    v101 = posDecod(point.x + x01Adder, point.y + y00Adder, point.z + z01Adder, *NX, *NY, *NZ);
    v110 = posDecod(point.x + x01Adder, point.y + y01Adder, point.z + z00Adder, *NX, *NY, *NZ);
    v111 = posDecod(point.x + x01Adder, point.y + y01Adder, point.z + z01Adder, *NX, *NY, *NZ);
    
    
        
    //    printf("\n x::  %lf < %lf < %lf scelto:<%lf, %lf>",xOld[cx0],xNew[punto],xOld[cx1],x00,x01);
    //    printf("\n y::  %lf < %lf < %lf scelto:<%lf, %lf>",yOld[cy0],yNew[punto],yOld[cy1],y00,y01);
    //    printf("\n z::  %lf < %lf < %lf scelto:<%lf, %lf>",zOld[cz0],zNew[punto],zOld[cz1],z00,z01);
    if(  fabs(x01-x00)>=5 || fabs(y01-y00)>=5 || fabs(z00-z01)>=5  ) {
      printf("\n\n punto da interpolare = %lf,%lf,%lf",xNew[punto],yNew[punto],zNew[punto]);
      printf("\n x00,x01=%lf,%lf  y00,y01=%lf,%lf   z00,z01=%lf,%lf",x00,x01,y00,y01,z00,z01);
      return;
    }

    
    valoreCalcolato = trilineareGeneralizzata(
      originalValues[v000],  //x0y0z0 (sample value)
      originalValues[v001],  //x0y0z1 (sample value)
      originalValues[v010],  //x0y1z0 (sample value)
      originalValues[v011],  //x0y1z1 (sample value)
      originalValues[v100],  //x1y0z0 (sample value)
      originalValues[v101],  //x1y0z1 (sample value)
      originalValues[v110],  //x1y1z0 (sample value)
      originalValues[v111],  //x1y1z1 (sample value)
      xOld[v000], //x0
      yOld[v000], //y0,
      zOld[v000], //z0,
      xOld[v000]-xOld[v100], //dx1x0,
      yOld[v000]-xOld[v010], //dy1y0,
      zOld[v000]-xOld[v001], //dz1z0,
      xNew[punto], yNew[punto], zNew[punto]);
              
              
    interpolatedValues[punto] = valoreCalcolato;              
                                       
/*
    printf("\n\n Punto da interpolare: <%lf,%lf,%lf>",xNew[punto], yNew[punto], zNew[punto]); 
    printf("\n (000=<%lf,%lf,%lf>,001=<%lf,%lf,%lf>, \n010=<%lf,%lf,%lf>,011=<%lf,%lf,%lf>, \n100=<%lf,%lf,%lf>,101=<%lf,%lf,%lf>, \n110=<%lf,%lf,%lf>,111=<%lf,%lf,%lf>) ",
          xOld[v000],yOld[v000],zOld[v000],xOld[v001],yOld[v001],zOld[v001],
          xOld[v010],yOld[v010],zOld[v010],xOld[v011],yOld[v011],zOld[v011],
          xOld[v100],yOld[v100],zOld[v100],xOld[v101],yOld[v101],zOld[v101],
          xOld[v110],yOld[v110],zOld[v110],xOld[v111],yOld[v111],zOld[v111]
    );
    printf("\n\n (000=%lf,001=%lf,010=%lf,011=%lf,100=%lf,101=%lf,110=%lf,111=%lf) ",originalValues[v000],originalValues[v001],originalValues[v010],originalValues[v011],originalValues[v100],originalValues[v101],originalValues[v110],originalValues[v111]);  
    printf("\n\n x0=%lf  y0=%lf   z0=%lf  ",xOld[v000],yOld[v000],zOld[v000] ); 
    printf("\n\n dx1x0=%lf   dy1y0=%lf   dz1z0=%lf",xOld[v000]-xOld[v100],yOld[v000]-yOld[v010],zOld[v000]-zOld[v001]);
    printf("\n\n VALORE = %lf",valoreCalcolato);
 */     

  }
  free(xOld);  free(yOld);  free(zOld);
  free(xOldN);  free(yOldN);  free(zOldN);
  return;
}
