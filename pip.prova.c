#include <stdio.h>
#include <stdlib.h>

void pointIn3DPolygon( double *xMesh, double *yMesh, double *zMesh, int *howMany, double *x, double *y, double *z, int *result) {

  double error=0.1;
  int ct;
  int iPoint=0;

	for(ct=0; ct<*howMany; ct++) {

		if( (zMesh[ct]-*z)*(zMesh[ct]-*z) < error ) {
			iPoint++;
		}


	}

  printf("\n Interesting Points:%d",iPoint);

  *result = iPoint;


  return;
}
