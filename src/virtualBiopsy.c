#include <stdlib.h>  
#include <stdio.h>

void new_virtualBiopsy(   double *esame, int *dim1, int *dim2, int *dim3, int *stepX, int *stepY, 
                   int *stepZ, int *lunghezza, int *esameCarot, int *uni    ){
  
  int i;      int j;    int k;
  int indice;   int indice2;
  int m;      int n;    int o;
  int h;      int indsomma;
  
  h=0;
  //preso il primo dx, dy, dz utile applico il carotaggio di quella data quantita' all'array paziente
  for ( k=*stepZ;  k<((*dim3)-(*stepZ));  k++ )
  {
    for(  j=*stepY;  j<((*dim2)-(*stepY));  j++  )
    {
      for(  i=*stepX;   i<((*dim1)-(*stepX));   i++  )
      {
        indice= ((k*(*dim1)*(*dim2))+(j*(*dim1))+i);
        
        // se il valore analizzato e' diverso da zero entra in questo ciclo e fa tutti expand grid 
        //nell'intorno del punto      
        if(   esame[  indice   ]  !=  0   )
        { 
          
          indsomma=0;
          
          for (m=-*stepX; m<=*stepX; m++) 
          {
            for (n=-*stepY; n<=*stepY; n++) 
            {
              for (o=-*stepZ; o<=*stepZ; o++) 
              {
                indice2 = ((k+o)*(*dim1)*(*dim2))  +  ((j+n)*(*dim1))  +  (i+m);
                if ( esame   [  indice2    ] !=0) 
                {
                  //per ogni valore dell'intorno diverso da zero incrementa la variabile indsomma
                  indsomma=indsomma + 1;
                }
              }
            }
          }
          // se indsomma e' uguale al numero di check fatti incrementa di 1 la variabile uni e in corrispondenza
          // dell'indice del centroide sovrascrive 1
          if (indsomma == *lunghezza) 
          {
            esameCarot[indice] = 1;
            *uni = *uni+1;
          }
        }
      }
    }
  }
}

void virtualBiopsy(   double *esame, int *dimensione1, int *dimensione2, int *dimensione3, int *nx, int *ny, int *nz,
                      int *expand, int *controllo  ,int *lunghezza, int *carotaggioVolume    ){
  int i;
  int j;
  int k;
  int indice;
  int m;
  int n;
  int o;
  int indsomma;

  
  for ( k=*nz;  k<*dimensione3;  k++ )
  {
    for(  j=*ny;  j<*dimensione2;  j++  )
    {
      for(  i=*nx;   i<*dimensione1;   i++  )
      
      {
        indice= (k*(*dimensione1)*(*dimensione2))+(j*(*dimensione1))+i;
              
        
        if(   esame[  indice   ]  !=  0   )
        {
          indsomma=0;
          for (m=-*nx; m<=*nx; m++) {
            for (n=-*ny; n<=*ny; n++) {
              for (o=-*nz; o<=* nz; o++) {
                if ( esame   [  ((k+o)*(*dimensione1)*(*dimensione2))  +  ((j+n)*(*dimensione1))  +  (i+m)    ] !=0) {
                  indsomma=indsomma + 1;
                 
                }
              }
            }
          }
        }
        if (indsomma == *controllo) {
          carotaggioVolume[indice] = 1;
          
        }
      }
    }
  }
  
}

// C function to evaluate grey matrix to 1/0 matrix
void greyMatrix(int *indX, int *indY, int *indZ, int *numCentroidi, int *nx, int *ny, int *nz,
                int *outputGrigi, double *exam, int *size1, int *size2, int *size3, int *control  ){               
  int ct;
  int l;
  int m;
  int n;
  int numcarot=0;
  int indice2;
  // examine only index of the point not equal to zero
  for(ct=0; ct<*numCentroidi; ct++)
  {
    // expand grid around point not equal to zero
    for (n=-*nz; n<=*nz; n++)
    {
      for (m=-*ny; m<=*ny; m++)
      {
        for (l=-*nx; l<=*nx; l++)
        {
          indice2= ( ((((indZ[ct])-1)+n)  *  (*size1)*(*size2))  +  ((((indY[ct])-1)+m)*(*size1)) 
                                                +  (  ((indX[ct])-1)  +l   ));
          // create grey matrix
          outputGrigi[numcarot]= exam   [ indice2 ];
          numcarot++;
        }
      }   
    } 
  }
}