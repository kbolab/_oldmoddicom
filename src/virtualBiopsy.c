#include <stdlib.h>  
#include <stdio.h>

void virtualBiopsy(   int *esame, int *dimensione1, int *dimensione2, int *dimensione3, int *nx, int *ny, int *nz,
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
