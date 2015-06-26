  // PREDATOR-PREY MODEL
// PRINTS DELTA (SHARK DEATH RATE) VS. COVER (FISH DENSITY)

#include <iostream>
#include "../../util/ran2_NR.h"
using namespace std;

/**** globals ****/
const int   EMPTY     = 0;
const int   FISH      = 1;
const int   SHARK     = 2;
const int   NEIGHBORS = 4;
const int   W         = 100;
const int   NP        = W*W;
const int   MXT       = 1000;
const int   RMX       = 200;
const float ALPHA1    = 0.5;	// fish growth rate
const float PF        = 0.4;
const float PS        = 0.1;
const int   X[]       = {1,  0, -1,  0};
const int   Y[]       = {0,  1,  0, -1};

/**** main routine ****/
main(){
   int i, j, k, kk, l, i1, j1, count, temp, rep, row[4], col[4], **grid, t, n;
   float r, delta, alpha2, ntot, cover; //alpha2/delta => shark growth/death rate
   seed();
   grid = new int* [W];
   for(i=0;i<W;i++) grid[i] = new int [W];

   for(delta=0.95;delta>=0.1;delta-=0.1){
      alpha2 = 1.0 - delta;
      
      for(n=1;n<=RMX;n++){

         // initialize grid
         for(i=0;i<W;i++)
            for(j=0;j<W;j++){
	       r = randm();
	       if(r<PF) 	 grid[i][j] = FISH;
	       else if(r>=1.-PS) grid[i][j] = SHARK;
	       else	         grid[i][j] = EMPTY;
            }

         for(t=1;t<=MXT;t++){
            for(l=0;l<NP;l++){		// asynchronous updating begins
	       i = int(W*randm());
	       j = int(W*randm());

      	       if(grid[i][j]==FISH){		// if FISH, grow
	          k = int(NEIGHBORS*randm());
                  i1 = (i + X[k] + W) % W;
                  j1 = (j + Y[k] + W) % W;
                  if(grid[i1][j1]==EMPTY && randm()<ALPHA1) 
	    		grid[i1][j1] = FISH;
      	       }
               else if(grid[i][j]==SHARK){	// if SHARK, eat+grow/die
                  count = 0;
                  for(k=0;k<NEIGHBORS;k++){
                     i1 = (i + X[k] + W) % W;
                     j1 = (j + Y[k] + W) % W;
                     if(grid[i1][j1]==FISH){
                        row[count] = i1; col[count] = j1; count++;
                     }
                  }
                  if(count>0){
                     k = int(count*randm());
                     i1 = row[k]; j1 = col[k]; grid[i1][j1] = SHARK;
                     if(randm()>=alpha2) grid[i][j] = EMPTY;
                  }
                  else if(randm()<delta) grid[i][j] = EMPTY;
               }

               /****  randm pair-wise mixing ****/
               for(rep=0;rep<2;rep++){
	          i = int(W*randm());
	          j = int(W*randm());
                  k = int(NEIGHBORS*randm());
                  i1 = (i + X[k] + W) % W;
                  j1 = (j + Y[k] + W) % W;
                  if(grid[i][j]!=grid[i1][j1]){
                     temp = grid[i][j];
                     grid[i][j] = grid[i1][j1];
                     grid[i1][j1] = temp;
                  }
               }
            }
         }

	 // computing cover
	 ntot = 0;
	 for(i=0;i<W;i++)
	    for(j=0;j<W;j++) if(grid[i][j]==FISH) ntot++;
	 cover += ntot/NP;;
      }
      cout << "delta = " << delta << ",  cover = " << cover/RMX << "\n";
   }

   return 0;
}

