// 
// 
// Function to compute core dynamics of pred/prey model
// 
// 

#include <Rcpp.h>
using namespace Rcpp;


// Generate random integer between min and max
// [[Rcpp::export]] // for debug purposes only
int randn(double min, double max) { 
  NumericVector result = runif(1, 0, 1);
  return floor( min + result[0] * (max - min) );
}

// Generate random number between 0 and 1 as scalar using runif function in 
// Rcpp. Chapuza.
// [[Rcpp::export]] // for debug purposes only
double randp() { 
  NumericVector result = runif(1,0,1);
  return result[0];
}

// [[Rcpp::export]]
IntegerMatrix predprey_core(IntegerMatrix grid, // grid matrix
                            int subs,               // number of updates
                            int NEIGHBORS,          // number of neighbors
                            double betaf,           // growth rate of fish
                            double betas,           // growth rate of shark
                            double delta           // death  rate of shark
                            ) { 
  int i, j, h, w, update, k, i1, j1, randnb, rep, tmp;
  bool found_prey;
  
  // Define states explicitely for clarity
  int FISH  = 0;
  int SHARK = 1;
  int EMPTY = 2;

  // Vector to sample a random value between 1 and -1
  IntegerVector X = IntegerVector::create(1,  0, -1,  0);
  IntegerVector Y = IntegerVector::create(0,  1,  0, -1);
  
  // Grid height/width
  h = grid.nrow();
  w = grid.ncol();
  
  update = 0;
  // Default value for subs is 1000 in original model
  while (update < subs) { 
    
    // Choose a random cell
    i = randn(0, h);
    j = randn(0, w);
    
    if ( grid(i,j) == FISH ) { // if fish, try to grow 
      // consider random neighbor
      int k  = randn(0, NEIGHBORS); 
      // The + h/w makes sure i1 or j1 is always positive
      int i1 = (X[k] + i + h) % h; // wrap around
      int j1 = (Y[k] + j + w) % w; 
      
      if ( grid(i1, j1) == EMPTY && randp() < betaf ) { 
        grid(i1, j1) = FISH; // grow
      }
      
    } else if ( grid(i,j) == SHARK ) { // if shark, then eat or die, and maybe grow
      
      // consider neighbors starting with a random one
      randnb = randn(0,3);
      k = 0;
      found_prey = false;
      while ( !found_prey && k < NEIGHBORS) { 
        // Rcout << "plip:" << (k + randnb) % NEIGHBORS << "\n";
        i1 = (X[(k + randnb) % NEIGHBORS] + i + h) % h;
        j1 = (Y[(k + randnb) % NEIGHBORS] + j + w) % w;
        if ( grid(i1, j1) == FISH ) { 
          // eat
          found_prey = true;
          grid(i1, j1) = SHARK; 
          // and maybe reproduce 
          if ( randp() >= betas ) { 
            grid(i, j) = EMPTY;
          } 
        }
        k++; // consider next neighbor
      }
      // pred has not eaten -> he can die of starvation
      if ( !found_prey && (randp() < delta) ) { 
        grid(i, j) = EMPTY;
      }
    }
    
    // Random pair-wise mixing
    // Note: we do this for two pair of cells as in the original code, but this
    //   is not mentioned in Pascual2002.
    for (rep=0;rep<2;rep++) { 
      i  = randn(0, h);
      j  = randn(0, w);
      k  = randn(0, NEIGHBORS);
      i1 = (i + X[k] + h) % h; 
      j1 = (j + Y[k] + w) % w;
      
      if ( grid(i, j) != grid(i1, j1) ) { 
        tmp = grid(i, j);
        grid(i1, j1) = grid(i, j);
        grid(i, j) = tmp;
      }
    }
    
    update++;
  }
  
  
  return grid;
}