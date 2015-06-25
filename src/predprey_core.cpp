// 
// 
// Function to compute core dynamics of pred/prey model
// 
// 

# include <cstdlib> // for rand()
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
  Rcout << grid(1,3);
  
  update = 0;
  // Default value for subs is 1000 in original model
  while (update < subs) { 
    Rcout << "plop\n";
    
    // Choose a random cell
    i = randn(0, h);
    j = randn(0, w);
    Rcout << "i:" << i << "j:" << j << "\n";
    
    if ( grid(i,j) == FISH ) { // if fish, try to grow 
      Rcout << "plap\n";
      // consider random neighbor
      int k  = randn(0, NEIGHBORS); 
      int i1 = (X[k] + i) % h; // wrap around
      int j1 = (Y[k] + j) % w;
      
      if ( grid(i1, j1) == EMPTY && randp() < betaf ) { 
        grid(i1, j1) = FISH; // grow
      }
      
    } else if ( grid(i,j) == SHARK ) { // if shark, then eat or die, and maybe grow
      
      // consider neighbors starting randomly
      randnb = randn(0,3);
      k = 0;
      found_prey = false;
      while ( !found_prey && k < NEIGHBORS) { 
        Rcout << "plip:" << (k + randnb) % NEIGHBORS << "\n";
        i1 = (X[(k + randnb) % NEIGHBORS] + i) % w;
        j1 = (Y[(k + randnb) % NEIGHBORS] + j) % h;
        if ( grid(i1, j1) == FISH ) { 
          // eat
          found_prey = true;
          grid(i1, j1) = SHARK; 
          // and maybe reproduce 
          if ( ( 1 - randp()) < betas) { 
            grid(i, j) = EMPTY; // otherwise grid(i,j) has shark = repro occured
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
      i1 = (i + X[k]) % h;
      j1 = (j + Y[k]) % w;
      Rcout << "plap. i1:" << i1 << " j1:" << j1 << "\n";
      
      if ( grid(i, j) != grid(i1, j1) ) { 
        tmp = grid(i, j);
        grid(i1, j1) = grid(i, j);
        grid(i, j) = tmp;
      }
    }
    
    update++;
  }
  
  Rcout << "returning\n";
  
  return grid;
}