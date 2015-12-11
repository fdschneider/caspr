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
  return floor( min + result[0] * ((max+1) - min) );
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
                            double betaf,           // growth rate of fish
                            double betas,           // growth rate of shark
                            double delta,           // death  rate of shark
                            double m               // density indep death rate of fish
                            ) { 
  int i, j, h, w, update, k, i1, j1, randnb, rep, tmp;
  bool found_prey;
  
  // Define states explicitely for clarity
  int FISH  = 1;
  int SHARK = 2;
  int EMPTY = 3;
  int NEIGHBORS = 4; // 4 the max number of neighbors
  
  // Adjust probabilities so they match to the number of substeps
  double betas_ = betas / subs;
  double betaf_ = betaf / subs;
  double delta_ = delta / subs;
  double m_     = m / subs;
  
  // Vector to sample a random value between 1 and -1
  IntegerVector X = IntegerVector::create(1,  0, -1,  0);
  IntegerVector Y = IntegerVector::create(0,  1,  0, -1);
  
  // Grid height/width
  h = grid.nrow();
  w = grid.ncol();
  
  // Loop over sub-iterations
  for ( int sub = 0; sub < subs; sub++) { 
    
    // We give a change to each cell to transition state
    update = 0;
    while (update < h*w ) { 
      
      // Choose a random cell
      i = randn(0, h);
      j = randn(0, w);
      
      if ( grid(i,j) == FISH ) { // if fish, try to grow 
        
        double rtest = randp();
        
        if ( rtest < betaf_ ) {
          
          // consider random neighbor
          int k  = randn(0, NEIGHBORS); 
          // The + h/w makes sure i1 or j1 is always positive
          int i1 = (X[k] + i + h) % h; // wrap around
          int j1 = (Y[k] + j + w) % w; 
          
          if ( grid(i1, j1) == EMPTY ) {
            grid(i1, j1) = FISH; // grow
          }
          
        } else if ( rtest > (1 / subs) - m_ ) {
          grid(i, j) = EMPTY;
        }
        
      // if shark, then eat or die, and maybe grow
      } else if ( grid(i,j) == SHARK ) { 
        
        // consider neighbors starting with a random one (0, 1, 2 or 3)
        randnb = randn(0, NEIGHBORS - 1);
        k = 0;
        found_prey = false;
        while ( !found_prey && k < NEIGHBORS) { 
          // Rcout << "plip:" << (k + randnb) % NEIGHBORS << "\n";
          i1 = (X[(k + randnb) % NEIGHBORS] + i + h) % h;
          j1 = (Y[(k + randnb) % NEIGHBORS] + j + w) % w;
          if ( grid(i1, j1) == FISH ) { 
            // eat
            found_prey = true;
            // Note that the probability of eating is "one" but we need to 
            //   divide it by subs
            if ( randp() < 1 / subs ) { 
              grid(i1, j1) = SHARK; 
              // and maybe reproduce 
              if ( randp() >= betas_ ) { 
    //             Rcout << "shark reproduces\n";
                grid(i, j) = EMPTY;
              } 
            }
          }
          k++; // consider next neighbor
        }
        // pred has not eaten -> he can die of starvation
  //      Rcout << "randp: " << randp() << "delta: " << delta;
        if ( !found_prey && (randp() < delta_) ) { 
  //         Rcout << "shark dies\n";
          grid(i, j) = EMPTY;
        }
      }
      
      // Random pair-wise mixing
      // Note: we do this for two pair of cells as in the original code, but this
      //   is not mentioned in Pascual2002.
      // Note2: [!!] This mixes two pairs of cell: 
      if ( randp() < 1 / subs ) {  
//        for (rep=0;rep<1;rep++) { 
          i  = randn(0, h);
          j  = randn(0, w);
          k  = randn(0, NEIGHBORS - 1); // (0, 1, 2 or 3)
          i1 = (i + X[k] + h) % h; 
          j1 = (j + Y[k] + w) % w;
          
          if ( grid(i, j) != grid(i1, j1) ) { 
            tmp = grid(i, j);
            grid(i1, j1) = grid(i, j);
            grid(i, j) = tmp;
          }
//        }
      }
      
      update++;
    }
  }
  
  return grid;
}