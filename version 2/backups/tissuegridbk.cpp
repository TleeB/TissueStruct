#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
//must add flags -I /opt/local/include
//using the Boost library
#include <boost/random.hpp>
#include <boost/format.hpp>

//Vectors to hold neighbors
#include <vector>

using namespace std;
int length_myo = 5;
int length_fib = 1;
const int NROWS = 5, NCOLS = 125;
int grid[NROWS][NCOLS]={0};


int main(){//Forgot the double parenthesis
  FILE *grid_file;
  FILE *neighbors_myo_file;
  FILE *neighbors_fib_file;
  FILE *stim_myo_file;
  FILE *stim_fib_file;

  grid_file = fopen("grid_file.dat", "w");
  neighbors_myo_file = fopen("neighbors_myo_file.dat","w");
  neighbors_fib_file= fopen("neighbors_fib_file.dat", "w");
  stim_myo_file = fopen("stim_myo_file.dat","w");
  stim_fib_file = fopen("stim_fib_file.dat","w");

  //Define boost library random number generator
  typedef boost::mt19937 RNGType;
  static RNGType rng(static_cast<unsigned int>(std::time(0)));

  //Initialize myocyte counter
  int myo_count = 0;
  int fib_count = 0;


  double alpha = .1;//Ratio of fibroblasts to myocytes in final distribution
  double prob = alpha/(1 + alpha);
  int delta_sq=1;//Initial increment of the search
  int stim_number;
  int stim_row_max = NROWS/2;//stimulate half the rows
  const int stim_max = 2;
  int stim_fib_index = 0;
  int stim_myo_index = 0;

  const int stim_length = NROWS*NCOLS;//just set file to make value
  int stim_myo_array[stim_length]={0};//split into stim myo and stim fib array
  int stim_fib_array[stim_length]={0};//split into stim myo and stim fib array  
////////////////////////////////////////////////////////////////////////////////////////////
//Using the boost library to create a binomial distribution with success probability (prob)
  boost::binomial_distribution<> binomial(1, prob);
  boost::variate_generator < RNGType, boost::binomial_distribution<> >  trial(rng, binomial);
/////////////////////////////////////////////////////////////////////////////////////////////
//1. Constructs a matrix NROWS by NCOLS with random fibroblast and myocyte placement. A single myocyte
//occupies 5 grid points and a single fibroblast occupies 1 grid point.
//2. Creates a list of cells to be stimulated, potentially more complex stimulus patterns (move to
//main code)
  for (int r=0; r<NROWS; r++){
    stim_number = 0;
    for (int c=0; c<NCOLS; c+=delta_sq){
      int is_fibroblast  = trial();
      if (is_fibroblast == 1){ //if current cell is a fibroblast
        fib_count--;//decrement fibroblast counter for negative values
        grid[r][c] = fib_count;//assign to current cell

        if( (stim_number < stim_max) && (r < stim_row_max)){//stimulates a subset of cols and rows
          stim_fib_array[stim_fib_index] = 1;
          stim_number++;
        }
        else{ 
          stim_fib_array[stim_fib_index] = 0;
        }
        stim_fib_index++;
        delta_sq = length_fib;//increment step based on cell type
      }
      else{ //if current cell is a myocyte
        myo_count++;//increment myocyte counter for positive values
        //grid[r][c] = myo_count;//assign current cell as myocyte
        for (int d=0;d<5;d++){
          grid[r][c+d] = myo_count;//assign current cell and its 4 neighbors to myocyte
        }
        delta_sq = length_myo;
        if ((stim_number < stim_max)&& (r < stim_row_max)){ 
          stim_myo_array[stim_myo_index] = 1;
          stim_number++;
          stim_myo_index++;
        }else{
          stim_myo_array[stim_myo_index] = 0;
          stim_myo_index++;
        }

      }
    }
  }
///////////////////////////////////////////////////////////////////////////////////////////////
  //Check distribution
  cout <<"Set alpha: "<< alpha <<endl;
  cout <<"Calc. Probability: "<< prob <<endl;
  cout << "# Fibroblasts: " << -fib_count << endl;
  cout << "# Myocytes: " << myo_count << endl;
  cout << "Ratio fib/myo: " << double(-fib_count)/double(myo_count) << endl;
//////////////////////////////////////////////////////////////////////////////////////////////

  //Write tab delimited grid to file
  for (int cc = 0; cc <NROWS; cc++){
    for (int bb = 0; bb < NCOLS; bb++){
      fprintf(grid_file, "%d ", grid[cc][bb]);
    }
    fprintf(grid_file, "\n");
  }
fclose(grid_file);


    //Write stim_fib_array to file
  for (int bb = 0; bb < stim_myo_index; bb++){
    fprintf(stim_myo_file, "%d ", stim_myo_array[bb]);
  }
fclose(stim_myo_file);

  //Write stim_fib_array to file
  for (int bb = 0; bb < stim_fib_index; bb++){
    fprintf(stim_fib_file, "%d ", stim_fib_array[bb]);
  }
fclose(stim_fib_file);

  //Determine neighbors and store into vector V
  // vector < vector<int> > V;

///////////////////////////////////////////////////////////////////////////////////////////////
//Define neighbors if each cell (fibroblast (1 grid point), myocyte (5 grid points)
//Decision tree:
  for (int r=0; r < NROWS; r++){
    for (int c=0; c < NCOLS; c+=delta_sq){
      vector <int> neighbor;//Create empty vector to store neighbors

      neighbor.push_back(grid[r][c]);//First entry in row is cell value (skipped in simulation)

        if (grid[r][c] < 0){//If current cell is a fibroblast

        if ((c-1) < 0){//If col# - 1 is less than zero its at the left edge of grid
          neighbor.push_back(grid[r][c]);//define its boundary conditions by setting neighbor to self
        }
        else{//If col# is not at left edge then it has a neighbor to the left
          neighbor.push_back(grid[r][c-1]);//use grid to find the value of its left neighbor
        }

        if ( (r-1) < 0){//If row# - 1 is less than zero its at the top edge of grid
          neighbor.push_back(grid[r][c]);//define its boundary condition by setting neighbor to self
        }
        else{//If row# is not at top edge then it has a neighbor above it
          neighbor.push_back(grid[r-1][c]);//use grid to find the value of its above neighbor
        }

        if ((r+1) > NROWS-1){//If row# + 1 is less than zero its at the bottom edge of grid
          neighbor.push_back(grid[r][c]);//define its boundary condition by setting neighbor to self
        }
        else{//If row# is not at bottom edge then it has a neighbor above it
          neighbor.push_back(grid[r+1][c]);
        }

        if ((c+1) > NCOLS-1){//If col# + 1 is greater than zero its at the right edge of grid
          neighbor.push_back(grid[r][c]);//define its boundary condition by setting neighbor to self
        }
        else{ //If col# is not at right edge then it has a neighbor to the right
          neighbor.push_back(grid[r][c+1]);
        }

        delta_sq = length_fib;//increment fibroblast
      }
      else{//If current cell is a myocyte
        if ((c-1) < 0){
          neighbor.push_back(grid[r][c]);
        }
        else{
          neighbor.push_back(grid[r][c-1]);
        }
        if ( (r-1) < 0){
        for (int i = 0; i < length_myo; i++){
            if ((c+i) < NCOLS){
              neighbor.push_back(grid[r][c+i]);
            }
          }
        }
        else{
          for (int i = 0; i < length_myo; i++){
            if ((c+i) < NCOLS){
              neighbor.push_back(grid[r-1][c+i]);
            }
          }
        }
        if ((r+1) > NROWS-1){
          for (int i = 0; i < length_myo; i++){
            if ((c+i) < NCOLS){
              neighbor.push_back(grid[r][c+i]);
            }
          }
        }
        else{
          for (int i = 0; i < length_myo; i++){
            if ((c+i) < NCOLS){
              neighbor.push_back(grid[r + 1][c+i]);
            }
          }
        }
        if ((c+length_myo) > NCOLS-1){//Cell is partially truncated or at the edge
            neighbor.push_back(grid[r][c]);
        }
        else{
            neighbor.push_back(grid[r][c+length_myo]);
        }

        delta_sq = length_myo;//Increment the length of a myocyte
      }
//////////////////////////////////////////////////////////////////////////////////////
      if (neighbor[0] < 0){
        //Print fibroblast neighbors into space delimited file
        for (int i=0; i<neighbor.size();i++){
          fprintf(neighbors_fib_file, "%d ", neighbor[i] );
        }
        fprintf(neighbors_fib_file, "\n");//return
      }
      else{
        //Print myocyte neighbors into space delimited file
        for (int i=0; i<neighbor.size();i++){
          fprintf(neighbors_myo_file, "%d ", neighbor[i] );
        }
        fprintf(neighbors_myo_file, "\n");

      }
      //  V.push_back(neighbor);
    }
  }

  fclose(neighbors_fib_file);
  fclose(neighbors_myo_file);
  return 0;
}

