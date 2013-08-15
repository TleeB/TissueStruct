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
const int NROWS = 10, NCOLS = 10;
int grid[NROWS][NCOLS]={0};


int main(){//Forgot the double parenthesis
  FILE *grid_file;
  FILE *neighbors_myo_file;
  FILE *neighbors_fib_file;
 // FILE *stim_myo_file;
 // FILE *stim_fib_file;

  grid_file = fopen("grid_file.dat", "w");
  neighbors_myo_file = fopen("neighbors_myo_file.dat","w");
  neighbors_fib_file= fopen("neighbors_fib_file.dat", "w");
//  stim_myo_file = fopen("stim_myo_file.dat","w");
//  stim_fib_file = fopen("stim_fib_file.dat","w");

  //Define boost library random number generator
  typedef boost::mt19937 RNGType;
  static RNGType rng(static_cast<unsigned int>(std::time(0)));

  //Initialize myocyte counter
  int myo_count = 0;
  int fib_count = 0;


  double alpha = .1;//Ratio of fibroblasts to myocytes in final distribution
  double prob = alpha/(1 + alpha);
  int delta_sq=1;//Initial increment of the search
  int cell_length, cell_marker;//Holder for variable cell lengths and marker
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
  for (int r=0; r<NROWS; r ++){
    for (int c=0; c < NCOLS; c += delta_sq){
      int is_fibroblast  = trial();
      if (is_fibroblast == 1){ //if current cell is a fibroblast
        fib_count--;//decrement fibroblast counter for negative values
        cell_marker = fib_count;//set the cell marker to fibroblast
        cell_length = length_fib;
      }
      else{ //if current cell is a myocyte
        myo_count++;//increment myocyte counter for positive values
        cell_marker = myo_count;
        cell_length = length_myo;
      }
      for (int d=0; d < cell_length; d++){
        grid[r][c + d] = cell_marker;//assign current cell with current marker 
      }
      delta_sq = cell_length;//increment step based on cell type
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


/*
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
*/

  //Determine neighbors and store into vector V
  // vector < vector<int> > V;

  ///////////////////////////////////////////////////////////////////////////////////////////////
  //Define neighbors if each cell (fibroblast (usually 1 grid point), myocyte (usually 5 grid points)
  //Decision tree:
  for (int r=0; r < NROWS; r++){
    for (int c=0; c < NCOLS; c+=delta_sq){
      vector <int> neighbor;//Create empty vector to store neighbors

      neighbor.push_back(grid[r][c]);//First entry in row is cell value (skipped in simulation)

      if (grid[r][c] < 0){//If current cell is a fibroblast
        cell_length = length_fib;
      }
      else{//If current cell is a myocyte
        cell_length = length_myo;
      }
      if ((c-1) < 0){
        neighbor.push_back(grid[r][c]);
      }
      else{
        neighbor.push_back(grid[r][c-1]);
      }
      if ( (r-1) < 0){
        for (int i = 0; i < cell_length; i++){
          if ((c+i) < NCOLS){
            neighbor.push_back(grid[r][c+i]);
          }
        }
      }
      else{
        for (int i = 0; i < cell_length; i++){
          if ((c+i) < NCOLS){
            neighbor.push_back(grid[r-1][c+i]);
          }
        }
      }
      if ((r+1) > NROWS-1){
        for (int i = 0; i < cell_length; i++){
          if ((c+i) < NCOLS){
            neighbor.push_back(grid[r][c+i]);
          }
        }
      }
      else{
        for (int i = 0; i < cell_length; i++){
          if ((c+i) < NCOLS){
            neighbor.push_back(grid[r + 1][c+i]);
          }
        }
      }
      if ((c+length_myo) > NCOLS-1){//Cell is partially truncated or at the edge
        neighbor.push_back(grid[r][c]);
      }
      else{
        neighbor.push_back(grid[r][c+cell_length]);
      }
      //////////////////////////////////////////////////////////////////////////////////////
      //Write current cells neighbors to file
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
      //////////////////////////////////////////////////////////////////////////////////////      
      //  V.push_back(neighbor);
      delta_sq = cell_length;//Increment the length of a myocyte
    }
  }
  fclose(neighbors_fib_file);
  fclose(neighbors_myo_file);
  return 0;
}

