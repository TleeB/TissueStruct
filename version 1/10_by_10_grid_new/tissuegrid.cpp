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
int graph[NROWS][NCOLS]={0};


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
  const int stim_max = 1;
  int stim_fib_index = 0;
  int stim_myo_index = 0;

  const int stim_length = NROWS*NCOLS;//just set file to make value
  int stim_myo_array[stim_length]={0};//split into stim myo and stim fib array
  int stim_fib_array[stim_length]={0};//split into stim myo and stim fib array  

  boost::binomial_distribution<> binomial(1, prob);
  boost::variate_generator < RNGType, boost::binomial_distribution<> >  trial(rng, binomial);

  for (int r=0; r<NROWS; r++){
    stim_number = 0;

    for (int c=0; c<NCOLS; c+=delta_sq){
      int is_fibroblast  = trial();
      if (is_fibroblast == 1){ //if current cell is a fibroblast
        fib_count--;
        graph[r][c] = fib_count;
        delta_sq = length_fib;
        if( (stim_number < stim_max)&& (r < stim_row_max)){//stimulate half the rows
          stim_fib_array[stim_fib_index] = 1;
          stim_number++;
          stim_fib_index++;
        }else{
          stim_fib_array[stim_fib_index] = 0;
          stim_fib_index++;
        }

      }
      else{
        myo_count++;
        graph[r][c] = myo_count;
        for (int d=1;d<=4;d++){
          graph[r][c+d] = myo_count;
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
  //Check distribution
  cout <<"Alpha: "<< alpha <<endl;
  cout <<"Probability: "<< prob <<endl;
  cout << "Fibroblasts: " << -fib_count << endl;
  cout << "Myocytes: " << myo_count << endl;
  cout << "Ratio: " << double(-fib_count)/double(myo_count) << endl;

  //Also print stim_fib_array
  for (int bb = 0; bb < stim_myo_index; bb++){
    fprintf(stim_myo_file, "%d ", stim_myo_array[bb]);
  }
  fclose(stim_myo_file);

  //Also print stim_fib_array
  for (int bb = 0; bb < stim_fib_index; bb++){
    fprintf(stim_fib_file, "%d ", stim_fib_array[bb]);
  }
  fclose(stim_fib_file);

  //Write tab delimited grid to file
  for (int cc = 0; cc <NROWS; cc++){
    for (int bb = 0; bb < NCOLS; bb++){
      fprintf(grid_file, "%d ", graph[cc][bb]);
    }
    fprintf(grid_file, "\n");
  }
  fclose(grid_file);


  //Determine neighbors and store into vector V
  // vector < vector<int> > V;
  int col_index, prev_cell, curr_cell;
  col_index = 0;


  for (int r=0; r<NROWS; r++){
    //delta_sq = 1;
    for (int c=0; c<NCOLS; c+=delta_sq){
      vector <int> neighbor;
      neighbor.push_back(graph[r][c]);
      if (graph[r][c] < 0){
        if ((c-1) > 0){
          neighbor.push_back(graph[r][c - 1]);
        }
        else{
          neighbor.push_back(graph[r][c]);
        }
        if ((r+1)< NROWS){ //double check that I am not going out of bounds
          neighbor.push_back(graph[r + 1][c]);
        }
        else{
          neighbor.push_back(graph[r][c]);
        }
        if ( (r-1) >= 0){
          neighbor.push_back(graph[r - 1][c]);
        }
        else{
          neighbor.push_back(graph[r][c]);
        }
        if ((c+1)< NCOLS){
          neighbor.push_back(graph[r][c + 1]);
        }
        else{
          neighbor.push_back(graph[r][c]);
        }
        delta_sq = length_fib;
      }
      else{
        if ((c-1) > 0){
          neighbor.push_back(graph[r][c - 1]);
        }
        else{
          neighbor.push_back(graph[r][c]);
        }
        for (int i = 0; i < length_myo; i++){
          if ((r+1)< NROWS){ //double check that I am not going out of bounds
            neighbor.push_back(graph[r + 1][c+i]);
          }
          else{
            neighbor.push_back(graph[r][c+i]);
          }
        }
        for (int i = 0; i < length_myo; i++){
          if ( (r-1) >= 0){
            neighbor.push_back(graph[r - 1][c+i]);
          }
          else{
            neighbor.push_back(graph[r][c+i]);
          }
        }
        if ((c+length_myo)< NCOLS){
          neighbor.push_back(graph[r][c + length_myo]);
        }
        else{
          neighbor.push_back(graph[r][c]);
        }

        delta_sq = length_myo;
      }

      if (neighbor[0] < 0){
        //Print neighbors into comma delimited file
        for (int i=0; i<neighbor.size();i++){
          fprintf(neighbors_fib_file, "%d ", neighbor[i] );
        }
        fprintf(neighbors_fib_file, "\n");
      }
      else{
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

