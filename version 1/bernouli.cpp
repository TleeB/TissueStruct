#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
//must add flags -I /opt/local/include
//using the Boost library
#include <boost/random.hpp>
#include <boost/format.hpp>
using namespace std;

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

vector<vector<int> > fileToVector(const char *name)
{
    vector<vector<int> > result;
    ifstream input (name);
    string lineData;

    while(getline(input, lineData))
    {
        double d;
        vector<int> row;
        stringstream lineStream(lineData);

        while (lineStream >> d)
            row.push_back(d);

        result.push_back(row);
    }

    return result;
}
void PDE( double *V[], vector<vector<int> > & neighbor ){
double G;
    for (int p = 0; p < neighbor.size(); p++){
      for (int N = 0; N < neighbor[p].size(); N++){
        if(N > 0){ //If neighbor is a myocyte
          G = 200.;  //If the neighbor is end to end then conductance should be 600
        }
        else{ //If neighbor is a fibroblast
          G = 3.;
        }
        V[p] = G * (V[p] - V[ neighbor[p][N] ] );

      }
    }
  }
int main(){
//vector < vector<int> > V;
const char *filename = "neighbors_fib_file.dat";
vector < vector<int> > V = fileToVector(filename);

for (int i =0; i<V.size(); i++){
  for (int j =0; j<V[i].size(); j++){

cout << V[i][j] << " ";
  }
  cout<<endl;
}
  /*
  typedef boost::mt19937 RNGType;
  static RNGType rng(static_cast<unsigned int>(std::time(0)));
  double alpha = .5;
  boost::binomial_distribution<> binomial(1, 0.09);
  boost::variate_generator < RNGType, boost::binomial_distribution<> >  trial(rng, binomial);

for (int i=0; i<50; i++){
cout << trial() <<"\t";
}
cout << endl;
*/
  return 0;
}
