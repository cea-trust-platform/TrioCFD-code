#include <keras_model.h>
#include <PrePostNN.h>

class TBNN
{
public:
  TBNN(string keras_model_file,string preproc_file);
  ~TBNN();

  vector<double> predict(vector<double> lambda, vector<vector<double>> T);

private:

// Neural network
  KerasModel _model;              // Neural network
  vector<double> _g;              // output of the Neural network
  PrePostNN *_ppNN;               // objet pre et post processing
  
// Lambda pre-processing
  vector<double> _plambda;        // pre-processed lambda

// T pre-processing
  vector<vector<double>> _pT;     // pre-processed T

// b post-processing
  vector<double> _pb;             // result of the prediction not post-processed
  vector<double> _b;              // result of the prediction post-processed

// process methods
  void process_lambda(vector<double> lambda); // calculate the pre-processed lambda
  void process_T(vector<vector<double>> T);   // calculate the pre-processed T
  void applyNN();                             // prédiction du réseau de neurones
  void process_b();                           // calculate the post-processed b
};
