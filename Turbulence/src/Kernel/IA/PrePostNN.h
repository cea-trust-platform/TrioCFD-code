// TRUST_NO_INDENT
#include <string>
#include <vector>

using namespace std;

enum pp_lambda {INDEFL,LNORM,LU,LUS};
enum pp_T {INDEFT,TF,TR};

class PrePostNN
{
public:
  PrePostNN(string filename);
  ~PrePostNN();

  void AllDisplay();

  vector<double> get_alpha() {return alpha;}
  vector<double> get_lmean() {return lmean;}
  vector<double> get_lmax() {return lmax;}
  vector<double> get_tfn() {return tfn;}
  double get_bsigma() {return bsigma;}
  double get_t_thresh() {return t_thresh;}
  vector<vector<double>> get_lambda_au() {return lambda_au;}
  vector<vector<double>> get_lambda_as() {return lambda_as;}
  enum pp_lambda get_ppl() {return ppl;}
  enum pp_T get_ppt() {return ppt;}
  vector<int> get_ilambda() {return ilambda;}
  vector<int> get_iT() {return iT;}

private:
  void display(string tag,vector<double> vec);
  void display(string tag,vector<int> vec);
  void display(string tag,vector<vector<double>> mat);
  string trim(const std::string& str, const std::string& whitespace = " \t");
  vector<double> ReadDataFromLine(string buffer,string tag,size_t npos);
  double ReadOneDataFromLine(string buffer,string tag,size_t npos);
  enum pp_T ReadPPTFromLine(string buffer,string tag,size_t npos);
  enum pp_lambda ReadPPLFromLine(string buffer,string tag,size_t npos);
  vector<int> ReadIndexFromLine(string buffer,string tag,size_t npos);
  vector<vector<double>> ReadDataFromSeveralLines(ifstream &f,int nblines);

  vector<double> alpha;
  vector<double> lmean;
  vector<double> lmax;
  vector<double> tfn;
  double bsigma;
  double t_thresh;
  vector<vector<double>> lambda_au;
  vector<vector<double>> lambda_as;
  enum pp_lambda ppl;
  enum pp_T ppt;
  vector<int> ilambda;
  vector<int> iT;
  
};
