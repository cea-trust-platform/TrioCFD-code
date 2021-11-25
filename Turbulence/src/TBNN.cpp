#include <TBNN.h>
#include <iostream>

int sgn(double v) {
  return (v > 0) - (v < 0);
}

TBNN::TBNN(string keras_model_file,string preproc_file)
{
  _ppNN = new PrePostNN(preproc_file);
  _model.LoadModel(keras_model_file);
}

TBNN::~TBNN()
{
  delete(_ppNN);
}

vector<double> TBNN::predict(vector<double> lambda, vector<vector<double>> T)
{
  process_lambda(lambda);
  process_T(T);
  applyNN();
  process_b();
  return(_b);
}

void TBNN::process_lambda(vector<double> lambda)
{
  vector<double> lc;                 // lambda centre
  vector<double> lcr;                // lambda centre reduit
  unsigned int nbl = lambda.size();  // nombre d'invariants lambda

  lc.resize(nbl);
  lcr.resize(nbl);
  _plambda.resize(nbl);

  switch(_ppNN->get_ppl())
  {
  case LNORM:
    // centrage des lambda
    if( _ppNN->get_lmean().size() == nbl )
      for(unsigned int i=0;i<nbl;i++)
        lc[i] = lambda[i] - _ppNN->get_lmean()[i];
    else
      for(unsigned int i=0;i<nbl;i++)
        lc[i] = lambda[i];
    // reduction des lambda
    for(unsigned int i=0;i<nbl;i++)
      _plambda[i] = lc[i] / _ppNN->get_lmax()[i];
    break;
  case LU:
    // puissance alpha
    if( _ppNN->get_alpha().size() == nbl )
      for(unsigned int i=0;i<nbl;i++)
        lc[i] = sgn(lambda[i]) * pow(std::fabs(lambda[i]),_ppNN->get_alpha()[i]);
    else
      for(unsigned int i=0;i<nbl;i++)
        lc[i] = lambda[i];
    // centrage des lambda
    if( _ppNN->get_lmean().size() == nbl )
      for(unsigned int i=0;i<nbl;i++)
        lc[i] = lc[i] - _ppNN->get_lmean()[i];
    // reduction des lambda
    for(unsigned int i=0;i<nbl;i++)
      lcr[i] = lc[i] / _ppNN->get_lmax()[i];
    // multiplication par transpose(au)
    for(unsigned int i=0;i<nbl;i++){
      _plambda[i] = 0.;
      for(unsigned int j=0;j<nbl;j++)
	_plambda[i] += _ppNN->get_lambda_au()[j][i] * lcr[j];
    }
    break;
  case LUS:
    // puissance alpha
    if( _ppNN->get_alpha().size() == nbl )
      for(unsigned int i=0;i<nbl;i++)
        lc[i] = sgn(lambda[i]) * pow(std::fabs(lambda[i]),_ppNN->get_alpha()[i]);
    else
      for(unsigned int i=0;i<nbl;i++)
        lc[i] = lambda[i];
    // centrage des lambda
    if( _ppNN->get_lmean().size() == nbl )
      for(unsigned int i=0;i<nbl;i++)
        lc[i] = lc[i] - _ppNN->get_lmean()[i];
    // reduction des lambda
    for(unsigned int i=0;i<nbl;i++)
      lcr[i] = lc[i] / _ppNN->get_lmax()[i];
    // multiplication par transpose(as)
    for(unsigned int i=0;i<nbl;i++){
      _plambda[i] = 0.;
      for(unsigned int j=0;j<nbl;j++)
	_plambda[i] += _ppNN->get_lambda_as()[j][i] * lcr[j];
    }
    break;
  default:
    // on lance une exception
    cerr << "Mauvaise methode de pre traitement des lambda" << endl;
    break;
  }
}

void TBNN::process_T(vector<vector<double>> T)
{
  unsigned int nbt = T.size();    // nombre de tenseurs T
  unsigned int nbb = T[0].size(); // taille de chacun des tenseurs T

  _pT.resize(nbt);
  for(unsigned i=0;i<nbt;i++)
    _pT[i].resize(nbb);

  // le premier tenseur T0 reste inchange
  for(unsigned int j=0;j<nbb;j++)
    _pT[0][j] = T[0][j];

  // pre process de T
  switch(_ppNN->get_ppt())
  {
  case TF:
    // on calcule la norme de Frobenius de chaque tenseur
    for(unsigned int i=1;i<nbt;i++){
      double normf = 0.;
      for(unsigned int j=0;j<nbb;j++)
	normf += T[i][j] * T[i][j];
      normf = sqrt(normf);
      for(unsigned int j=0;j<nbb;j++)
	_pT[i][j] = T[i][j] / (normf + _ppNN->get_t_thresh());
    }
    break;
  case TR:
    // on divise le tenseur Ti par la norme globale
    for(unsigned int i=1;i<nbt;i++)
      for(unsigned int j=0;j<nbb;j++)
	_pT[i][j] = T[i][j] / _ppNN->get_tfn()[i-1];
    break;
  default:
    // on lance une exception
    cerr << "Mauvaise methode de pre traitement des tenseurs T" << endl;
    break;
  }
}

void TBNN::process_b()
{
  vector<int> iT = _ppNN->get_iT();
  unsigned int nbt = iT.size();
  unsigned int nbb = _pT[0].size();

  // calcul de _pb a partir de _g et de _pT
  _pb.resize(nbb);
  for(unsigned int i=0;i<nbb;i++){
    _pb[i] = 0;
    for(unsigned int j=0;j<nbt;j++)
      _pb[i] += _g[j] * _pT[iT[j]][i];
  }

  // post process de b
  _b.resize(nbb);
  for(unsigned int i=0;i<nbb;i++)
    _b[i] = _ppNN->get_bsigma() * _pb[i];
}

void TBNN::applyNN()
{
  vector<int> il = _ppNN->get_ilambda();
  vector<int> iT = _ppNN->get_iT();
  Tensor in(il.size());
  Tensor out;
  unsigned int nbt = iT.size();

  _g.resize(nbt);

  // l'entree du reseau est determinee par les indices des lambdas definis dans le vecteur il
  for(unsigned int i=0;i<il.size();i++)
    in.data_[i] = (float)_plambda[il[i]];

  // on fait la prediction a l'aide du reseau de neurones
  _model.Apply(&in,&out);

  // on stocke les sorties dans _g
  for(unsigned int i=0;i<nbt;i++) _g[i] = out(i);

}
