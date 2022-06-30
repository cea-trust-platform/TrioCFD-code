/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/
/////////////////////////////////////////////////////////////////////////////
//
// File      : Fourier_trans.h
// Directory : $NEW_ALGO_QC_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////
#ifndef Fourier_trans_H
#define Fourier_trans_H
#include <FixedVector.h>
#include <IJK_Field.h>
#include <Objet_U.h>
#include <TRUSTArrays.h>
#include <Redistribute_Field.h>
#include <fftw3.h>
#include <TRUSTTab.h>

// NBR DE TERMES INTER  !!!
// #define NVAL_INTER 32

// NBR DE TERMES DE SORTIE  !!!
// #define N_RES 21

class IJK_Grid_Geometry;
class Fourier_trans : public Objet_U
{
  Declare_instanciable(Fourier_trans);
public:
  void sauvegarde() const;
  void reprise();
  void postraiter(Nom&) const;
  void update(const FixedVector<IJK_Field_double, 3>& vitesse,
              const IJK_Field_double& pression,
              //const IJK_Field_double &temperature,
              const IJK_Field_double& masse_vol,
              const IJK_Field_double& champ_mu);

  void initialize( const IJK_Splitting&  /* splitting d'origine */
                   , const IJK_Splitting&  /* splitting de post-traitement*/
                   , VECT(ArrOfDouble) /* tableau des v mo */
                   , ArrOfDouble rhomoy
                   , ArrOfDouble numoy
                   , const double T_KMAX
                   , const double T_KMIN
                   , int avec_reprise );

  int N_echantillons() const
  {
    return N_echantillons_;
  }

protected:

  void Traitement_spectral(const IJK_Field_double& entree, IJK_Field_double& reel, IJK_Field_double& imag, int sens); // 1 pour direct -1 pour inverse
  void Traitement_spectral_klayer_direct(const IJK_Field_double& entree, IJK_Field_double& reel, IJK_Field_double& imag, int k);
  void Traitement_spectral_klayer_inverse(IJK_Field_double& reel, IJK_Field_double& imag, int k);
  void Ecrire_entete(Sortie& os, int flag_valeur_instantanee = 0) const;
  void postraiter_TFi(Sortie&, int flag_valeur_instantanee = 0) const;
  void postraiter_TFx(Sortie&, int flag_valeur_instantanee = 0) const;
  void postraiter_TFy(Sortie&, int flag_valeur_instantanee = 0) const;
  FixedVector<Redistribute_Field, 3> redistribute_to_post_splitting_faces_;
  Redistribute_Field redistribute_to_post_splitting_elem_;

  /* ici penser a fixer le resultat en fonction du nombre de param que l'on veux changer */
  FixedVector<IJK_Field_double, 27> resultat_ ;
  IJK_Splitting Post_splitting_;
  VECT(ArrOfDouble) vit_moy_;
  ArrOfDouble rho_moy_;
  ArrOfDouble nu_moy_;
  double TCL_kmax_;
  double TCL_kmin_;

  fftw_plan plan;
  fftw_plan inverseplan;
  fftw_complex* in;
  fftw_complex* out;

  // Z coordinates of statistics points
  ArrOfDouble elem_coord_;
  //taille des mailles;
  double dx_;
  double dy_;
  ArrOfDouble tab_dz_;
  // Last instantaneous value of the space average (only on processor 0)
  VECT(ArrOfDouble) moyenne_spatiale_instantanee_;
  // Temporal integral of statistics variables
  VECT(ArrOfDouble) integrale_temporelle_;
  // Integration number
  int N_echantillons_;

  /* ATTENTION Les tableaux ont ete deplacer ici pour pouvoir etre initialiser avant le calcl car cela coute aussi cher que la TF */
  FixedVector <IJK_Field_double, 86> Avant_TF ;
  FixedVector <IJK_Field_double, 86> Reel_TF ;
  FixedVector <IJK_Field_double, 86> Imag_TF ;
  FixedVector<IJK_Field_double, 3> vitesse ;
  IJK_Field_double champ_mu ;
  IJK_Field_double masse_vol ;
  IJK_Field_double pression ;
  /*
   IJK_Field_double Div_U;
   IJK_Field_double DUiDx;
   IJK_Field_double DUjDy;
   IJK_Field_double DUkDz;
  */
  VECT(ArrOfDouble) TFx_ ;
  VECT(ArrOfDouble) TFy_ ;
  VECT(ArrOfDouble) TFi_ ;
  ArrOfDouble tri_ki;
  IntTab ref_ki;
};
#endif
