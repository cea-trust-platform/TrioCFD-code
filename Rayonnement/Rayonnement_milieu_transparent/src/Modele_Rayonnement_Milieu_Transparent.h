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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Modele_Rayonnement_Milieu_Transparent.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement/src
// Version:     /main/19
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modele_Rayonnement_Milieu_Transparent_included
#define Modele_Rayonnement_Milieu_Transparent_included

#define Sigma_DEFAULT 5.67e-8

#include <Vect_Face_Rayonnante.h>
#include <Modele_Rayonnement_base.h>
#include <Matrice_Morse.h>

class Schema_Temps_base;
class Discretisation_base;

class Modele_Rayonnement_Milieu_Transparent: public Modele_Rayonnement_base
{
  Declare_instanciable_sans_constructeur(Modele_Rayonnement_Milieu_Transparent);


public:

  void discretiser(const Discretisation_base&, const Domaine& ) override;
  void mettre_a_jour (double temps) override;
  void preparer_calcul () override;
  double flux_radiatif( int num_face_global) const; // 0 < face < nb_faces_de_bord
  inline int nb_faces_rayonnantes() const;
  inline int nb_faces_totales() const;
  inline int ordre_matrice_fac_forme() const;
  inline const Face_Rayonnante& face_rayonnante(int ) const;
  inline  Face_Rayonnante& face_rayonnante(int );
  inline Modele_Rayonnement_Milieu_Transparent(double Sigma = Sigma_DEFAULT);
  void lire_fichiers(Nom& nom1, Nom& nom2);
  void lire_fichiers(Nom& nom1, Nom& nom2, Nom& nom3);
  inline int processeur_rayonnant();
  inline int processeur_rayonnant() const;
  inline void associer_processeur_rayonnant(int proc);
  inline double relaxation()  const;
  inline const Nom&  nom_pb_rayonnant() const;
  inline  Nom&  nom_pb_rayonnant() ;
private :
  int nb_faces_rayonnantes_;
  int nb_faces_totales_;
  int ordre_mat_forme_;
  //IntTab temporaire;
  VECT(Face_Rayonnante) les_faces_rayonnantes;
  DoubleTab les_facteurs_de_forme;
  DoubleTab matrice_rayo;
  DoubleTab les_flux_radiatifs;
  mutable IntVect corres;
  // void calculer_temperatures(Cond_lim_base&);
  void calculer_temperatures() override;
  void calculer_radiosites() override;
  void calculer_flux_radiatifs() override;
  void imprimer_flux_radiatifs(Sortie& ) const override;
  double temps_; // on garde le temps pour les impressions
  mutable int deja_imprime_;
protected :
  double SIGMA;
  int inversion_debut_;
  int lire_matrice_inv_;
  int fic_mat_ray_inv_bin_;
  Nom nom_fic_mat_ray_inv_;
  int processeur_rayonnant_;
  double relaxation_;
  Nom nom_pb_rayonnant_;

};


inline void Modele_Rayonnement_Milieu_Transparent::associer_processeur_rayonnant(int proc)
{
  processeur_rayonnant_ = proc;
}

inline int Modele_Rayonnement_Milieu_Transparent::processeur_rayonnant()
{
  return processeur_rayonnant_;
}

inline int Modele_Rayonnement_Milieu_Transparent::processeur_rayonnant() const
{
  return processeur_rayonnant_;
}

inline Modele_Rayonnement_Milieu_Transparent::Modele_Rayonnement_Milieu_Transparent(double Sigma) : deja_imprime_(0),SIGMA(Sigma),relaxation_(1) {}

inline const Face_Rayonnante& Modele_Rayonnement_Milieu_Transparent::face_rayonnante(int j) const
{
  return les_faces_rayonnantes[j];
}

inline Face_Rayonnante& Modele_Rayonnement_Milieu_Transparent::face_rayonnante(int j)
{
  return les_faces_rayonnantes[j];
}

inline int Modele_Rayonnement_Milieu_Transparent::nb_faces_rayonnantes() const
{
  return nb_faces_rayonnantes_;
}

inline int Modele_Rayonnement_Milieu_Transparent::nb_faces_totales() const
{
  return nb_faces_totales_;
}

inline int Modele_Rayonnement_Milieu_Transparent::ordre_matrice_fac_forme() const
{
  return  ordre_mat_forme_;
}
inline double Modele_Rayonnement_Milieu_Transparent::relaxation() const
{
  return  relaxation_;
}
inline const Nom& Modele_Rayonnement_Milieu_Transparent::nom_pb_rayonnant() const
{
  return  nom_pb_rayonnant_;
}
inline  Nom& Modele_Rayonnement_Milieu_Transparent::nom_pb_rayonnant()
{
  return  nom_pb_rayonnant_;
}
#endif
