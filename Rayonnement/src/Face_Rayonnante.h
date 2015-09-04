/****************************************************************************
* Copyright (c) 2015, CEA
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
// File:        Face_Rayonnante.h
// Directory:   $TRUST_ROOT/src/Rayonnement
// Version:     /main/14
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Face_Rayonnante_included
#define Face_Rayonnante_included

#include <VECT_Ensemble_Faces.h>

class Face_Rayonnante : public Objet_U
{
  Declare_instanciable(Face_Rayonnante);

public :

  inline double T_face_rayo() const ;
  inline double emissivite() const;
  inline double flux_radiatif() const;
  inline double surface_rayo() const;
  inline const Nom& nom_bord_rayo() const;
  inline const Nom& nom_bord_rayo_lu() const;

  inline const IntVect& num_faces(int j) const;
  inline const Ensemble_Faces_base& ensembles_faces_bord(int j) const;
  inline Ensemble_Faces_base& ensembles_faces_bord(int j);
  double calculer_temperature();
  //double calculer_flux_radiatif();
  double imprimer_flux_radiatif(Sortie&,Sortie&, Sortie&  ) const;
  int chercher_ensemble_faces(const Nom& ) const;
  inline int nb_ensembles_faces() const;
  inline void mettre_a_jour_flux_radiatif(double J);
  void ecrire_temperature_bord() const;
private :

  VECT(Ensemble_Faces_base) Les_ensembles_faces_bord;
  IntVect num_faces_;
  double T_face_rayo_;
  double emissivite_;
  double surf_;
  Nom nom_bord_rayo_,nom_bord_rayo_lu_;
  //double radiosite;
  double flux_radiatif_;
  int nb_ensembles_faces_;
};

inline const IntVect& Face_Rayonnante::num_faces(int j) const
{
  return num_faces_;
}

inline double Face_Rayonnante::T_face_rayo() const
{
  return T_face_rayo_;
}

inline double Face_Rayonnante::emissivite() const
{
  return emissivite_;
}

inline int Face_Rayonnante::nb_ensembles_faces() const
{
  return nb_ensembles_faces_;
}
inline double Face_Rayonnante::flux_radiatif() const
{
  return flux_radiatif_;
}

inline double Face_Rayonnante::surface_rayo() const
{
  return  surf_;
}

inline const Nom& Face_Rayonnante::nom_bord_rayo() const
{
  return nom_bord_rayo_;
}
inline const Nom& Face_Rayonnante::nom_bord_rayo_lu() const
{
  return nom_bord_rayo_lu_;
}

inline void Face_Rayonnante::mettre_a_jour_flux_radiatif(double J)
{
  flux_radiatif_ = J;
}

inline const Ensemble_Faces_base& Face_Rayonnante::ensembles_faces_bord(int j) const
{
  return Les_ensembles_faces_bord[j];
}

inline Ensemble_Faces_base& Face_Rayonnante::ensembles_faces_bord(int j)
{
  return Les_ensembles_faces_bord[j];
}
#endif
