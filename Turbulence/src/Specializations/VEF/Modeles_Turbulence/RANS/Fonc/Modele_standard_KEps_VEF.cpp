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
// File:        Modele_standard_KEps_VEF.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Modeles_Turbulence/RANS/Fonc
//
//////////////////////////////////////////////////////////////////////////////
#include <Modele_standard_KEps_VEF.h>

#include <Domaine_VEF.h>
#include <Domaine_Cl_VEF.h>
#include <Champ_Uniforme.h>
#include <Scatter.h>
#include <Champ_P1NC.h>
#include <Champ_P0_VEF.h>
#include <Discretisation_base.h>
#include <Check_espace_virtuel.h>
#include <LecFicDiffuse.h>
#include <EcritureLectureSpecial.h>

Implemente_instanciable(Modele_standard_KEps_VEF,"Modele_standard_KEps_VEF",Modele_Lam_Bremhorst_VEF);

// XD standard_KEps Lam_Bremhorst standard_KEps -1 Model described in ' E. Baglietto , CFD and DNS methodologies development for fuel bundle simulaions, Nuclear Engineering and Design, 1503--1510 (236), 2006. '

// printOn et readOn

Sortie& Modele_standard_KEps_VEF::printOn(Sortie& s ) const
{
  return s;
}

Entree& Modele_standard_KEps_VEF::readOn(Entree& is )
{
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades_depuis(is);
  if (is_Reynolds_stress_isotrope_ == 0)
    is_Cmu_constant_ = 0;
  return is;
}

void Modele_standard_KEps_VEF::lire_distance_paroi( )
{

// Pas besoin de la distance a la proi dans ce modele

}

DoubleTab& Modele_standard_KEps_VEF::Calcul_D(DoubleTab& D,const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,
                                              const DoubleTab& vitesse,const DoubleTab& K_eps_Bas_Re, const Champ_Don& ch_visco ) const
{
  D = 0;
  return D;
}

DoubleTab& Modele_standard_KEps_VEF::Calcul_E(DoubleTab& E,const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis, const DoubleTab& transporte,const DoubleTab& K_eps_Bas_Re,const Champ_Don& ch_visco, const DoubleTab& visco_turb ) const
{
  E = 0;
  return E;
}

/* DoubleTab& Modele_standard_KEps_VEF::Calcul_F1( DoubleTab& F1, const Domaine_dis& domaine_dis) const
{
  const Domaine_VEF& le_dom = ref_cast(Domaine_VEF,domaine_dis.valeur());
  int nb_faces = le_dom.nb_faces();
  for (int num_face=0; num_face <nb_faces; num_face ++ )
    F1[num_face] = 1.;
  return F1;
}
*/
DoubleTab& Modele_standard_KEps_VEF::Calcul_F1( DoubleTab& F1, const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis, const DoubleTab& P, const DoubleTab& K_eps_Bas_Re,const Champ_base& ch_visco) const
{
  const Domaine_VEF& le_dom = ref_cast(Domaine_VEF,domaine_dis.valeur());
  int nb_faces = le_dom.nb_faces();
  int num_face;
  // Calcul de F1
  for (num_face=0; num_face <nb_faces; num_face ++ )
    {
      F1[num_face] = 1.;
    }
  return F1;
}

DoubleTab& Modele_standard_KEps_VEF::Calcul_F2( DoubleTab& F2, DoubleTab& Deb, const Domaine_dis& domaine_dis,const DoubleTab& K_eps_Bas_Re,const Champ_base& ch_visco ) const
{
  const Domaine_VEF& le_dom = ref_cast(Domaine_VEF,domaine_dis.valeur());
  int nb_faces = le_dom.nb_faces();
  int num_face;

  for (num_face=0; num_face<nb_faces  ; num_face++)
    {
      F2[num_face] = 1.;
    }
  //Cerr<<F2.mp_min_vect()<<" F2 "<<F2.mp_max_vect()<<finl;
  return F2;
}

DoubleTab&  Modele_standard_KEps_VEF::Calcul_Fmu( DoubleTab& Fmu,const Domaine_dis& domaine_dis,const Domaine_Cl_dis& domaine_Cl_dis,const DoubleTab& K_eps_Bas_Re,const Champ_Don& ch_visco ) const
{
  const Domaine_VEF& le_dom = ref_cast(Domaine_VEF,domaine_dis.valeur());
  int nb_faces = le_dom.nb_faces();
  int num_face;

  for (num_face=0; num_face <nb_faces; num_face ++ )
    {
      Fmu(num_face) = 1.;
    }

  return Fmu;
}

DoubleTab&  Modele_standard_KEps_VEF::Calcul_Fmu_BiK( DoubleTab& Fmu,const Domaine_dis& domaine_dis,const Domaine_Cl_dis& domaine_Cl_dis,const DoubleTab& K_Bas_Re,const DoubleTab& eps_Bas_Re,const Champ_Don& ch_visco ) const
{
  const Domaine_VEF& le_dom = ref_cast(Domaine_VEF,domaine_dis.valeur());
  int nb_faces = le_dom.nb_faces();
  int num_face;

  for (num_face=0; num_face <nb_faces; num_face ++ )
    {
      Fmu(num_face) = 1.;
    }

  return Fmu;
}

DoubleTab& Modele_standard_KEps_VEF::Calcul_F2_BiK( DoubleTab& F2, DoubleTab& Deb, const Domaine_dis& domaine_dis,const DoubleTab& K_Bas_Re,const DoubleTab& eps_Bas_Re,const Champ_base& ch_visco ) const
{
  const Domaine_VEF& le_dom = ref_cast(Domaine_VEF,domaine_dis.valeur());
  int nb_faces = le_dom.nb_faces();
  int num_face;

  for (num_face=0; num_face<nb_faces  ; num_face++)
    {
      F2[num_face] = 1.;
    }
  //Cerr<<F2.mp_min_vect()<<" F2 "<<F2.mp_max_vect()<<finl;
  return F2;
}

DoubleTab& Modele_standard_KEps_VEF::Calcul_F1_BiK( DoubleTab& F1, const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis, const DoubleTab& P, const DoubleTab& K_Bas_Re, const DoubleTab& eps_Bas_Re,const Champ_base& ch_visco) const
{
  const Domaine_VEF& le_dom = ref_cast(Domaine_VEF,domaine_dis.valeur());
  int nb_faces = le_dom.nb_faces();
  int num_face;
  // Calcul de F1
  for (num_face=0; num_face <nb_faces; num_face ++ )
    {
      F1[num_face] = 1.;
    }
  return F1;
}

DoubleTab& Modele_standard_KEps_VEF::Calcul_E_BiK(DoubleTab& E,const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis, const DoubleTab& transporte,const DoubleTab& K_Bas_Re,const DoubleTab& eps_Bas_Re,const Champ_Don& ch_visco, const DoubleTab& visco_turb ) const
{
  E = 0;
  return E;
}

DoubleTab& Modele_standard_KEps_VEF::Calcul_D_BiK(DoubleTab& D,const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,
                                                  const DoubleTab& vitesse,const DoubleTab& K_Bas_Re,const DoubleTab& eps_Bas_Re, const Champ_Don& ch_visco ) const
{
  D = 0;
  return D;
}






