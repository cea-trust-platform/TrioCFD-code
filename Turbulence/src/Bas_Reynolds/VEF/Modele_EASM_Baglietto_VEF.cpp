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
// File:        Modele_EASM_Baglietto_VEF.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Bas_Reynolds/src/VEF
// Version:     /main/20
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_EASM_Baglietto_VEF.h>
#include <Zone_VEF.h>
#include <Zone_Cl_VEF.h>
#include <Periodique.h>
#include <Champ_Uniforme.h>
#include <Scatter.h>
#include <Champ_P1NC.h>
#include <Champ_P0_VEF.h>
#include <Discretisation_base.h>
#include <Check_espace_virtuel.h>
#include <LecFicDiffuse.h>
#include <EcritureLectureSpecial.h>

Implemente_instanciable(Modele_EASM_Baglietto_VEF,"Modele_EASM_Baglietto_VEF",Modele_Lam_Bremhorst_VEF);
// XD EASM_Baglietto Lam_Bremhorst EASM_Baglietto -1 Model described in ' E. Baglietto and H. Ninokata ,  A turbulence model study for simulating flow inside tight lattice rod bundles, Nuclear Engineering and Design, 773--784 (235), 2005. '

// printOn et readOn

Sortie& Modele_EASM_Baglietto_VEF::printOn(Sortie& s ) const
{
  return s;
}

Entree& Modele_EASM_Baglietto_VEF::readOn(Entree& is )
{
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades_depuis(is);
  if (is_Reynolds_stress_isotrope_ == 0)
    is_Cmu_constant_ = 0;
  return is;
}


DoubleTab& Modele_EASM_Baglietto_VEF::Calcul_D(DoubleTab& D,const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,
                                               const DoubleTab& vitesse,const DoubleTab& K_eps_Bas_Re, const Champ_Don& ch_visco ) const
{
  D=0;
  return D;
}


DoubleTab& Modele_EASM_Baglietto_VEF::Calcul_E(DoubleTab& E,const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis, const DoubleTab& transporte,const DoubleTab& K_eps_Bas_Re,const Champ_Don& ch_visco, const DoubleTab& visco_turb ) const
{
  double kkk = CNL1;
  Cerr<<"kkk = " << kkk <<finl;
  E=0;
  return E;
}

/*DoubleTab& Modele_EASM_Baglietto_VEF::Calcul_F1( DoubleTab& F1, const Zone_dis& zone_dis) const
{
  const Zone_VEF& la_zone = ref_cast(Zone_VEF,zone_dis.valeur());
  int nb_faces = la_zone.nb_faces();
  for (int num_face=0; num_face <nb_faces; num_face ++ )
    F1[num_face] = 1.;
  return F1;
}
*/
DoubleTab& Modele_EASM_Baglietto_VEF::Calcul_F1( DoubleTab& F1, const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis, const DoubleTab& P, const DoubleTab& K_eps_Bas_Re,const Champ_base& ch_visco) const
{

  double visco=-1;
  const DoubleTab& tab_visco=ch_visco.valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco);
  if (is_visco_const)
    visco=tab_visco(0,0);
  const Zone_VEF& la_zone = ref_cast(Zone_VEF,zone_dis.valeur());
  const Zone_Cl_VEF& zone_Cl_VEF = ref_cast(Zone_Cl_VEF,zone_Cl_dis.valeur());
  const DoubleTab& wall_length = BR_wall_length_.valeurs();
  DoubleTab wall_length_face(0);
  la_zone.creer_tableau_faces(wall_length_face);
  DoubleTab Pderive(0);
  la_zone.creer_tableau_faces(Pderive);
  int nb_faces = la_zone.nb_faces();
  const Conds_lim& les_cl = zone_Cl_VEF.les_conditions_limites();
  int nb_cl=les_cl.size();
  const IntTab& face_voisins = la_zone.face_voisins();
  int num_face;
  double Rey,Re;
  /*
  	for (num_face=0; num_face <nb_faces; num_face ++ )
  		F1[num_face] = 1.;
  	return F1;
  */
  // Calcul de la distance a la paroi aux faces
  for (int n_bord=0; n_bord<nb_cl; n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
      int ndeb = le_bord.num_premiere_face();
      int nfin = ndeb + le_bord.nb_faces();

      if (sub_type(Periodique,la_cl.valeur()))
        {
          for (num_face=ndeb; num_face<nfin; num_face++)
            {
              int elem1 = face_voisins(num_face,0);
              int elem2 = face_voisins(num_face,1);
              wall_length_face(num_face) = 0.5*wall_length(elem1) + 0.5*wall_length(elem2);
            }
        }
      /* else if (sub_type(Pa,la_cl.valeur()))
      {
        for (int num_face=ndeb; num_face<nfin; num_face++)
        {
      	  int elem1 = face_voisins(num_face,0);
      	  int elem2 = face_voisins(num_face,1);
      			  wall_length_face(num_face) = 0.;
        }
      }*/
      else
        {
          for (num_face=ndeb; num_face<nfin; num_face++)
            {
              int elem1 = face_voisins(num_face,0);
              wall_length_face(num_face) = wall_length(elem1);
            }
        }
    }
  int n0 = la_zone.premiere_face_int();
  for (num_face=n0; num_face<nb_faces; num_face++)
    {
      int elem0 = la_zone.face_voisins(num_face,0);
      int elem1 = la_zone.face_voisins(num_face,1);
      wall_length_face(num_face) = 0.5*wall_length(elem0)+0.5*wall_length(elem1);
    }
  // Calcul de la distance a la paroi aux faces
  /*    for (num_face=0; num_face< la_zone.premiere_face_int(); num_face++)
      {
    	  int elem0 = la_zone.face_voisins(num_face,0);
    	  if (elem0 != -1)
    		  wall_length_face(num_face) = wall_length(elem0);
    	  else
    	  {
    		  elem0 = la_zone.face_voisins(num_face,1);
    		  wall_length_face(num_face) = wall_length(elem0);
    	  }
      }

      for (; num_face<nb_faces; num_face++)
      {
    	  int elem0 = la_zone.face_voisins(num_face,0);
    	  int elem1 = la_zone.face_voisins(num_face,1);
    	  wall_length_face(num_face) = 0.5*wall_length(elem0)+0.5*wall_length(elem1);
      }
  */
  for (num_face=0; num_face <nb_faces; num_face ++ )
    {
      Re = (K_eps_Bas_Re(num_face,0)*K_eps_Bas_Re(num_face,0))/(visco*K_eps_Bas_Re(num_face,1)+DMINFLOAT);
      Rey = wall_length_face(num_face)*sqrt(K_eps_Bas_Re(num_face,0))/visco;
      Pderive(num_face) = 1.33*(1. - 0.3*exp(-Re*Re))
                          * (P(num_face) + 2.*visco*K_eps_Bas_Re(num_face,0)/(wall_length_face(num_face)*wall_length_face(num_face)+DMINFLOAT))
                          * exp(-0.00375*Rey*Rey);
    }

  // Calcul de F1
  for (num_face=0; num_face <nb_faces; num_face ++ )
    if (P(num_face) > 0)
      F1[num_face] = 1.+ Pderive(num_face) / (P(num_face));
    else
      F1[num_face] = 1.;
  return F1;
}

DoubleTab& Modele_EASM_Baglietto_VEF::Calcul_F2( DoubleTab& F2, DoubleTab& Deb, const Zone_dis& zone_dis,const DoubleTab& K_eps_Bas_Re,const Champ_base& ch_visco ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco.valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco);
  if (is_visco_const)
    visco=tab_visco(0,0);
  const Zone_VEF& la_zone = ref_cast(Zone_VEF,zone_dis.valeur());
  int nb_faces = la_zone.nb_faces();
  int num_face;
  double Re;

  /*
    for (num_face=0; num_face<nb_faces  ; num_face++)
  	  F2[num_face] = 1.;
    return F2;
  */

  for (num_face=0; num_face<nb_faces  ; num_face++)
    {
      if (!is_visco_const)
        {
          int elem0 = la_zone.face_voisins(num_face,0);
          int elem1 = la_zone.face_voisins(num_face,1);
          if (elem1!=-1)
            {
              visco = tab_visco(elem0)*la_zone.volumes(elem0)+tab_visco(elem1)*la_zone.volumes(elem1);
              visco /= la_zone.volumes(elem0) + la_zone.volumes(elem1);
            }
          else
            visco =  tab_visco(elem0);
        }
      if (visco>0 && K_eps_Bas_Re(num_face,1)>0)
        {
          Re = (K_eps_Bas_Re(num_face,0)*K_eps_Bas_Re(num_face,0))/(visco*K_eps_Bas_Re(num_face,1));
          F2[num_face] = 1. - (0.3*exp(-Re*Re));
        }
      else
        F2[num_face] = 1.;
    }
  //Cerr<<F2.mp_min_vect()<<" F2 "<<F2.mp_max_vect()<<finl;
  return F2;
}
/*
  DoubleTab& Modele_EASM_Baglietto_VEF::Calcul_F2( DoubleTab& F2, DoubleTab& D, const Zone_dis& zone_dis,const DoubleTab& K_eps_Bas_Re, const DoubleTab& tab_visco ) const
  {
  const Zone_VEF& la_zone = ref_cast(Zone_VEF,zone_dis.valeur());
  int nb_faces = la_zone.nb_faces();
  int num_face,elem0,elem1;
  double Re,nulam;

  for (num_face=0; num_face<nb_faces  ; num_face++)
  {
  elem0 = la_zone.face_voisins(num_face,0);
  elem1 = la_zone.face_voisins(num_face,1);
  if (elem1!=-1)
  {
  nulam = tab_visco(elem0)*la_zone.volumes(elem0)+tab_visco(elem1)*la_zone.volumes(elem1);
  nulam /= la_zone.volumes(elem0) + la_zone.volumes(elem1);
  }
  else
  nulam =  tab_visco(elem0);

  if (K_eps_Bas_Re(num_face,1)>0)
  {
  Re = (K_eps_Bas_Re(num_face,0)*K_eps_Bas_Re(num_face,0))/(nulam*K_eps_Bas_Re(num_face,1));
  F2[num_face] = 1. - (0.3*exp(-Re*Re));
  }
  else
  {
  F2[num_face] = 1.;
  }
  }
  return F2;
  }

*/
DoubleTab&  Modele_EASM_Baglietto_VEF::Calcul_Fmu( DoubleTab& Fmu,const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis, const DoubleTab& K_eps_Bas_Re,const Champ_Don& ch_visco ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco.valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco.valeur());
  if (is_visco_const)
    visco=tab_visco(0,0);
  const Zone_VEF& la_zone = ref_cast(Zone_VEF,zone_dis.valeur());
  const Zone_Cl_VEF& zone_Cl_VEF = ref_cast(Zone_Cl_VEF,zone_Cl_dis.valeur());
  int nb_faces = la_zone.nb_faces();
  int num_face;
  double Rey;
  const DoubleTab& wall_length = BR_wall_length_.valeurs();
  DoubleTab wall_length_face(0);
  la_zone.creer_tableau_faces(wall_length_face);
  const Conds_lim& les_cl = zone_Cl_VEF.les_conditions_limites();
  int nb_cl=les_cl.size();
  const IntTab& face_voisins = la_zone.face_voisins();
  Fmu = 0;

  // Calcul de la distance a la paroi aux faces
  for (int n_bord=0; n_bord<nb_cl; n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
      int ndeb = le_bord.num_premiere_face();
      int nfin = ndeb + le_bord.nb_faces();

      if (sub_type(Periodique,la_cl.valeur()))
        {
          for (num_face=ndeb; num_face<nfin; num_face++)
            {
              int elem1 = face_voisins(num_face,0);
              int elem2 = face_voisins(num_face,1);
              wall_length_face(num_face) = 0.5*wall_length(elem1) + 0.5*wall_length(elem2);
            }
        }
      /* else if (sub_type(Pa,la_cl.valeur()))
      {
        for (int num_face=ndeb; num_face<nfin; num_face++)
        {
      	  int elem1 = face_voisins(num_face,0);
      	  int elem2 = face_voisins(num_face,1);
      			  wall_length_face(num_face) = 0.;
        }
      }*/
      else
        {
          for (num_face=ndeb; num_face<nfin; num_face++)
            {
              int elem1 = face_voisins(num_face,0);
              wall_length_face(num_face) = wall_length(elem1);
            }
        }
    }
  int n0 = la_zone.premiere_face_int();
  for (num_face=n0; num_face<nb_faces; num_face++)
    {
      int elem0 = la_zone.face_voisins(num_face,0);
      int elem1 = la_zone.face_voisins(num_face,1);
      wall_length_face(num_face) = 0.5*wall_length(elem0)+0.5*wall_length(elem1);
    }
//  Cerr<<wall_length_face.mp_min_vect()<<" wall_length " <<wall_length_face.mp_max_vect()<<finl;
  /*for (num_face=0; num_face< la_zone.premiere_face_int(); num_face++)
   {
    int elem0 = la_zone.face_voisins(num_face,0);
    if (elem0 != -1)
  	  wall_length_face(num_face) = wall_length(elem0);
    else
    {
  	  elem0 = la_zone.face_voisins(num_face,1);
  	  wall_length_face(num_face) = wall_length(elem0);
    }
   }

   for (; num_face<nb_faces; num_face++)
   {
    int elem0 = la_zone.face_voisins(num_face,0);
    int elem1 = la_zone.face_voisins(num_face,1);
    wall_length_face(num_face) = 0.5*wall_length(elem0)+0.5*wall_length(elem1);
   }*/

// Calcul de Fmu
  for (num_face=0; num_face<nb_faces  ; num_face++)
    {
      if (!is_visco_const)
        {
          int elem0 = la_zone.face_voisins(num_face,0);
          int elem1 = la_zone.face_voisins(num_face,1);
          if (elem1!=-1)
            {
              visco = tab_visco(elem0)*la_zone.volumes(elem0)+tab_visco(elem1)*la_zone.volumes(elem1);
              visco /= la_zone.volumes(elem0) + la_zone.volumes(elem1);
            }
          else
            visco =  tab_visco(elem0);
        }
      if (visco>0)
        {
          Rey = wall_length_face(num_face)*sqrt(K_eps_Bas_Re(num_face,0))/visco;
          Fmu[num_face] = 1.0 - exp(-0.029*pow(Rey,0.5) - 0.00011*Rey*Rey);
        }
      else
        Fmu(num_face) = 1.;
    }
  return Fmu;
}

