/****************************************************************************
* Copyright (c) 2019, CEA
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
// File:        Source_Transport_K_Eps_Bas_Reynolds_W_VDF_Elem.cpp
// Directory:   $TRUST_ROOT/src/VDF/Turbulence
// Version:     /main/22
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Transport_K_Eps_Bas_Reynolds_W_VDF_Elem.h>
#include <Modele_turbulence_scal_Fluctuation_Temperature_W_Bas_Re.h>
#include <Convection_Diffusion_Temperature.h>
#include <Probleme_base.h>
#include <IntTrav.h>
#include <Entree_fluide_temperature_imposee.h>
#include <Champ_Uniforme.h>
#include <Champ_Face.h>
#include <Zone_VDF.h>
#include <Zone_Cl_VDF.h>
#include <Modele_turbulence_hyd_K_Eps_Bas_Reynolds.h>
#include <DoubleTrav.h>
#include <Fluide_Quasi_Compressible.h>

Implemente_instanciable_sans_constructeur(Source_Transport_K_Eps_Bas_Reynolds_W_VDF_Elem,"Source_Transport_K_Eps_Bas_Reynolds_W_VDF_P0_VDF",Source_Transport_K_Eps_VDF_Elem);





//// printOn
//

Sortie& Source_Transport_K_Eps_Bas_Reynolds_W_VDF_Elem::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}



//// readOn
//

Entree& Source_Transport_K_Eps_Bas_Reynolds_W_VDF_Elem::readOn(Entree& is )
{
  Motcle accolade_ouverte("{");
  Motcle accolade_fermee("}");
  Motcle motlu;

  is >> motlu;
  if (motlu != accolade_ouverte)
    {
      Cerr << "On attendait { pour commencer a lire les constantes de Source_Transport_K_Eps_Bas_Reynolds_W" << finl;
      exit();
    }
  Cerr << "Lecture des constantes de Source_Transport_K_Eps_Bas_Reynolds_W" << finl;
  Motcles les_mots(2);
  {
    les_mots[0] = "C1_eps";
    les_mots[1] = "C2_eps";
  }
  is >> motlu;
  while (motlu != accolade_fermee)
    {
      int rang=les_mots.search(motlu);
      switch(rang)
        {
        case 0 :
          {
            is >> C1;
            break;
          }
        case 1 :
          {
            is >> C2;
            break;
          }
        default :
          {
            Cerr << "On ne comprend pas le mot cle : " << motlu << "dans Source_Transport_K_Eps_Bas_Reynolds_W" << finl;
            exit();
          }
        }

      is >> motlu;
    }
  return is;
}






/////////////////////////////////////////////////////////////////////////////
//
//            Implementation des fonctions de la classe
//
//             Source_Transport_K_Eps_Bas_Reynolds_W_VDF_Elem
//
/////////////////////////////////////////////////////////////////////////////

void Source_Transport_K_Eps_Bas_Reynolds_W_VDF_Elem::associer_pb(const Probleme_base& pb )
{
  eq_hydraulique = pb.equation(0);
  eqn_keps_bas_re = ref_cast(Transport_K_Eps_Bas_Reynolds,equation());
}

void Source_Transport_K_Eps_Bas_Reynolds_W_VDF_Elem::associer_zones(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis)
{
  la_zone_VDF = ref_cast(Zone_VDF, zone_dis.valeur());
  la_zone_Cl_VDF = ref_cast(Zone_Cl_VDF, zone_Cl_dis.valeur());
}

DoubleTab& Source_Transport_K_Eps_Bas_Reynolds_W_VDF_Elem::ajouter(DoubleTab& resu) const
{

  const Zone_Cl_dis& zcl=eq_hydraulique->zone_Cl_dis();
  const Zone_dis& z = eq_hydraulique->zone_dis();
  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const Zone_Cl_VDF& zone_Cl_VDF = la_zone_Cl_VDF.valeur();
  const  DoubleTab& K_eps_Bas_Re = eqn_keps_bas_re->inconnue().valeurs();
  const DoubleTab& visco_turb = eqn_keps_bas_re->modele_turbulence().viscosite_turbulente().valeurs();
  const Modele_turbulence_hyd_K_Eps_Bas_Reynolds& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_Bas_Reynolds,eqn_keps_bas_re->modele_turbulence());
  const Modele_Fonc_Bas_Reynolds& mon_modele_fonc = mod_turb.associe_modele_fonction();
  const Fluide_base& fluide = ref_cast(Fluide_base,eq_hydraulique->milieu());
  const Champ_Don& ch_visco_cin = fluide.viscosite_cinematique();
  /*
    const DoubleTab& tab_visco = ch_visco_cin->valeurs();
    double visco;
    if (sub_type(Champ_Uniforme,ch_visco_cin.valeur()))
    {
    visco = std::max(tab_visco(0,0),DMINFLOAT);
    }
    else {assert(0);exit();visco=-1;}
  */
  const DoubleTab& vit = eq_hydraulique->inconnue().valeurs();
  const DoubleVect& volumes = zone_VDF.volumes();
  const DoubleVect& porosite_vol = zone_VDF.porosite_elem();
  int nb_elem = zone_VDF.nb_elem();
  //  int nb_elem_tot = zone_VDF.nb_elem_tot();
  DoubleTrav P(nb_elem);
  DoubleTrav D(nb_elem);
  DoubleTrav E(nb_elem);
  DoubleTrav F1(nb_elem);
  DoubleTrav F2(nb_elem);
  // Cerr << "le modele fonc est " << mon_modele_fonc.valeur().que_suis_je() << finl;

  mon_modele_fonc.Calcul_D(D,z,zcl,vit,K_eps_Bas_Re,ch_visco_cin);
  mon_modele_fonc.Calcul_E(E,z,zcl,vit,K_eps_Bas_Re,ch_visco_cin,visco_turb);
  //mon_modele_fonc.Calcul_F1(F1,z);
  mon_modele_fonc.Calcul_F2(F2,D,z,K_eps_Bas_Re,ch_visco_cin );
  if (axi)
    {
      Champ_Face& vitesse = ref_cast_non_const(Champ_Face,eq_hydraulique->inconnue().valeur());
      calculer_terme_production_K_Axi(zone_VDF,vitesse,P,K_eps_Bas_Re,visco_turb);
    }
  else
    {
      Champ_Face& vitesse = ref_cast_non_const(Champ_Face,eq_hydraulique->inconnue().valeur());
      calculer_terme_production_K(zone_VDF,zone_Cl_VDF ,P,K_eps_Bas_Re,vit,vitesse,visco_turb);
    }

  for (int elem=0; elem<nb_elem; elem++)
    {
      resu(elem,0) += (P(elem)-K_eps_Bas_Re(elem,1)-D(elem))*volumes(elem)*porosite_vol(elem);

      //   Cerr << "K = " << K_eps_Bas_Re(elem,0) << finl;
      if (K_eps_Bas_Re(elem,0) >= 10.e-10)
        resu(elem,1) += (C1*F1(elem)*P(elem)- C2*F2(elem)*K_eps_Bas_Re(elem,1))*volumes(elem)*porosite_vol(elem)
                        *K_eps_Bas_Re(elem,1)/K_eps_Bas_Re(elem,0)+E(elem)*volumes(elem)*porosite_vol(elem);
    }
  //   Cerr << "resu = " << resu(elem,1) << finl;
  return resu;
}

DoubleTab& Source_Transport_K_Eps_Bas_Reynolds_W_VDF_Elem::calculer(DoubleTab& resu) const
{
  resu = 0.;
  return ajouter(resu);
}

/////////////////////////////////////////////////////////////////////////////
//
//            Implementation des fonctions de la classe
//
//             Source_Transport_K_Eps_Bas_Reynolds_anisotherme_W_VDF_Elem
//
/////////////////////////////////////////////////////////////////////////////
//


