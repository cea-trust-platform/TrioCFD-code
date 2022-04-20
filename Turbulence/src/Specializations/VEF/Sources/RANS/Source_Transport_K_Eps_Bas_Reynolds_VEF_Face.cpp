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
// File:        Source_Transport_K_Eps_Bas_Reynolds_VEF_Face.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Sources
//
//////////////////////////////////////////////////////////////////////////////


#include <Source_Transport_K_Eps_Bas_Reynolds_VEF_Face.h>
#include <Convection_Diffusion_Temperature.h>
#include <Modele_turbulence_scal_base.h>
#include <Fluide_base.h>
#include <Probleme_base.h>
#include <Champ_Uniforme.h>
#include <Zone_VEF.h>
#include <Modele_turbulence_hyd_K_Eps_Bas_Reynolds.h>
#include <Champ_P1NC.h>
#include <Debog.h>
#include <TRUSTTrav.h>
#include <Fluide_Quasi_Compressible.h>

Implemente_instanciable_sans_constructeur(Source_Transport_K_Eps_Bas_Reynolds_VEF_Face,"Source_Transport_K_Eps_Bas_Reynolds_VEF_P1NC",Source_base);



//// printOn
//

Sortie& Source_Transport_K_Eps_Bas_Reynolds_VEF_Face::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}


//// readOn
//

Entree& Source_Transport_K_Eps_Bas_Reynolds_VEF_Face::readOn(Entree& is )
{
  Motcle accolade_ouverte("{");
  Motcle accolade_fermee("}");
  Motcle motlu;

  is >> motlu;
  if (motlu != accolade_ouverte)
    {
      Cerr << "On attendait { pour commencer a lire les constantes de Source_Transport_K_Eps_Bas_Reynolds" << finl;
      exit();
    }
  Cerr << "Lecture des constantes de Source_Transport_K_Eps_Bas_Reynolds" << finl;
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
            Cerr << "On ne comprend pas le mot cle : " << motlu << "dans Source_Transport_K_Eps_Bas_Reynolds" << finl;
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
//             Source_Transport_K_Eps_Bas_Reynolds_VEF_Face
//
/////////////////////////////////////////////////////////////////////////////

void Source_Transport_K_Eps_Bas_Reynolds_VEF_Face::associer_pb(const Probleme_base& pb )
{
  eq_hydraulique = pb.equation(0);
  eqn_keps_bas_re = ref_cast(Transport_K_Eps_Bas_Reynolds,equation());
}

void Source_Transport_K_Eps_Bas_Reynolds_VEF_Face::associer_zones(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis)
{
  la_zone_VEF = ref_cast(Zone_VEF, zone_dis.valeur());
  la_zone_Cl_VEF = ref_cast(Zone_Cl_VEF, zone_Cl_dis.valeur());
}

DoubleTab& Source_Transport_K_Eps_Bas_Reynolds_VEF_Face::ajouter(DoubleTab& resu) const
{
  Debog::verifier("Source_Transport_K_Eps_Bas_Reynolds_VEF_Face::ajouter resu 0",resu);
  const Zone_Cl_dis& zcl_keps=eqn_keps_bas_re->zone_Cl_dis();
  const Zone_dis& zone_dis_keps =eqn_keps_bas_re ->zone_dis();
  const Zone_VEF& zone_VEF = la_zone_VEF.valeur();
  const Zone_Cl_VEF& zone_Cl_VEF = la_zone_Cl_VEF.valeur();
  const DoubleTab& K_eps_Bas_Re = eqn_keps_bas_re->inconnue().valeurs();
  const Modele_turbulence_hyd_K_Eps_Bas_Reynolds& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_Bas_Reynolds,eqn_keps_bas_re->modele_turbulence());
  const DoubleTab& visco_turb = mod_turb.viscosite_turbulente().valeurs();
  const Modele_Fonc_Bas_Reynolds& mon_modele_fonc = mod_turb.associe_modele_fonction();
  const Fluide_base& fluide = ref_cast(Fluide_base,eq_hydraulique->milieu());
  const Champ_Don& ch_visco_cin = fluide.viscosite_cinematique();
  //const DoubleTab& K_eps = mon_eq_transport_K_Eps->inconnue().valeurs();
  /*
    const DoubleTab& tab_visco = ch_visco_cin->valeurs();
    double visco=-1;
    if (sub_type(Champ_Uniforme,ch_visco_cin.valeur()))
    {
    visco = std::max(tab_visco(0,0),DMINFLOAT);
    }
  */
  const DoubleTab& vit = eq_hydraulique->inconnue().valeurs();
  const DoubleVect& vol_ent = zone_VEF.volumes_entrelaces();

  int nb_faces = zone_VEF.nb_faces();
  //  int nb_faces_tot = zone_VEF.nb_faces_tot();
  DoubleTrav P(nb_faces);
  DoubleTrav D(vol_ent);
  DoubleTrav E(vol_ent);
  DoubleTrav F1(nb_faces);
  DoubleTrav F2(nb_faces);

  mon_modele_fonc.Calcul_D(D,zone_dis_keps,zcl_keps,vit,K_eps_Bas_Re,ch_visco_cin);
  D.echange_espace_virtuel();
  mon_modele_fonc.Calcul_E(E,zone_dis_keps,zcl_keps,vit,K_eps_Bas_Re,ch_visco_cin,visco_turb);
  //const Champ_base& ch_visco_cin_ou_dyn =ref_cast(Op_Diff_K_Eps_base, equation().operateur(0).l_op_base()).diffusivite();
  //mon_modele_fonc.Calcul_F1(F1,zone_dis_keps,zcl_keps, P, K_eps,ch_visco_cin_ou_dyn);
  mon_modele_fonc.Calcul_F2(F2,D,zone_dis_keps,K_eps_Bas_Re,ch_visco_cin);

  calculer_terme_production_K(zone_VEF,zone_Cl_VEF,P,K_eps_Bas_Re,vit,visco_turb);
  // Ajout des termes sources

  // for (int num_face=0; num_face<nb_faces; num_face++)
  //     {
  //       resu(num_face,0) += (P(num_face)-K_eps_Bas_Re(num_face,1)-D(num_face))*vol_ent(num_face);
  //        if (K_eps_Bas_Re(num_face,0) >= 1.e-10)
  //          resu(num_face,1) += ((C1*F1(num_face)*P(num_face)- C2*F2(num_face)*K_eps_Bas_Re(num_face,1))*K_eps_Bas_Re(num_face,1)/K_eps_Bas_Re(num_face,0)+E(num_face))*vol_ent(num_face);
  //     }


  Debog::verifier("Source_Transport_K_Eps_Bas_Reynolds_VEF_Face::ajouter P 0",P);
  Debog::verifier("Source_Transport_K_Eps_Bas_Reynolds_VEF_Face::ajouter D 0",D);
  Debog::verifier("Source_Transport_K_Eps_Bas_Reynolds_VEF_Face::ajouter E 0",E);
  Debog::verifier("Source_Transport_K_Eps_Bas_Reynolds_VEF_Face::ajouter F1 0",F1);
  Debog::verifier("Source_Transport_K_Eps_Bas_Reynolds_VEF_Face::ajouter F2 0",F2);


  for (int num_face=0; num_face<nb_faces; num_face++)
    {
      if (K_eps_Bas_Re(num_face,0) >= 1.e-20 && K_eps_Bas_Re(num_face,1) > 1.e-20)
        {
          resu(num_face,0) += (P(num_face)-K_eps_Bas_Re(num_face,1)-D(num_face))*vol_ent(num_face);

          resu(num_face,1) += ((C1*F1(num_face)*P(num_face)- C2*F2(num_face)*K_eps_Bas_Re(num_face,1))*K_eps_Bas_Re(num_face,1)/K_eps_Bas_Re(num_face,0)+E(num_face))*vol_ent(num_face);
        }
      else
        {
          resu(num_face,0) += 0.;
          resu(num_face,1) += 0.;
        }
    }


  return resu;
}

DoubleTab& Source_Transport_K_Eps_Bas_Reynolds_VEF_Face::calculer(DoubleTab& resu) const
{
  resu = 0.;
  return ajouter(resu);
}

