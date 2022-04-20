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
// File:        Source_Transport_K_Eps_Realisable_VEF_Face.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Sources
//
//////////////////////////////////////////////////////////////////////////////


#include <Source_Transport_K_Eps_Realisable_VEF_Face.h>
#include <Convection_Diffusion_Temperature.h>
#include <Convection_Diffusion_Concentration.h>
#include <Modele_turbulence_scal_base.h>
#include <Fluide_base.h>
#include <Probleme_base.h>
#include <Champ_Uniforme.h>
#include <Zone_VEF.h>
#include <Modele_turbulence_hyd_K_Eps_Realisable.h>
#include <Champ_P1NC.h>
#include <Debog.h>
#include <TRUSTTrav.h>
#include <Fluide_Quasi_Compressible.h>
#include <Modele_Shih_Zhu_Lumley_VEF.h>
#include <Pb_Hydraulique_Turbulent.h>
#include <Pb_Hydraulique_Concentration_Turbulent.h>
#include <Pb_Thermohydraulique_Turbulent_QC.h>
#include <Pb_Thermohydraulique_Turbulent.h>
#include <Pb_Thermohydraulique_Concentration_Turbulent.h>
#include <Constituant.h>

Implemente_instanciable_sans_constructeur(Source_Transport_K_Eps_Realisable_VEF_Face,"Source_Transport_K_Eps_Realisable_VEF_P1NC",Source_base);



//// printOn
//

Sortie& Source_Transport_K_Eps_Realisable_VEF_Face::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}





// readOn


Entree& Source_Transport_K_Eps_Realisable_VEF_Face::readOn(Entree& is )
{
  Motcle accolade_ouverte("{");
  Motcle accolade_fermee("}");
  Motcle motlu;

  is >> motlu;
  if (motlu != accolade_ouverte)
    {
      Cerr << "On attendait { pour commencer a lire les constantes de Source_Transport_K_Eps_Realisable" << finl;
      exit();
    }
  Cerr << "Lecture des constantes de Source_Transport_K_Eps_Realisable" << finl;
  Motcles les_mots(1);
  {
    les_mots[0] = "C2_eps";
  }
  is >> motlu;
  while (motlu != accolade_fermee)
    {
      int rang=les_mots.search(motlu);
      switch(rang)
        {
        case 0 :
          {
            is >> C2_;
            break;
          }
        default :
          {
            Cerr << "On ne comprend pas le mot cle : " << motlu << "dans Source_Transport_K_Eps_Realisable" << finl;
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
//             Source_Transport_K_Eps_Realisable_VEF_Face
//
/////////////////////////////////////////////////////////////////////////////

void Source_Transport_K_Eps_Realisable_VEF_Face::associer_pb(const Probleme_base& pb )
{
  eq_hydraulique = pb.equation(0);
  eqn_keps_Rea = ref_cast(Transport_K_Eps_Realisable,equation());
}

void Source_Transport_K_Eps_Realisable_VEF_Face::associer_zones(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis)
{
  la_zone_VEF = ref_cast(Zone_VEF, zone_dis.valeur());
  la_zone_Cl_VEF = ref_cast(Zone_Cl_VEF, zone_Cl_dis.valeur());
}

DoubleTab& Source_Transport_K_Eps_Realisable_VEF_Face::ajouter(DoubleTab& resu) const
{
  Debog::verifier("Source_Transport_K_Eps_Realisable_VEF_Face::ajouter resu 0",resu);

  const Zone_VEF& zone_VEF                               = la_zone_VEF.valeur();
  const Zone_Cl_VEF& zone_Cl_VEF                         = la_zone_Cl_VEF.valeur();
  const DoubleTab& K_eps_Rea                             = eqn_keps_Rea->inconnue().valeurs();
  const Modele_turbulence_hyd_K_Eps_Realisable& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_Realisable,eqn_keps_Rea->modele_turbulence());
  const DoubleTab& visco_turb                            = mod_turb.viscosite_turbulente().valeurs();
  const Modele_Fonc_Realisable_base& mon_modele_fonc     = mod_turb.associe_modele_fonction();
  const Fluide_base& fluide                    = ref_cast(Fluide_base,eq_hydraulique->milieu());

  const Champ_Don& ch_visco_cin = fluide.viscosite_cinematique();
  const DoubleTab& tab_visco    = ch_visco_cin->valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco_cin.valeur());
  double visco=-1;
  if (is_visco_const)
    {
      visco = std::max(tab_visco(0,0),DMINFLOAT);
    }

  const DoubleTab& vit = eq_hydraulique->inconnue().valeurs();
  const DoubleVect& vol_ent = zone_VEF.volumes_entrelaces();

  DoubleTab vitesse_filtree(vit);
  ref_cast(Champ_P1NC,eq_hydraulique->inconnue().valeur()).filtrer_L2(vitesse_filtree);

  double LeK_MIN = mod_turb.get_LeK_MIN();
  double LeEPS_MIN = mod_turb.get_LeEPS_MIN();

  int nb_faces = zone_VEF.nb_faces();
  //  int nb_faces_tot = zone_VEF.nb_faces_tot();
  DoubleTrav P(nb_faces);

  const DoubleTab&  C1 = mon_modele_fonc.get_C1( );
  const DoubleTab&  S  = mon_modele_fonc.get_S( );
//   C1.echange_espace_virtuel();

//   calculer_terme_production_K(zone_VEF,zone_Cl_VEF,P,K_eps_Rea,vit,visco_turb);
  calculer_terme_production_K(zone_VEF,zone_Cl_VEF,P,K_eps_Rea,vitesse_filtree,visco_turb);

  // Ajout des termes sources

  Debog::verifier("Source_Transport_K_Eps_Realisable_VEF_Face::ajouter P 0",P);
  Debog::verifier("Source_Transport_K_Eps_Realisable_VEF_Face::ajouter C1 0",C1);

  for (int num_face=0; num_face<nb_faces; num_face++)
    {
      if (!is_visco_const)
        {
          int elem0 = zone_VEF.face_voisins(num_face,0);
          int elem1 = zone_VEF.face_voisins(num_face,1);
          if (elem1!=-1)
            {
              visco = tab_visco(elem0)*zone_VEF.volumes(elem0)+tab_visco(elem1)*zone_VEF.volumes(elem1);
              visco /= zone_VEF.volumes(elem0) + zone_VEF.volumes(elem1);
            }
          else
            {
              visco =  tab_visco(elem0);
            }
        }

      assert(visco>0.);

      resu(num_face,0) += ( P(num_face)-K_eps_Rea(num_face,1) )*vol_ent(num_face);

      if ( ( K_eps_Rea(num_face,0) >= LeK_MIN ) and ( K_eps_Rea(num_face,1) >= LeEPS_MIN ) )
        {
          resu(num_face,1) += K_eps_Rea(num_face,1)*( C1(num_face)*S(num_face)  - ( C2_*K_eps_Rea(num_face,1)/( K_eps_Rea(num_face,0) + sqrt(  visco*K_eps_Rea(num_face,1) ) ) ) )*vol_ent(num_face);
        }
    }

  return resu;
}


DoubleTab& Source_Transport_K_Eps_Realisable_VEF_Face::calculer(DoubleTab& resu) const
{
  resu = 0.;
  return ajouter(resu);
}

void Source_Transport_K_Eps_Realisable_VEF_Face::mettre_a_jour(double temps)
{
  const Zone_Cl_dis& zcl_keps                            = eqn_keps_Rea->zone_Cl_dis();
  const Zone_dis& zone_dis_keps                          = eqn_keps_Rea ->zone_dis();
  const DoubleTab& K_eps_Rea                             = eqn_keps_Rea->inconnue().valeurs();

  Modele_turbulence_hyd_K_Eps_Realisable& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_Realisable,eqn_keps_Rea->modele_turbulence());
  Modele_Fonc_Realisable_base& mon_modele_fonc     = mod_turb.associe_modele_fonction();
  const DoubleTab& visco_turb                      = mod_turb.viscosite_turbulente().valeurs();

  const DoubleTab& vit                                   = eq_hydraulique->inconnue().valeurs();
  double epsilon_minimum                                 = eqn_keps_Rea.valeur().modele_turbulence().get_LeEPS_MIN();

//   DoubleTab vitesse_filtree(vit);
//   ref_cast(Champ_P1NC,eq_hydraulique->inconnue().valeur()).filtrer_L2(vitesse_filtree);

  const Champ_Don ch_visco      = ref_cast(Fluide_base,eqn_keps_Rea->milieu()).viscosite_cinematique();
  const Champ_Don& ch_visco_cin = ref_cast(Fluide_base,eqn_keps_Rea->milieu()).viscosite_cinematique();
  const DoubleTab& tab_visco    = ch_visco_cin->valeurs();

  /*Paroi*/
  Nom lp=mod_turb.loi_paroi().valeur().que_suis_je();
  DoubleTab visco_tab(visco_turb.dimension_tot(0));
  assert(sub_type(Champ_Uniforme,ch_visco_cin.valeur()));
  visco_tab = tab_visco(0,0);
  const int idt =  eq_hydraulique->schema_temps().nb_pas_dt();
  const DoubleTab& tab_paroi = mod_turb.loi_paroi().valeur().Cisaillement_paroi();

//   mon_modele_fonc.Contributions_Sources_Paroi(zone_dis_keps,zcl_keps,vitesse_filtree,K_eps_Rea,epsilon_minimum,visco_tab,visco_turb,tab_paroi,idt);

  mon_modele_fonc.Contributions_Sources_Paroi(zone_dis_keps,zcl_keps,vit,K_eps_Rea,epsilon_minimum,visco_tab,visco_turb,tab_paroi,idt);

  Calcul_Production_K_VEF::mettre_a_jour(temps);
}

void Source_Transport_K_Eps_Realisable_VEF_Face::contribuer_a_avec(const DoubleTab& a, Matrice_Morse& matrice) const
{
  const  DoubleTab& K_eps_Rea                            = eqn_keps_Rea->inconnue().valeurs();
  int    size                                            = K_eps_Rea.dimension(0);
  const Modele_turbulence_hyd_K_Eps_Realisable& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_Realisable,eqn_keps_Rea->modele_turbulence());
  double LeK_MIN                                         = mod_turb.get_LeK_MIN();
  double LeEPS_MIN                                       = mod_turb.get_LeEPS_MIN();

  const Fluide_base& fluide                    = ref_cast(Fluide_base,eq_hydraulique->milieu());

  const Champ_Don& ch_visco_cin = fluide.viscosite_cinematique();
  const DoubleTab& tab_visco    = ch_visco_cin->valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco_cin.valeur());
  double visco=-1;
  if (is_visco_const)
    {
      visco = std::max(tab_visco(0,0),DMINFLOAT);
    }

  const Zone_VEF& zone_VEF = la_zone_VEF.valeur();
  const DoubleVect& porosite_face = zone_VEF.porosite_face();
  // on implicite le -eps et le -eps^2/k
  const DoubleVect& volumes_entrelaces=zone_VEF.volumes_entrelaces();
  for (int face=0; face<size; face++)
    {
      if (!is_visco_const)
        {
          int elem0 = zone_VEF.face_voisins(face,0);
          int elem1 = zone_VEF.face_voisins(face,1);
          if (elem1!=-1)
            {
              visco = tab_visco(elem0)*zone_VEF.volumes(elem0)+tab_visco(elem1)*zone_VEF.volumes(elem1);
              visco /= zone_VEF.volumes(elem0) + zone_VEF.volumes(elem1);
            }
          else
            {
              visco =  tab_visco(elem0);
            }
        }

      assert(visco>0.);


      // -eps*vol  donne +vol dans la bonne case
      if ( ( K_eps_Rea(face,0) >= LeK_MIN ) and ( K_eps_Rea(face,1) >= LeEPS_MIN ) )
        {
          double coef_k=porosite_face(face)*volumes_entrelaces(face)*K_eps_Rea(face,1)/( K_eps_Rea(face,0) + sqrt(  visco*K_eps_Rea(face,1) ) );
          matrice(face*2,face*2)+=coef_k;
          double coef_eps=C2_*K_eps_Rea(face,1)/( K_eps_Rea(face,0) + sqrt( visco*K_eps_Rea(face,1) ) )*volumes_entrelaces(face)*porosite_face(face);
          matrice(face*2+1,face*2+1)+=coef_eps;
        }
    }
}

