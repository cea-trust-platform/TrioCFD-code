/****************************************************************************
* Copyright (c) 2021, CEA
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
// File:        Source_Transport_K_Realisable_VEF_Face.cpp
// Directory:   $TRUST_ROOT/src/VEF/Turbulence
// Version:     /main/20
//
//////////////////////////////////////////////////////////////////////////////


#include <Source_Transport_K_Realisable_VEF_Face.h>
#include <Convection_Diffusion_Temperature.h>
#include <Convection_Diffusion_Concentration.h>
#include <Modele_turbulence_scal_base.h>
#include <Fluide_base.h>
#include <Probleme_base.h>
#include <Champ_Uniforme.h>
#include <Zone_VEF.h>
#include <Modele_turbulence_hyd_K_Eps_Realisable_Bicephale.h>
#include <Champ_P1NC.h>
#include <Debog.h>
#include <DoubleTrav.h>
#include <Fluide_Quasi_Compressible.h>
#include <Modele_Shih_Zhu_Lumley_VEF.h>
#include <Pb_Hydraulique_Turbulent.h>
#include <Pb_Hydraulique_Concentration_Turbulent.h>
#include <Pb_Thermohydraulique_Turbulent_QC.h>
#include <Pb_Thermohydraulique_Turbulent.h>
#include <Pb_Thermohydraulique_Concentration_Turbulent.h>
#include <Constituant.h>
#include <Param.h>

Implemente_instanciable(Source_Transport_K_Realisable_VEF_Face,"Source_Transport_K_Realisable_VEF_P1NC",Source_base);

Implemente_instanciable(Source_Transport_K_Realisable_anisotherme_VEF_Face,"Source_Transport_K_Realisable_anisotherme_VEF_P1NC",Source_Transport_K_Realisable_VEF_Face);
Implemente_instanciable(Source_Transport_K_Realisable_aniso_concen_VEF_Face,"Source_Transport_K_Realisable_aniso_concen_VEF_P1NC",Source_Transport_K_Realisable_VEF_Face);
Implemente_instanciable(Source_Transport_K_Realisable_aniso_therm_concen_VEF_Face,"Source_Transport_K_Realisable_aniso_therm_concen_VEF_P1NC",Source_Transport_K_Realisable_VEF_Face);


//// printOn
//

Sortie& Source_Transport_K_Realisable_VEF_Face::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}

Sortie& Source_Transport_K_Realisable_anisotherme_VEF_Face::printOn(Sortie& s) const
{
  return s << que_suis_je() ;
}

Sortie& Source_Transport_K_Realisable_aniso_concen_VEF_Face::printOn(Sortie& s) const
{
  return s << que_suis_je() ;
}

Sortie& Source_Transport_K_Realisable_aniso_therm_concen_VEF_Face::printOn(Sortie& s) const
{
  return s << que_suis_je() ;
}

// readOn


Entree& Source_Transport_K_Realisable_VEF_Face::readOn(Entree& is )
{
  Param param(que_suis_je());
  param.lire_avec_accolades(is);
  return is ;
}

Entree& Source_Transport_K_Realisable_anisotherme_VEF_Face::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.lire_avec_accolades(is);
  return is ;
}

Entree& Source_Transport_K_Realisable_aniso_concen_VEF_Face::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.lire_avec_accolades(is);
  return is ;
}

Entree& Source_Transport_K_Realisable_aniso_therm_concen_VEF_Face::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.lire_avec_accolades(is);
  return is ;
}


/////////////////////////////////////////////////////////////////////////////
//
//            Implementation des fonctions de la classe
//
//             Source_Transport_K_Realisable_VEF_Face
//
/////////////////////////////////////////////////////////////////////////////

void Source_Transport_K_Realisable_VEF_Face::associer_pb(const Probleme_base& pb )
{
  eq_hydraulique = pb.equation(0);
  eqn_k_Rea      = ref_cast(Transport_K_ou_Eps_Realisable,equation());
  eqn_eps_Rea    = ref_cast(Transport_K_ou_Eps_Realisable,eqn_eps_Rea.valeur().modele_turbulence().eqn_transp_Eps());
}

void Source_Transport_K_Realisable_VEF_Face::associer_zones(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis)
{
  la_zone_VEF = ref_cast(Zone_VEF, zone_dis.valeur());
  la_zone_Cl_VEF = ref_cast(Zone_Cl_VEF, zone_Cl_dis.valeur());
}

DoubleTab& Source_Transport_K_Realisable_VEF_Face::ajouter(DoubleTab& resu) const
{
  Debog::verifier("Source_Transport_K_Realisable_VEF_Face::ajouter resu 0",resu);

  const Zone_VEF& zone_VEF                               = la_zone_VEF.valeur();
  const Zone_Cl_VEF& zone_Cl_VEF                         = la_zone_Cl_VEF.valeur();

  const DoubleTab&                                           K_Rea = eqn_k_Rea->inconnue().valeurs();
  const DoubleTab&                                         eps_Rea = eqn_eps_Rea->inconnue().valeurs();
  const Modele_turbulence_hyd_K_Eps_Realisable_Bicephale& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_Realisable_Bicephale,eqn_k_Rea->modele_turbulence());

  const DoubleTab& visco_turb                            = mod_turb.viscosite_turbulente().valeurs();

  const DoubleTab& vit = eq_hydraulique->inconnue().valeurs();
  const DoubleVect& vol_ent = zone_VEF.volumes_entrelaces();

  DoubleTab vitesse_filtree(vit);
  ref_cast(Champ_P1NC,eq_hydraulique->inconnue().valeur()).filtrer_L2(vitesse_filtree);

  int nb_faces = zone_VEF.nb_faces();
  //  int nb_faces_tot = zone_VEF.nb_faces_tot();
  DoubleTrav P(nb_faces);

//   calculer_terme_production_K(zone_VEF,zone_Cl_VEF,P,K_eps_Rea,vit,visco_turb);
  calculer_terme_production_K_BiK(zone_VEF,zone_Cl_VEF,P,K_Rea,eps_Rea,vitesse_filtree,visco_turb);

  // Ajout des termes sources

  Debog::verifier("Source_Transport_K_Realisable_VEF_Face::ajouter P 0",P);

  for (int num_face=0; num_face<nb_faces; num_face++)
    {
      resu(num_face) += ( P(num_face)-eps_Rea(num_face) )*vol_ent(num_face);
    }

  return resu;
}


DoubleTab& Source_Transport_K_Realisable_VEF_Face::calculer(DoubleTab& resu) const
{
  resu = 0.;
  return ajouter(resu);
}

void Source_Transport_K_Realisable_VEF_Face::mettre_a_jour(double temps)
{
  const Zone_Cl_dis& zcl_keps                            = eqn_k_Rea->zone_Cl_dis();
  const Zone_dis& zone_dis_keps                          = eqn_k_Rea ->zone_dis();
  const DoubleTab& K_Rea                                 = eqn_k_Rea->inconnue().valeurs();
  const DoubleTab& eps_Rea                               = eqn_eps_Rea->inconnue().valeurs();

  Modele_turbulence_hyd_K_Eps_Realisable_Bicephale& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_Realisable_Bicephale,eqn_k_Rea->modele_turbulence());
  Modele_Fonc_Realisable_base& mon_modele_fonc               = mod_turb.associe_modele_fonction();
  const DoubleTab& visco_turb                                = mod_turb.viscosite_turbulente().valeurs();

  const DoubleTab& vit                                   = eq_hydraulique->inconnue().valeurs();
  double epsilon_minimum                                 = eqn_k_Rea.valeur().modele_turbulence().get_LeEPS_MIN();

//   DoubleTab vitesse_filtree(vit);
//   ref_cast(Champ_P1NC,eq_hydraulique->inconnue().valeur()).filtrer_L2(vitesse_filtree);

  const Champ_Don ch_visco      = ref_cast(Fluide_base,eqn_k_Rea->milieu()).viscosite_cinematique();
  const Champ_Don& ch_visco_cin = ref_cast(Fluide_base,eqn_k_Rea->milieu()).viscosite_cinematique();
  const DoubleTab& tab_visco    = ch_visco_cin->valeurs();

  /*Paroi*/
  Nom lp=mod_turb.loi_paroi().valeur().que_suis_je();
  DoubleTab visco_tab(visco_turb.dimension_tot(0));
  assert(sub_type(Champ_Uniforme,ch_visco_cin.valeur()));
  visco_tab = tab_visco(0,0);
  const int idt =  eq_hydraulique->schema_temps().nb_pas_dt();
  const DoubleTab& tab_paroi = mod_turb.loi_paroi().valeur().Cisaillement_paroi();

  mon_modele_fonc.Contributions_Sources_Paroi_BiK(zone_dis_keps,zcl_keps,vit,K_Rea,eps_Rea,epsilon_minimum,visco_tab,visco_turb,tab_paroi,idt);

  Calcul_Production_K_VEF::mettre_a_jour(temps);
}

void Source_Transport_K_Realisable_VEF_Face::contribuer_a_avec(const DoubleTab& a, Matrice_Morse& matrice) const
{
  const DoubleTab& K_Rea                                 = eqn_k_Rea->inconnue().valeurs();
  const DoubleTab& eps_Rea                               = eqn_eps_Rea->inconnue().valeurs();
  int    size                                            = K_Rea.dimension(0);
  const Modele_turbulence_hyd_K_Eps_Realisable_Bicephale& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_Realisable_Bicephale,eqn_k_Rea->modele_turbulence());
  double LeK_MIN                                         = mod_turb.get_LeK_MIN();
  double LeEPS_MIN                                       = mod_turb.get_LeEPS_MIN();

  const Fluide_base& fluide                    = ref_cast(Fluide_base,eq_hydraulique->milieu());

  const Champ_Don& ch_visco_cin = fluide.viscosite_cinematique();
  const DoubleTab& tab_visco    = ch_visco_cin->valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco_cin.valeur());
  double visco=-1;
  if (is_visco_const)
    {
      visco = max(tab_visco(0,0),DMINFLOAT);
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
      if ( ( K_Rea(face) >= LeK_MIN ) and ( eps_Rea(face) >= LeEPS_MIN ) )
        {
          matrice(face,face)+=porosite_face(face)*volumes_entrelaces(face)*eps_Rea(face)/( K_Rea(face) + sqrt(  visco*eps_Rea(face) ) );
        }
    }
}

DoubleTab& Source_Transport_K_Realisable_anisotherme_VEF_Face::ajouter(DoubleTab& resu) const
{
  Source_Transport_K_Realisable_VEF_Face::ajouter(resu);
  //
  // Plutost que de calculer P, on appelle Source_Transport_K_Realisable_VEF_Face::ajouter(resu)
  // et on ajoute directement G
  //
  const Zone_VEF& zone_VEF = la_zone_VEF.valeur();
  const Zone_Cl_VEF& zcl_VEF_th = ref_cast(Zone_Cl_VEF,eq_thermique->zone_Cl_dis().valeur());
  const DoubleTab& scalaire = eq_thermique->inconnue().valeurs();
  const Modele_turbulence_scal_base& le_modele_scalaire = ref_cast(Modele_turbulence_scal_base,eq_thermique->get_modele(TURBULENCE).valeur());
  const DoubleTab& alpha_turb = le_modele_scalaire.diffusivite_turbulente().valeurs();
  const DoubleTab& g = gravite->valeurs();
  const Champ_Don& ch_beta = beta_t.valeur();
  const DoubleVect& volumes_entrelaces = zone_VEF.volumes_entrelaces();

  int nb_face = zone_VEF.nb_faces();
  DoubleTrav G(nb_face);

  // C'est l'objet de type zone_Cl_dis de l'equation thermique
  // qui est utilise dans le calcul de G
  calculer_terme_destruction_K_gen(zone_VEF,zcl_VEF_th,G,scalaire,alpha_turb,ch_beta,g,0);

  for (int face=0; face<nb_face; face++)
    {
      resu(face) += G(face)*volumes_entrelaces(face);
    }
  return resu;
}

DoubleTab& Source_Transport_K_Realisable_aniso_concen_VEF_Face::ajouter(DoubleTab& resu) const
{
  Source_Transport_K_Realisable_VEF_Face::ajouter(resu);
  //
  // Plutost que de calculer P, on appelle Source_Transport_K_Realisable_VEF_Face::ajouter(resu)
  // et on ajoute directement G
  // On en profite pour faire des tests sur LeK_MIN
  //
  const Zone_VEF& zone_VEF = la_zone_VEF.valeur();
  const Zone_Cl_VEF& zcl_VEF_co = ref_cast(Zone_Cl_VEF,eq_concentration->zone_Cl_dis().valeur());
  const DoubleTab& concen = eq_concentration->inconnue().valeurs();
  const Modele_turbulence_scal_base& le_modele_scalaire =
    ref_cast(Modele_turbulence_scal_base,eq_concentration->get_modele(TURBULENCE).valeur());
  const DoubleTab& alpha_turb = le_modele_scalaire.diffusivite_turbulente().valeurs();
  const DoubleVect& g = gravite->valeurs();
  const Champ_Don& ch_beta_concen = beta_c.valeur();

  const DoubleVect& volumes_entrelaces = zone_VEF.volumes_entrelaces();
  int nb_face = zone_VEF.nb_faces();

  DoubleTrav G(nb_face);

  int nb_consti = eq_concentration->constituant().nb_constituants();

  calculer_terme_destruction_K_gen(zone_VEF,zcl_VEF_co,G,concen,alpha_turb,ch_beta_concen,g,nb_consti);

  for (int face=0; face<nb_face; face++)
    {
      resu(face) += G(face)*volumes_entrelaces(face);
    }
  return resu;
}

DoubleTab& Source_Transport_K_Realisable_aniso_therm_concen_VEF_Face::ajouter(DoubleTab& resu) const
{
  Source_Transport_K_Realisable_VEF_Face::ajouter(resu);
  //
  // Plutost que de calculer P, on appelle Source_Transport_K_Realisable_VEF_Face::ajouter(resu)
  // et on ajoute directement G
  // On en profite pour faire des tests sur LeK_MIN
  //
  const Zone_VEF& zone_VEF = la_zone_VEF.valeur();
  const Zone_Cl_VEF& zcl_VEF_th = ref_cast(Zone_Cl_VEF,eq_thermique->zone_Cl_dis().valeur());
  const Zone_Cl_VEF& zcl_VEF_co = ref_cast(Zone_Cl_VEF,eq_concentration->zone_Cl_dis().valeur());
  const DoubleTab& temper = eq_thermique->inconnue().valeurs();
  const DoubleTab& concen = eq_concentration->inconnue().valeurs();
  const Modele_turbulence_scal_base& le_modele_scalaire =
    ref_cast(Modele_turbulence_scal_base,eq_thermique->get_modele(TURBULENCE).valeur());
  const DoubleTab& alpha_turb = le_modele_scalaire.diffusivite_turbulente().valeurs();
  const DoubleVect& g = gravite->valeurs();
  const Champ_Don& ch_beta_temper = beta_t.valeur();
  const Champ_Don& ch_beta_concen = beta_c.valeur();

  const DoubleVect& volumes_entrelaces = zone_VEF.volumes_entrelaces();
  int nb_face = zone_VEF.nb_faces();

  DoubleTrav P(nb_face);
  DoubleTrav G_t(nb_face);
  DoubleTrav G_c(nb_face);

  int nb_consti = eq_concentration->constituant().nb_constituants();

  calculer_terme_destruction_K_gen(zone_VEF,zcl_VEF_th,G_t,temper,alpha_turb,ch_beta_temper,g,0);
  calculer_terme_destruction_K_gen(zone_VEF,zcl_VEF_co,G_c,concen,alpha_turb,ch_beta_concen,g,nb_consti);

  for (int face=0; face<nb_face; face++)
    {
      resu(face) += ( G_t(face)+G_c(face) ) *volumes_entrelaces(face);
    }
  return resu;
}

DoubleTab& Source_Transport_K_Realisable_anisotherme_VEF_Face::calculer(DoubleTab& resu) const
{
  resu=0;
  return ajouter(resu);
}

DoubleTab& Source_Transport_K_Realisable_aniso_concen_VEF_Face::calculer(DoubleTab& resu) const
{
  resu = 0;
  return ajouter(resu);
}

DoubleTab& Source_Transport_K_Realisable_aniso_therm_concen_VEF_Face::calculer(DoubleTab& resu) const
{
  resu = 0;
  return ajouter(resu);
}

void Source_Transport_K_Realisable_anisotherme_VEF_Face::associer_pb(const Probleme_base& pb)
{
  if (pb.nombre_d_equations()<2)
    {
      Cerr<<"The K_Eps_Realisable source term "<<que_suis_je()<<" cannot be activated"<<finl;
      Cerr<<"for a "<<pb.que_suis_je()<<" problem."<<finl;
    }
  const Equation_base& eqn = pb.equation(1);
  const Milieu_base& milieu = eqn.milieu();
  const Fluide_base& fluide = ref_cast(Fluide_base,milieu);

  if (sub_type(Fluide_Quasi_Compressible,fluide))
    {
      Cerr<<"The K_Eps_Realisable source term "<<que_suis_je()<<" cannot be activated"<<finl;
      Cerr<<"with a "<<milieu.que_suis_je()<<" medium."<<finl;
      exit();
    }
  Source_Transport_K_Realisable_VEF_Face::associer_pb(pb);
  const Convection_Diffusion_Temperature& eqn_th = ref_cast(Convection_Diffusion_Temperature,eqn);
  eq_thermique = eqn_th;
  beta_t=fluide.beta_t();
  gravite = fluide.gravite();
}

void Source_Transport_K_Realisable_aniso_concen_VEF_Face::associer_pb(const Probleme_base& pb)
{
  if (pb.nombre_d_equations()<2)
    {
      Cerr<<"The K_Eps_Realisable source term "<<que_suis_je()<<" cannot be activated"<<finl;
      Cerr<<"for a "<<pb.que_suis_je()<<" problem."<<finl;
    }

  const Equation_base& eqn = pb.equation(1);
  const Milieu_base& milieu = pb.equation(0).milieu();
  const Fluide_base& fluide = ref_cast(Fluide_base,milieu);

  if (sub_type(Fluide_Quasi_Compressible,fluide))
    {
      Cerr<<"The K_Eps_Realisable source term "<<que_suis_je()<<" cannot be activated"<<finl;
      Cerr<<"with a "<<milieu.que_suis_je()<<" medium."<<finl;
      exit();
    }

  Source_Transport_K_Realisable_VEF_Face::associer_pb(pb);

  const Convection_Diffusion_Concentration& eqn_c = ref_cast(Convection_Diffusion_Concentration,eqn);
  eq_concentration = eqn_c;
  if (!fluide.beta_c().non_nul())
    {
      Cerr << "You forgot to define beta_co field in the fluid." << finl;
      Cerr << "It is mandatory when using the K-Eps model (buoyancy effects)." << finl;
      Cerr << "If you don't want buoyancy effects, then specify: beta_co champ_uniforme 1 0." << finl;
      exit();
    }
  beta_c = fluide.beta_c();
  gravite = fluide.gravite();
}

void Source_Transport_K_Realisable_aniso_therm_concen_VEF_Face::associer_pb(const Probleme_base& pb)
{
  if (pb.nombre_d_equations()<3)
    {
      Cerr<<"The K_Eps_Realisable source term "<<que_suis_je()<<" cannot be activated"<<finl;
      Cerr<<"for a "<<pb.que_suis_je()<<" problem."<<finl;
    }
  const Equation_base& eqn_therm = pb.equation(1);
  const Equation_base& eqn_conc = pb.equation(2);
  const Milieu_base& milieu = eqn_therm.milieu();
  const Fluide_base& fluide = ref_cast(Fluide_base,milieu);

  if (sub_type(Fluide_Quasi_Compressible,fluide))
    {
      Cerr<<"The K_Eps_Realisable source term "<<que_suis_je()<<" cannot be activated"<<finl;
      Cerr<<"with a "<<milieu.que_suis_je()<<" medium."<<finl;
      exit();
    }
  Source_Transport_K_Realisable_VEF_Face::associer_pb(pb);

  const Convection_Diffusion_Temperature& eqn_th =
    ref_cast(Convection_Diffusion_Temperature,eqn_therm);
  eq_thermique = eqn_th;
  const Convection_Diffusion_Concentration& eqn_c =
    ref_cast(Convection_Diffusion_Concentration,eqn_conc);
  eq_concentration = eqn_c;
  beta_t = fluide.beta_t();
  if (!fluide.beta_c().non_nul())
    {
      Cerr << "You forgot to define beta_co field in the fluid." << finl;
      Cerr << "It is mandatory when using the K-Eps model (buoyancy effects)." << finl;
      Cerr << "If you don't want buoyancy effects, then specify: beta_co champ_uniforme 1 0." << finl;
      exit();
    }
  beta_c = fluide.beta_c();
  gravite = fluide.gravite();
}

