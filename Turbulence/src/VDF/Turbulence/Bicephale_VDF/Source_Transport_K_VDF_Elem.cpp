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
// File:        Source_Transport_K_VDF_Elem.cpp
// Directory:   $TRUST_ROOT/src/VDF/Turbulence
// Version:     /main/41
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Transport_K_VDF_Elem.h>
#include <Transport_K_ou_Eps.h>
#include <Convection_Diffusion_Temperature.h>
#include <Convection_Diffusion_Concentration.h>
#include <Modele_turbulence_scal_base.h>
#include <Probleme_base.h>
#include <IntTrav.h>
#include <Champ_Uniforme.h>
#include <Zone_VDF.h>
#include <Champ_Face.h>
#include <Zone_Cl_VDF.h>
#include <Fluide_Quasi_Compressible.h>
#include <Debog.h>
#include <Modele_turbulence_hyd_K_Eps_Bicephale.h>
#include <DoubleTrav.h>
#include <Pb_Hydraulique_Turbulent.h>
#include <Pb_Hydraulique_Concentration_Turbulent.h>
#include <Pb_Thermohydraulique_Turbulent_QC.h>
#include <Pb_Thermohydraulique_Turbulent.h>
#include <Pb_Thermohydraulique_Concentration_Turbulent.h>
#include <Param.h>
#include <Constituant.h>

Implemente_instanciable(Source_Transport_K_VDF_Elem,"Source_Transport_K_VDF_P0_VDF",Source_base);
Implemente_instanciable(Source_Transport_K_anisotherme_VDF_Elem,"Source_Transport_K_anisotherme_VDF_P0_VDF",Source_Transport_K_VDF_Elem);
Implemente_instanciable(Source_Transport_K_concen_VDF_Elem,"Source_Transport_K_aniso_concen_VDF_P0_VDF",Source_Transport_K_VDF_Elem);
Implemente_instanciable(Source_Transport_K_aniso_therm_concen_VDF_Elem,"Source_Transport_K_aniso_therm_concen_VDF_P0_VDF",Source_Transport_K_VDF_Elem);

//// printOn
//

Sortie& Source_Transport_K_VDF_Elem::printOn(Sortie& s) const
{
  return s << que_suis_je() ;
}

Sortie& Source_Transport_K_anisotherme_VDF_Elem::printOn(Sortie& s) const
{
  return s << que_suis_je() ;
}

Sortie& Source_Transport_K_concen_VDF_Elem::printOn(Sortie& s) const
{
  return s << que_suis_je() ;
}

Sortie& Source_Transport_K_aniso_therm_concen_VDF_Elem::printOn(Sortie& s) const
{
  return s << que_suis_je() ;
}


//// readOn
//
inline void error(const Nom& source, const Nom& problem)
{
  Cerr << "Error ! You can't use the " << source << " source term for the K equation of the problem: " << problem << finl;
  Cerr << "Check the reference manual. It is may be another source term Source_Transport_K_.... which should be used." << finl;
  Process::exit();
}

Entree& Source_Transport_K_VDF_Elem::readOn(Entree& is)
{
  const Probleme_base& problem = mon_equation->probleme();
  if (!sub_type(Pb_Hydraulique_Turbulent,problem) && !sub_type(Pb_Thermohydraulique_Turbulent_QC,problem)) error(que_suis_je(),problem.que_suis_je());
  Param param(que_suis_je());
  param.lire_avec_accolades(is);
  return is ;
}

Entree& Source_Transport_K_anisotherme_VDF_Elem::readOn(Entree& is)
{
  const Probleme_base& problem = mon_equation->probleme();
  if (!sub_type(Pb_Thermohydraulique_Turbulent,problem)) error(que_suis_je(),problem.que_suis_je());
  Param param(que_suis_je());
  param.lire_avec_accolades(is);
  return is ;
}

Entree& Source_Transport_K_concen_VDF_Elem::readOn(Entree& is)
{
  const Probleme_base& problem = mon_equation->probleme();
  if (!sub_type(Pb_Hydraulique_Concentration_Turbulent,problem)) error(que_suis_je(),problem.que_suis_je());
  Param param(que_suis_je());
  param.lire_avec_accolades(is);
  return is ;
}

Entree& Source_Transport_K_aniso_therm_concen_VDF_Elem::readOn(Entree& is)
{
  const Probleme_base& problem = mon_equation->probleme();
  if (!sub_type(Pb_Thermohydraulique_Concentration_Turbulent,problem)) error(que_suis_je(),problem.que_suis_je());
  Param param(que_suis_je());
  param.lire_avec_accolades(is);
  return is ;
}



/////////////////////////////////////////////////////////////////////////////
//
//            Implementation des fonctions de la classe
//
//             Source_Transport_K_VDF_Elem
//
/////////////////////////////////////////////////////////////////////////////

void Source_Transport_K_VDF_Elem::associer_zones(const Zone_dis& zone_dis,
                                                 const Zone_Cl_dis&  )
{
  la_zone_VDF = ref_cast(Zone_VDF, zone_dis.valeur());
}

// remplit les references eq_hydraulique et mon_eq_transport_K_Eps
void Source_Transport_K_VDF_Elem::associer_pb(const Probleme_base& pb)
{
  eq_hydraulique = pb.equation(0);
  mon_eq_transport_K   = ref_cast(Transport_K_ou_Eps,equation());
  mon_eq_transport_Eps = ref_cast(Transport_K_ou_Eps, mon_eq_transport_K.valeur().modele_turbulence().eqn_transp_Eps());
}

DoubleTab& Source_Transport_K_VDF_Elem::ajouter(DoubleTab& resu) const
{

  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const Zone_Cl_VDF& zcl_VDF = ref_cast(Zone_Cl_VDF,eq_hydraulique->zone_Cl_dis().valeur());
  const DoubleTab& K   = mon_eq_transport_K->inconnue().valeurs();
  const DoubleTab& Eps = mon_eq_transport_Eps->inconnue().valeurs();
  const DoubleTab& visco_turb = mon_eq_transport_K->modele_turbulence().viscosite_turbulente().valeurs();
  const DoubleTab& vit = eq_hydraulique->inconnue().valeurs();
  const DoubleVect& volumes = zone_VDF.volumes();
  const DoubleVect& porosite_vol = zone_VDF.porosite_elem();

  int nb_elem = zone_VDF.nb_elem();

  // Ajout d'un espace virtuel au tableau P
  DoubleVect P;

  zone_VDF.zone().creer_tableau_elements(P);

  if (axi)
    {
      const Champ_Face& vitesse = ref_cast(Champ_Face,eq_hydraulique->inconnue().valeur());
      calculer_terme_production_K_BiK_Axi(zone_VDF,vitesse,P,K,visco_turb);
    }
  else
    {
      const Champ_Face& vitesse = ref_cast(Champ_Face,eq_hydraulique->inconnue().valeur());
      calculer_terme_production_K_BiK(zone_VDF,zcl_VDF,P,K,vit,vitesse,visco_turb);
    }
  const Modele_Fonc_Bas_Reynolds& mon_modele_fonc=ref_cast(Modele_turbulence_hyd_K_Eps_Bicephale,mon_eq_transport_K->modele_turbulence()).associe_modele_fonction();
  int is_modele_fonc=(mon_modele_fonc.non_nul());
  //is_modele_fonc=0;
  if (is_modele_fonc)
    {

      DoubleTab& D=ref_cast_non_const(DoubleTab,mon_modele_fonc.valeur().get_champ("D").valeurs());
      DoubleTab& E=ref_cast_non_const(DoubleTab,mon_modele_fonc.valeur().get_champ("E").valeurs());
      DoubleTab& F1=ref_cast_non_const(DoubleTab,mon_modele_fonc.valeur().get_champ("F1").valeurs());
      DoubleTab& F2=ref_cast_non_const(DoubleTab,mon_modele_fonc.valeur().get_champ("F2").valeurs());
      const Fluide_Incompressible& fluide=ref_cast(Fluide_Incompressible,eq_hydraulique.valeur().milieu());
      const Champ_Don& ch_visco_cin = fluide.viscosite_cinematique();
      // const DoubleTab& tab_visco = ch_visco_cin->valeurs();
      const Zone_Cl_dis& zcl_k=mon_eq_transport_K->zone_Cl_dis();
      const Zone_dis& zone_dis_k =mon_eq_transport_K->zone_dis();
      mon_modele_fonc.Calcul_D_BiK(D,zone_dis_k,zcl_k,vit,K, Eps,ch_visco_cin);
      mon_modele_fonc.Calcul_E_BiK(E,zone_dis_k,zcl_k,vit,K, Eps,ch_visco_cin,visco_turb);

      D.echange_espace_virtuel();
      E.echange_espace_virtuel();
      const Champ_base& ch_visco_cin_ou_dyn =ref_cast(Op_Diff_K_Eps_base, equation().operateur(0).l_op_base()).diffusivite();
      DoubleTab P_tab;
      P_tab.ref(P);
      mon_modele_fonc.Calcul_F1_BiK(F1,zone_dis_k,zcl_k, P_tab, K,Eps,ch_visco_cin_ou_dyn);


      mon_modele_fonc.Calcul_F2_BiK(F2,D,zone_dis_k,K, Eps, ch_visco_cin_ou_dyn  );
      Debog::verifier("D",D);
      Debog::verifier("E",E);
      Debog::verifier("F2",F2);
      Debog::verifier("F1",F1);
      Debog::verifier("avt",resu);
      for (int elem=0; elem<nb_elem; elem++)
        {
          //if (K(elem) > 1.e-10 && Eps(elem) > 1.e-10)
          {
            // resu(elem) += (P(elem)-eps_sur_k(elem)*K(elem))*volumes(elem)*porosite_vol(elem)-D(elem);
            resu(elem) += (P(elem)-Eps(elem)-D(elem))*volumes(elem)*porosite_vol(elem);
            // on en profite pour diviser D  par le volume
          }
        }
      Debog::verifier("ap",resu);

      // int elem=0;
      //Cerr<<" ici "<<D(elem)<< " "<<Eps(elem)<<finl;
    }
  else
    {
      for (int elem=0; elem<nb_elem; elem++)
        {
          resu(elem) += (P(elem)-Eps(elem))*volumes(elem)*porosite_vol(elem);
        }
    }

  //Debog::verifier("Source_Transport_K_VDF_Elem::ajouter resu",resu);
  resu.echange_espace_virtuel();
  return resu;
}

DoubleTab& Source_Transport_K_VDF_Elem::calculer(DoubleTab& resu) const
{
  resu =0.;
  return ajouter(resu);
}

/////////////////////////////////////////////////////////////////////////////
//
//            Implementation des fonctions de la classe
//
//             Source_Transport_K_anisotherme_VDF_Elem
//
/////////////////////////////////////////////////////////////////////////////
//

// remplit les references eq_thermique et beta
void Source_Transport_K_anisotherme_VDF_Elem::associer_pb(const Probleme_base& pb)
{
  if (pb.nombre_d_equations()<2)
    {
      Cerr<<"The K source term "<<que_suis_je()<<" cannot be activated"<<finl;
      Cerr<<"for a "<<pb.que_suis_je()<<" problem."<<finl;
    }
  const Equation_base& eqn = pb.equation(1);
  const Milieu_base& milieu = eqn.milieu();
  const Fluide_Incompressible& fluide = ref_cast(Fluide_Incompressible,milieu);

  if (sub_type(Fluide_Quasi_Compressible,fluide))
    {
      Cerr<<"The K source term "<<que_suis_je()<<" cannot be activated"<<finl;
      Cerr<<"with a "<<milieu.que_suis_je()<<" medium."<<finl;
      exit();
    }
  Source_Transport_K_VDF_Elem::associer_pb(pb);

  const Convection_Diffusion_Temperature& eqn_th =
    ref_cast(Convection_Diffusion_Temperature,eqn);
  eq_thermique = eqn_th;

  beta_t=fluide.beta_t();
  gravite = fluide.gravite();
}

DoubleTab& Source_Transport_K_anisotherme_VDF_Elem::ajouter(DoubleTab& resu) const
{
  Source_Transport_K_VDF_Elem::ajouter(resu);

  // Modifs VB : plutot que de calculer P, on appelle Source_Transport_K_VDF_Elem::ajouter(resu)
  // et on ajoute directement G
  //

  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const Zone_Cl_VDF& zcl_VDF_th = ref_cast(Zone_Cl_VDF,eq_thermique->zone_Cl_dis().valeur());
  const DoubleTab& scalaire = eq_thermique->inconnue().valeurs();
  const Modele_turbulence_scal_base& le_modele_scalaire =
    ref_cast(Modele_turbulence_scal_base,eq_thermique->get_modele(TURBULENCE).valeur());
  DoubleTab alpha_turb(le_modele_scalaire.diffusivite_turbulente().valeurs());
  const DoubleTab& g = gravite->valeurs();
  const Champ_Don& ch_beta = beta_t.valeur();
  const DoubleVect& volumes = zone_VDF.volumes();
  const DoubleVect& porosite_vol = zone_VDF.porosite_elem();

  // Ajout d'un espace virtuel au tableau G
  DoubleVect G;
  zone_VDF.zone().creer_tableau_elements(G);

  //else
  //  calculer_terme_production_K(zone_VDF,P,K_eps,vit,visco_turb);

  // C'est l'objet de type zone_Cl_dis de l'equation thermique
  // qui est utilise dans le calcul de G

  const DoubleTab& tab_beta = ch_beta.valeurs();
  if (sub_type(Champ_Uniforme,ch_beta.valeur()))
    calculer_terme_destruction_K(zone_VDF,zcl_VDF_th,G,scalaire,alpha_turb,tab_beta(0,0),g);
  else
    calculer_terme_destruction_K(zone_VDF,zcl_VDF_th,G,scalaire,alpha_turb,tab_beta,g);

  int nb_elem = zone_VDF.nb_elem();
  for (int elem=0; elem<nb_elem; elem++)
    {
      resu(elem) += G(elem)*volumes(elem)*porosite_vol(elem);
    }

  //Debog::verifier("Source_Transport_K_VDF_Elem::ajouter : " ,resu);
  return resu;
}

DoubleTab& Source_Transport_K_anisotherme_VDF_Elem::calculer(DoubleTab& resu) const
{
  resu=0;
  return ajouter(resu);
}

/////////////////////////////////////////////////////////////////////////////
//
//            Implementation des fonctions de la classe
//
//             Source_Transport_K_concen_VDF_Elem
//
/////////////////////////////////////////////////////////////////////////////
//

// remplit les references
void Source_Transport_K_concen_VDF_Elem::associer_pb(const Probleme_base& pb)
{
  if (pb.nombre_d_equations()<2)
    {
      Cerr<<"The K source term "<<que_suis_je()<<" cannot be activated"<<finl;
      Cerr<<"for a "<<pb.que_suis_je()<<" problem."<<finl;
    }

  const Equation_base& eqn = pb.equation(1);
  const Milieu_base& milieu = pb.equation(0).milieu();
  const Fluide_Incompressible& fluide = ref_cast(Fluide_Incompressible,milieu);

  if (sub_type(Fluide_Quasi_Compressible,fluide))
    {
      Cerr<<"The K source term "<<que_suis_je()<<" cannot be activated"<<finl;
      Cerr<<"with a "<<milieu.que_suis_je()<<" medium."<<finl;
      exit();
    }
  Source_Transport_K_VDF_Elem::associer_pb(pb);

  const Convection_Diffusion_Concentration& eqn_c =
    ref_cast(Convection_Diffusion_Concentration,eqn);
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

DoubleTab& Source_Transport_K_concen_VDF_Elem::ajouter(DoubleTab& resu) const
{
  Source_Transport_K_VDF_Elem::ajouter(resu);
  //
  //// Modifs VB : plutot que de calculer P, on appelle Source_Transport_K_VDF_Elem::ajouter(resu)
  //// et on ajoute directement G
  ////
  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const Zone_Cl_VDF& zcl_VDF_co = ref_cast(Zone_Cl_VDF,eq_concentration->zone_Cl_dis().valeur());
  const DoubleTab& concen = eq_concentration->inconnue().valeurs();
  const Modele_turbulence_scal_base& le_modele_scalaire =
    ref_cast(Modele_turbulence_scal_base,eq_concentration->get_modele(TURBULENCE).valeur());
  const DoubleTab& diffu_turb  = le_modele_scalaire.conductivite_turbulente().valeurs();
  const Champ_Uniforme& ch_beta_concen = ref_cast(Champ_Uniforme, beta_c->valeur());
  const DoubleVect& g = gravite->valeurs();
  const DoubleVect& volumes = zone_VDF.volumes();
  const DoubleVect& porosite_vol = zone_VDF.porosite_elem();
  int nb_elem = zone_VDF.nb_elem();
  int nb_consti = eq_concentration->constituant().nb_constituants();

  // Ajout d'un espace virtuel au tableau G
  DoubleVect G;
  zone_VDF.zone().creer_tableau_elements(G);

  if (nb_consti == 1)
    {
      double d_beta_c = ch_beta_concen(0,0);
      calculer_terme_destruction_K(zone_VDF,zcl_VDF_co,G,
                                   concen,diffu_turb,d_beta_c,g);
    }
  else
    {
      const DoubleVect& d_beta_c = ch_beta_concen.valeurs();
      calculer_terme_destruction_K(zone_VDF,zcl_VDF_co,G,
                                   concen,diffu_turb,d_beta_c,g,
                                   nb_consti);
    }

  for (int elem=0; elem<nb_elem; elem++)
    {
      resu(elem) += G(elem)*volumes(elem)*porosite_vol(elem);
    }
  return resu;
}

DoubleTab& Source_Transport_K_concen_VDF_Elem::calculer(DoubleTab& resu) const
{
  resu=0;
  return ajouter(resu);
}

/////////////////////////////////////////////////////////////////////////////
//
//            Implementation des fonctions de la classe
//
//             Source_Transport_K_aniso_therm_concen_VDF_Elem
//
/////////////////////////////////////////////////////////////////////////////

// remplit les references eq_thermique, eq_concentration, beta_t_, beta_c_
void Source_Transport_K_aniso_therm_concen_VDF_Elem::associer_pb(const Probleme_base& pb)
{
  if (pb.nombre_d_equations()<3)
    {
      Cerr<<"The K source term "<<que_suis_je()<<" cannot be activated"<<finl;
      Cerr<<"for a "<<pb.que_suis_je()<<" problem."<<finl;
    }
  const Equation_base& eqn_therm = pb.equation(1);
  const Equation_base& eqn_conc = pb.equation(2);
  const Milieu_base& milieu = eqn_therm.milieu();
  const Fluide_Incompressible& fluide = ref_cast(Fluide_Incompressible,milieu);

  if (sub_type(Fluide_Quasi_Compressible,fluide))
    {
      Cerr<<"The K source term "<<que_suis_je()<<" cannot be activated"<<finl;
      Cerr<<"with a "<<milieu.que_suis_je()<<" medium."<<finl;
      exit();
    }
  Source_Transport_K_VDF_Elem::associer_pb(pb);

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

DoubleTab& Source_Transport_K_aniso_therm_concen_VDF_Elem::ajouter(DoubleTab& resu) const
{

  Source_Transport_K_VDF_Elem::ajouter(resu);
  //
  //
  //// Modifs VB : plutot que de calculer P, on appelle Source_Transport_K_VDF_Elem::ajouter(resu)
  //// et on ajoute directement G
  ////
  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const Zone_Cl_VDF& zcl_VDF_th = ref_cast(Zone_Cl_VDF,eq_thermique->zone_Cl_dis().valeur());
  const Zone_Cl_VDF& zcl_VDF_co = ref_cast(Zone_Cl_VDF,eq_concentration->zone_Cl_dis().valeur());
  const DoubleTab& temper = eq_thermique->inconnue().valeurs();
  const DoubleTab& concen = eq_concentration->inconnue().valeurs();
  const Modele_turbulence_scal_base& le_modele_scalaire =
    ref_cast(Modele_turbulence_scal_base,eq_thermique->get_modele(TURBULENCE).valeur());
  DoubleTab alpha_turb(le_modele_scalaire.conductivite_turbulente().valeurs());
  double rhocp = eq_thermique->milieu().capacite_calorifique().valeurs()(0, 0) * eq_thermique->milieu().masse_volumique().valeurs()(0, 0);
  alpha_turb /= rhocp;
  const Modele_turbulence_scal_base& le_modele_scal_co =
    ref_cast(Modele_turbulence_scal_base,eq_concentration->get_modele(TURBULENCE).valeur());
  const DoubleTab& diffu_turb = le_modele_scal_co.conductivite_turbulente().valeurs();
  const DoubleVect& g = gravite->valeurs();
  const Champ_Don& ch_beta_temper = beta_t.valeur();
  if (!beta_c->non_nul())
    {
      Cerr << finl << "Le champ beta_co n'a pas ete renseigne pour le fluide dans le jeu de donnees." << finl;
      exit();
    }
  const Champ_Uniforme& ch_beta_concen = ref_cast(Champ_Uniforme, beta_c->valeur());

  const DoubleVect& volumes = zone_VDF.volumes();
  const DoubleVect& porosite_vol = zone_VDF.porosite_elem();
  int nb_elem = zone_VDF.nb_elem();
  int nb_consti = eq_concentration->constituant().nb_constituants();

  // Ajout d'un espace virtuel au tableaux Gt et Gc
  DoubleVect G_t;
  DoubleVect G_c;
  zone_VDF.zone().creer_tableau_elements(G_t);
  zone_VDF.zone().creer_tableau_elements(G_c);
  //Cerr<<"!!!!!>>>!!!! nb elem de G_t :"<<G_t.dimension(0);
  //Cerr<<"!!!!!>>>!!!! nb elem de P :"<<P.dimension(0);
  //  DoubleTrav P(nb_elem_tot);
  //   DoubleTrav G_t(nb_elem_tot);
  //   DoubleTrav G_c(nb_elem_tot);

  //else
  //calculer_terme_production_K(zone_VDF,P,K_eps,vit,visco_turb);

  const DoubleTab& tab_beta_t = ch_beta_temper.valeurs();
  if (sub_type(Champ_Uniforme,ch_beta_temper.valeur()))
    calculer_terme_destruction_K(zone_VDF,zcl_VDF_th,G_t,temper,alpha_turb,tab_beta_t(0,0),g);
  else
    calculer_terme_destruction_K(zone_VDF,zcl_VDF_th,G_t,temper,alpha_turb,tab_beta_t,g);

  if (nb_consti == 1)
    {
      double d_beta_c = ch_beta_concen(0,0);
      calculer_terme_destruction_K(zone_VDF,zcl_VDF_co,G_c,
                                   concen,diffu_turb,d_beta_c,g);
    }
  else
    {
      const DoubleVect& d_beta_c = ch_beta_concen.valeurs();
      calculer_terme_destruction_K(zone_VDF,zcl_VDF_co,G_c,
                                   concen,diffu_turb,d_beta_c,g,
                                   nb_consti);
    }

  double G_sum;

  for (int elem=0; elem<nb_elem; elem++)
    {
      G_sum = G_t(elem)+G_c(elem) ;

      resu(elem) += G_sum *volumes(elem)*porosite_vol(elem);
    }
  return resu;
}

DoubleTab& Source_Transport_K_aniso_therm_concen_VDF_Elem::calculer(DoubleTab& resu) const
{
  resu = 0;
  return ajouter(resu);
}


void  Source_Transport_K_VDF_Elem::contribuer_a_avec(const DoubleTab& a,  Matrice_Morse& matrice) const
{
  const DoubleTab& K   = mon_eq_transport_K->inconnue().valeurs();
  const DoubleTab& Eps = mon_eq_transport_Eps->inconnue().valeurs();

  int size=K.dimension(0);

  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const DoubleVect& porosite = zone_VDF.porosite_elem();
  // on implicite le -eps et le -eps^2/k

  const Modele_Fonc_Bas_Reynolds& mon_modele_fonc=ref_cast(Modele_turbulence_hyd_K_Eps_Bicephale,mon_eq_transport_K->modele_turbulence()).associe_modele_fonction();
  int is_modele_fonc=(mon_modele_fonc.non_nul());
  //is_modele_fonc=0;
  DoubleTab F2;
  if (is_modele_fonc)
    {
      DoubleTrav D(0);
      F2.resize(K.dimension_tot(0));
      const Zone_dis& zone_dis_k =mon_eq_transport_K->zone_dis();

      const Champ_base& ch_visco_cin_ou_dyn =ref_cast(Op_Diff_K_Eps_base, equation().operateur(0).l_op_base()).diffusivite();

      mon_modele_fonc.Calcul_F2_BiK(F2,D,zone_dis_k,K,Eps, ch_visco_cin_ou_dyn  );
    }

  {
    const DoubleVect& volumes=zone_VDF.volumes();
    for (int c=0; c<size; c++)
      {
        // -eps*vol  donne +vol dans la bonne case
        if (K(c)>DMINFLOAT)
          {
            double coef_k=porosite(c)*volumes(c)*Eps(c)/K(c);
            matrice(c,c)+=coef_k;
          }
      }
  }
}



void Source_Transport_K_VDF_Elem::mettre_a_jour(double temps)
{

  Calcul_Production_K_VDF::mettre_a_jour(temps);

}
