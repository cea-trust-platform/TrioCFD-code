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
// File:        Source_Transport_K_Eps_Realisable_VDF_Elem.cpp
// Directory:   $TRUST_ROOT/src/VDF/Turbulence
// Version:     /main/23
//
//////////////////////////////////////////////////////////////////////////////


#include <Source_Transport_K_Eps_Realisable_VDF_Elem.h>
#include <Convection_Diffusion_Temperature.h>
#include <Convection_Diffusion_Concentration.h>
#include <Modele_turbulence_scal_base.h>
#include <Fluide_Incompressible.h>
#include <Probleme_base.h>
#include <Champ_Uniforme.h>
#include <Champ_Face.h>
#include <Zone_VDF.h>
#include <Zone_Cl_VDF.h>
#include <Modele_turbulence_hyd_K_Eps_Realisable.h>
#include <DoubleTrav.h>
#include <Fluide_Quasi_Compressible.h>
#include <Convection_Diffusion_Temperature_Turbulent.h>
#include <Ref_Transport_K_Eps_Realisable.h>
#include <Constituant.h>

Implemente_instanciable_sans_constructeur(Source_Transport_K_Eps_Realisable_VDF_Elem,"Source_Transport_K_Eps_Realisable_VDF_P0_VDF",Source_base);

Implemente_instanciable_sans_constructeur(Source_Transport_K_Eps_Realisable_anisotherme_VDF_Elem,"Source_Transport_K_Eps_Realisable_anisotherme_VDF_P0_VDF",Source_Transport_K_Eps_Realisable_VDF_Elem);
Implemente_instanciable_sans_constructeur(Source_Transport_K_Eps_Realisable_aniso_concen_VDF_Elem,"Source_Transport_K_Eps_Realisable_aniso_concen_VDF_P0_VDF",Source_Transport_K_Eps_Realisable_anisotherme_VDF_Elem);

Implemente_instanciable_sans_constructeur(Source_Transport_K_Eps_Realisable_aniso_therm_concen_VDF_Elem,"Source_Transport_K_Eps_Realisable_aniso_therm_concen_VDF_P0_VDF",Source_Transport_K_Eps_Realisable_VDF_Elem);


//// printOn
//

Sortie& Source_Transport_K_Eps_Realisable_VDF_Elem::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}

Sortie& Source_Transport_K_Eps_Realisable_anisotherme_VDF_Elem::printOn(Sortie& s) const
{
  return s << que_suis_je() ;
}
Sortie& Source_Transport_K_Eps_Realisable_aniso_concen_VDF_Elem::printOn(Sortie& s) const
{
  return s << que_suis_je() ;
}

Sortie& Source_Transport_K_Eps_Realisable_aniso_therm_concen_VDF_Elem::printOn(Sortie& s) const
{
  return s << que_suis_je() ;
}

//// readOn
//

Entree& Source_Transport_K_Eps_Realisable_VDF_Elem::readOn(Entree& is )
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

Entree& Source_Transport_K_Eps_Realisable_anisotherme_VDF_Elem::readOn(Entree& is)
{
  Motcle accolade_ouverte("{");
  Motcle accolade_fermee("}");
  Motcle motlu;

  is >> motlu;
  if (motlu != accolade_ouverte)
    {
      Cerr << "On attendait { pour commencer a lire les constantes de Source_Transport_K_Eps_Realisable_anisotherme" << finl;
      exit();
    }
  Cerr << "Lecture des constantes de Source_Transport_K_Eps_Realisable_anisotherme" << finl;
  Motcles les_mots(2);
  {
    les_mots[0] = "C2_eps";
    les_mots[1] = "C3_eps";
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
        case 1 :
          {
            is >> C3_;
            break;
          }
        default :
          {
            Cerr << "On ne comprend pas le mot cle : " << motlu << "dans Source_Transport_K_Eps_Realisable_anisotherme" << finl;
            exit();
          }
        }

      is >> motlu;
    }
  return is;
}
Entree& Source_Transport_K_Eps_Realisable_aniso_concen_VDF_Elem::readOn(Entree& is)
{
  Motcle accolade_ouverte("{");
  Motcle accolade_fermee("}");
  Motcle motlu;

  is >> motlu;
  if (motlu != accolade_ouverte)
    {
      Cerr << "On attendait { pour commencer a lire les constantes de Source_Transport_K_Eps_Realisable_aniso_concen (peut-estre choisi automatiquement en raison de la nature du probleme)" << finl;
      exit();
    }
  Cerr << "Lecture des constantes de Source_Transport_K_Eps_Realisable_aniso_concen (peut-estre choisi automatiquement en raison de la nature du probleme)" << finl;
  Motcles les_mots(2);
  {
    les_mots[0] = "C2_eps";
    les_mots[1] = "C3_eps";
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
        case 1 :
          {
            is >> C3_;
            break;
          }
        default :
          {
            Cerr << "On ne comprend pas le mot cle : " << motlu << "dans Source_Transport_K_Eps_Realisable_aniso_concen (peut-estre choisi automatiquement en raison de la nature du probleme)" << finl;
            exit();
          }
        }

      is >> motlu;
    }
  return is;
}

Entree& Source_Transport_K_Eps_Realisable_aniso_therm_concen_VDF_Elem::readOn(Entree& is)
{
  Motcle accolade_ouverte("{");
  Motcle accolade_fermee("}");
  Motcle motlu;

  is >> motlu;
  if (motlu != accolade_ouverte)
    {
      Cerr << "On attendait { pour commencer a lire les constantes de Source_Transport_K_Eps_Realisable_aniso_therm_concen (peut-estre choisi automatiquement en raison de la nature du probleme)" << finl;
      exit();
    }
  Cerr << "Lecture des constantes de Source_Transport_K_Eps_Realisable_aniso_therm_concen (peut-estre choisi automatiquement en raison de la nature du probleme)" << finl;
  Motcles les_mots(2);
  {
    les_mots[0] = "C2_eps";
    les_mots[1] = "C3_eps";
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
        case 1 :
          {
            is >> C3_;
            break;
          }
        default :
          {
            Cerr << "On ne comprend pas le mot cle : " << motlu << "dans Source_Transport_K_Eps_Realisable_aniso_therm_concen (peut-estre choisi automatiquement en raison de la nature du probleme)" << finl;
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
//             Source_Transport_K_Eps_Realisable_VDF_Elem
//
/////////////////////////////////////////////////////////////////////////////

void Source_Transport_K_Eps_Realisable_VDF_Elem::associer_pb(const Probleme_base& pb )
{
  eq_hydraulique = pb.equation(0);
  eqn_keps_Rea = ref_cast(Transport_K_Eps_Realisable,equation());
}

void Source_Transport_K_Eps_Realisable_VDF_Elem::associer_zones(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis)
{
  la_zone_VDF = ref_cast(Zone_VDF, zone_dis.valeur());
  la_zone_Cl_VDF = ref_cast(Zone_Cl_VDF, zone_Cl_dis.valeur());
}

DoubleTab& Source_Transport_K_Eps_Realisable_VDF_Elem::ajouter(DoubleTab& resu) const
{
  const Zone_VDF&                               zone_VDF = ref_cast(Zone_VDF,eq_hydraulique->zone_dis().valeur());
  const Zone_Cl_VDF&                             zcl_VDF = ref_cast(Zone_Cl_VDF,eq_hydraulique->zone_Cl_dis().valeur());
  const DoubleTab&                             K_eps_Rea = eqn_keps_Rea->inconnue().valeurs();
  const Modele_turbulence_hyd_K_Eps_Realisable& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_Realisable,eqn_keps_Rea->modele_turbulence());
  const DoubleTab&                            visco_turb = mod_turb.viscosite_turbulente().valeurs();
  const Modele_Fonc_Realisable_base&     mon_modele_fonc = mod_turb.associe_modele_fonction();
  const Fluide_Incompressible&                    fluide = ref_cast(Fluide_Incompressible,eq_hydraulique->milieu());
  const DoubleTab&                                  vit  = eq_hydraulique->inconnue().valeurs();

  const Champ_Don& ch_visco_cin = fluide.viscosite_cinematique();
  const DoubleTab& tab_visco    = ch_visco_cin->valeurs();

  int is_visco_const=sub_type(Champ_Uniforme,ch_visco_cin.valeur());

  double visco=-1;
  if (is_visco_const)
    {
      visco = max(tab_visco(0,0),DMINFLOAT);
    }

  const DoubleVect& volumes = zone_VDF.volumes();
  const DoubleVect& porosite_vol = zone_VDF.porosite_elem();
  int nb_elem = zone_VDF.nb_elem();
  //  int nb_elem_tot = zone_VDF.nb_elem_tot();
  DoubleTrav P(visco_turb);
  DoubleTrav C1(nb_elem);
  DoubleTrav S(nb_elem);

  C1 = mon_modele_fonc.get_C1( );
  S  = mon_modele_fonc.get_S( );

  double   LeK_MIN = mod_turb.get_LeK_MIN();
  double LeEPS_MIN = mod_turb.get_LeEPS_MIN();

  // Rq : la distinction entre zone_cl est importante pour les deux equations pour l imposition des conditions aux limites!!!!!!
  if (axi)
    {
      Champ_Face& vitesse = ref_cast_non_const(Champ_Face,eq_hydraulique->inconnue().valeur());

      calculer_terme_production_K_Axi(zone_VDF,vitesse,P,K_eps_Rea,visco_turb);
    }
  else
    {
      Champ_Face& vitesse = ref_cast_non_const(Champ_Face,eq_hydraulique->inconnue().valeur());
      calculer_terme_production_K(zone_VDF,zcl_VDF,P,K_eps_Rea,vit,vitesse,visco_turb);
    }

  P.echange_espace_virtuel();

  for (int elem=0; elem<nb_elem; elem++)
    {
      if (!is_visco_const)
        {
          visco =  tab_visco(elem);
        }

      assert(visco>0.);

      resu(elem,0) += ( P(elem)-K_eps_Rea(elem,1) )*volumes(elem)*porosite_vol(elem);

      if ( ( K_eps_Rea(elem,0) >= LeK_MIN ) and ( K_eps_Rea(elem,1) >= LeEPS_MIN ) )
        {
          resu(elem,1) += K_eps_Rea(elem,1)*( C1(elem)*S(elem) - ( C2_*K_eps_Rea(elem,1)/( K_eps_Rea(elem,0) + sqrt( visco*K_eps_Rea(elem,1) ) ) ) )*volumes(elem)*porosite_vol(elem);
        }
    }

  return resu;
}

DoubleTab& Source_Transport_K_Eps_Realisable_VDF_Elem::calculer(DoubleTab& resu) const
{
  resu = 0.;
  return ajouter(resu);
}

void Source_Transport_K_Eps_Realisable_VDF_Elem::mettre_a_jour(double temps)
{
  const Zone_Cl_dis&   zcl_keps  = eqn_keps_Rea->zone_Cl_dis();
  const Zone_dis& zone_dis_keps  = eqn_keps_Rea ->zone_dis();
  const DoubleTab&    K_eps_Rea  = eqn_keps_Rea->inconnue().valeurs();
  const DoubleTab&          vit  = eq_hydraulique->inconnue().valeurs();

  Modele_turbulence_hyd_K_Eps_Realisable&       mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_Realisable,eqn_keps_Rea->modele_turbulence());
  Modele_Fonc_Realisable_base&           mon_modele_fonc = mod_turb.associe_modele_fonction();
  double                                 epsilon_minimum = eqn_keps_Rea.valeur().modele_turbulence().get_LeEPS_MIN();

  mon_modele_fonc.Contributions_Sources(zone_dis_keps,zcl_keps,vit,K_eps_Rea,epsilon_minimum);

  Calcul_Production_K_VDF::mettre_a_jour(temps);
}

void  Source_Transport_K_Eps_Realisable_VDF_Elem::contribuer_a_avec(const DoubleTab& a,  Matrice_Morse& matrice) const
{
  const DoubleTab& val=equation().inconnue().valeurs();
  int size=val.dimension(0);

  const Zone_VDF&   zone_VDF = la_zone_VDF.valeur();
  const DoubleVect& porosite = zone_VDF.porosite_elem();

  const Fluide_Incompressible&  fluide = ref_cast(Fluide_Incompressible,eq_hydraulique->milieu());
  const Champ_Don&        ch_visco_cin = fluide.viscosite_cinematique();
  const DoubleTab&           tab_visco = ch_visco_cin->valeurs();

  int is_visco_const=sub_type(Champ_Uniforme,ch_visco_cin.valeur());

  double visco=-1;
  if (is_visco_const)
    {
      visco = max(tab_visco(0,0),DMINFLOAT);
    }

  // on implicite le -eps et le -eps^2/k

  const DoubleVect& volumes=zone_VDF.volumes();
  for (int c=0; c<size; c++)
    {
      if (!is_visco_const)
        {
          visco =  tab_visco(c);
        }

      assert(visco>0.);

      // -eps*vol  donne +vol dans la bonne case
      if (val(c,0)>DMINFLOAT)

        {

          double coef_k=porosite(c)*volumes(c)*val(c,1)/( val(c,0) + sqrt( visco*val(c,1) ) );
          matrice(c*2,c*2)+=coef_k;

          double coef_eps=C2_*coef_k;

          matrice(c*2+1,c*2+1)+=coef_eps;
        }
    }

}

/////////////////////////////////////////////////////////////////////////////
//
//            Implementation des fonctions de la classe
//
//             Source_Transport_K_Eps_Realisable_anisotherme_VDF_Elem
//
/////////////////////////////////////////////////////////////////////////////
//

// remplit les references eq_thermique et beta

void Source_Transport_K_Eps_Realisable_anisotherme_VDF_Elem::associer_pb(const Probleme_base& pb)
{
  if (pb.nombre_d_equations()<2)
    {
      Cerr<<"The K_Eps_Realisable source term "<<que_suis_je()<<" cannot be activated"<<finl;
      Cerr<<"for a "<<pb.que_suis_je()<<" problem."<<finl;
    }
  const Equation_base& eqn = pb.equation(1);
  const Milieu_base& milieu = eqn.milieu();
  const Fluide_Incompressible& fluide = ref_cast(Fluide_Incompressible,milieu);

  if (sub_type(Fluide_Quasi_Compressible,fluide))
    {
      Cerr<<"The K_Eps_Realisable source term "<<que_suis_je()<<" cannot be activated"<<finl;
      Cerr<<"with a "<<milieu.que_suis_je()<<" medium."<<finl;
      exit();
    }
  Source_Transport_K_Eps_Realisable_VDF_Elem::associer_pb(pb);

  const Convection_Diffusion_Temperature& eqn_th =
    ref_cast(Convection_Diffusion_Temperature,eqn);
  eq_thermique = eqn_th;

  beta_t=fluide.beta_t();
  gravite = fluide.gravite();
}

DoubleTab& Source_Transport_K_Eps_Realisable_anisotherme_VDF_Elem::ajouter(DoubleTab& resu) const
{
  Source_Transport_K_Eps_Realisable_VDF_Elem::ajouter(resu);
  //
  //// Plutost que de calculer P, on appelle Source_Transport_K_Eps_Realisable_VDF_Elem::ajouter(resu)
  //// et on ajoute directement G
  ////
  //
  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const Zone_Cl_VDF& zcl_VDF_th = ref_cast(Zone_Cl_VDF,eq_thermique->zone_Cl_dis().valeur());
  const DoubleTab& K_eps = eqn_keps_Rea->inconnue().valeurs();
  const DoubleTab& scalaire = eq_thermique->inconnue().valeurs();
  const Modele_turbulence_scal_base& le_modele_scalaire =
    ref_cast(Modele_turbulence_scal_base,eq_thermique->get_modele(TURBULENCE).valeur());
  const DoubleTab& alpha_turb =        le_modele_scalaire.diffusivite_turbulente().valeurs();
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

  double C1_loc=1.44; // C1 value is not a constant in Realizable K-Epsilon model but here, we take the default value of C1 used in standard K-Epsilon, as proposed by litterature
  double C3_loc ;
  //const Mod_turb_hyd_RANS& mod_turb_RANS = ref_cast(Mod_turb_hyd_RANS,eq_hydraulique->modele_turbulence().valeur());
  //double LeK_MIN = mod_turb_RANS.get_LeK_MIN() ;
  double LeK_MIN = eqn_keps_Rea->modele_turbulence().get_LeK_MIN() ;
  int nb_elem = zone_VDF.nb_elem();
  for (int elem=0; elem<nb_elem; elem++)
    {
      resu(elem,0) += G(elem)*volumes(elem)*porosite_vol(elem);

      if (K_eps(elem,0) >= LeK_MIN)
        {
          C3_loc = C3_ ;
          if ( G(elem) > 0. ) C3_loc = 0. ;

          resu(elem,1) += C1_loc*(1-C3_loc)*G(elem)*volumes(elem)*porosite_vol(elem)
                          * K_eps(elem,1)/K_eps(elem,0);
        }
    }

  //Debog::verifier("Source_Transport_K_Eps_Realisable_anisotherme_VDF_Elem::ajouter : " ,resu);
  return resu;
}

DoubleTab& Source_Transport_K_Eps_Realisable_anisotherme_VDF_Elem::calculer(DoubleTab& resu) const
{

  resu = 0.;
  return ajouter(resu);

}
/////////////////////////////////////////////////////////////////////////////
//
//            Implementation des fonctions de la classe
//
//             Source_Transport_K_Eps_Realisable_aniso_concen_VDF_Elem
//
/////////////////////////////////////////////////////////////////////////////
//
// remplit les references
void Source_Transport_K_Eps_Realisable_aniso_concen_VDF_Elem::associer_pb(const Probleme_base& pb)
{
  if (pb.nombre_d_equations()<2)
    {
      Cerr<<"The K_Eps_Realisable source term "<<que_suis_je()<<" cannot be activated"<<finl;
      Cerr<<"for a "<<pb.que_suis_je()<<" problem."<<finl;
    }

  const Equation_base& eqn = pb.equation(1);
  const Milieu_base& milieu = pb.equation(0).milieu();
  const Fluide_Incompressible& fluide = ref_cast(Fluide_Incompressible,milieu);

  if (sub_type(Fluide_Quasi_Compressible,fluide))
    {
      Cerr<<"The K_Eps_Realisable source term "<<que_suis_je()<<" cannot be activated"<<finl;
      Cerr<<"with a "<<milieu.que_suis_je()<<" medium."<<finl;
      exit();
    }
  Source_Transport_K_Eps_Realisable_VDF_Elem::associer_pb(pb);

  const Convection_Diffusion_Concentration& eqn_c =
    ref_cast(Convection_Diffusion_Concentration,eqn);
  eq_concentration = eqn_c;
  if (!fluide.beta_c().non_nul())
    {
      Cerr << "You forgot to define beta_co field in the fluid." << finl;
      Cerr << "It is mandatory when using the K-Eps_Realisable model (buoyancy effects)." << finl;
      Cerr << "If you don't want buoyancy effects, then specify: beta_co champ_uniforme 1 0." << finl;
      exit();
    }
  beta_c = fluide.beta_c();
  gravite = fluide.gravite();
}


DoubleTab& Source_Transport_K_Eps_Realisable_aniso_concen_VDF_Elem::ajouter(DoubleTab& resu) const
{
  Source_Transport_K_Eps_Realisable_VDF_Elem::ajouter(resu);
  //
  //// Plutost que de calculer P, on appelle Source_Transport_K_Eps_Realisable_VDF_Elem::ajouter(resu)
  //// et on ajoute directement G
  ////
  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const Zone_Cl_VDF& zcl_VDF_co = ref_cast(Zone_Cl_VDF,eq_concentration->zone_Cl_dis().valeur());
  const DoubleTab& K_eps = eqn_keps_Rea->inconnue().valeurs();
  const DoubleTab& concen = eq_concentration->inconnue().valeurs();
  const Modele_turbulence_scal_base& le_modele_scalaire =
    ref_cast(Modele_turbulence_scal_base,eq_concentration->get_modele(TURBULENCE).valeur());
  const DoubleTab& diffu_turb  = le_modele_scalaire.diffusivite_turbulente().valeurs();
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

  double C1_loc=1.44; // C1 value is not a constant in Realizable K-Epsilon model but here, we take the default value of C1 used in standard K-Epsilon, as proposed by litterature
  double C3_loc;
  double LeK_MIN = eqn_keps_Rea->modele_turbulence().get_LeK_MIN();
  //const Mod_turb_hyd_RANS& mod_turb_RANS = ref_cast(Mod_turb_hyd_RANS,eq_hydraulique->modele_turbulence().valeur());
  //double LeK_MIN = mod_turb_RANS.get_LeK_MIN() ;
  for (int elem=0; elem<nb_elem; elem++)
    {

      resu(elem,0) += G(elem)*volumes(elem)*porosite_vol(elem);

      if (K_eps(elem,0) >= LeK_MIN)
        {
          C3_loc = C3_ ;
          if ( G(elem) > 0. ) C3_loc = 0. ;
          resu(elem,1) += C1_loc*(1.-C3_loc)*G(elem) *volumes(elem)*porosite_vol(elem)
                          * K_eps(elem,1)/K_eps(elem,0);
        }
    }
  return resu;
}

DoubleTab& Source_Transport_K_Eps_Realisable_aniso_concen_VDF_Elem::calculer(DoubleTab& resu) const
{

  resu = 0.;
  return ajouter(resu);

}

/////////////////////////////////////////////////////////////////////////////
//
//            Implementation des fonctions de la classe
//
//             Source_Transport_K_Eps_Realisable_aniso_therm_concen_VDF_Elem
//
/////////////////////////////////////////////////////////////////////////////

// remplit les references eq_thermique, eq_concentration, beta_t_, beta_c_
void Source_Transport_K_Eps_Realisable_aniso_therm_concen_VDF_Elem::associer_pb(const Probleme_base& pb)
{
  if (pb.nombre_d_equations()<3)
    {
      Cerr<<"The K_Eps source term "<<que_suis_je()<<" cannot be activated"<<finl;
      Cerr<<"for a "<<pb.que_suis_je()<<" problem."<<finl;
    }
  const Equation_base& eqn_therm = pb.equation(1);
  const Equation_base& eqn_conc = pb.equation(2);
  const Milieu_base& milieu = eqn_therm.milieu();
  const Fluide_Incompressible& fluide = ref_cast(Fluide_Incompressible,milieu);

  if (sub_type(Fluide_Quasi_Compressible,fluide))
    {
      Cerr<<"The K_Eps_Realisable source term "<<que_suis_je()<<" cannot be activated"<<finl;
      Cerr<<"with a "<<milieu.que_suis_je()<<" medium."<<finl;
      exit();
    }
  Source_Transport_K_Eps_Realisable_VDF_Elem::associer_pb(pb);

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
      Cerr << "It is mandatory when using the K-Eps_Realisable model (buoyancy effects)." << finl;
      Cerr << "If you don't want buoyancy effects, then specify: beta_co champ_uniforme 1 0." << finl;
      exit();
    }
  beta_c = fluide.beta_c();
  gravite = fluide.gravite();
}

DoubleTab& Source_Transport_K_Eps_Realisable_aniso_therm_concen_VDF_Elem::ajouter(DoubleTab& resu) const
{

  Source_Transport_K_Eps_Realisable_VDF_Elem::ajouter(resu);
  //
  //
  //// Plutost que de calculer P, on appelle Source_Transport_K_Eps_Realisable_VDF_Elem::ajouter(resu)
  //// et on ajoute directement G
  ////
  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const Zone_Cl_VDF& zcl_VDF_th = ref_cast(Zone_Cl_VDF,eq_thermique->zone_Cl_dis().valeur());
  const Zone_Cl_VDF& zcl_VDF_co = ref_cast(Zone_Cl_VDF,eq_concentration->zone_Cl_dis().valeur());
  const DoubleTab& K_eps = eqn_keps_Rea->inconnue().valeurs();
  const DoubleTab& temper = eq_thermique->inconnue().valeurs();
  const DoubleTab& concen = eq_concentration->inconnue().valeurs();
  const Modele_turbulence_scal_base& le_modele_scalaire =
    ref_cast(Modele_turbulence_scal_base,eq_thermique->get_modele(TURBULENCE).valeur());
  const DoubleTab& alpha_turb =        le_modele_scalaire.diffusivite_turbulente().valeurs();
  const Modele_turbulence_scal_base& le_modele_scal_co =
    ref_cast(Modele_turbulence_scal_base,eq_concentration->get_modele(TURBULENCE).valeur());
  const DoubleTab& diffu_turb =        le_modele_scal_co.diffusivite_turbulente().valeurs();
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

  double C1_loc=1.44; // C1 value is not a constant in Realizable K-Epsilon model but here, we take the default value of C1 used in standard K-Epsilon, as proposed by litterature
  double C3_loc,G_sum;
  //const Mod_turb_hyd_RANS& mod_turb_RANS = ref_cast(Mod_turb_hyd_RANS,eq_hydraulique->modele_turbulence().valeur());
  //double LeK_MIN = mod_turb_RANS.get_LeK_MIN() ;
  double LeK_MIN = eqn_keps_Rea->modele_turbulence().get_LeK_MIN();

  for (int elem=0; elem<nb_elem; elem++)
    {

      G_sum = G_t(elem)+G_c(elem) ;

      resu(elem,0) += G_sum *volumes(elem)*porosite_vol(elem);

      if (K_eps(elem,0) >= LeK_MIN)
        {
          C3_loc = C3_ ;
          if ( G_sum > 0. ) C3_loc = 0. ;
          resu(elem,1) += C1_loc*(1.-C3_loc)*G_sum *volumes(elem)*porosite_vol(elem)
                          * K_eps(elem,1)/K_eps(elem,0);
        }
    }
  return resu;
}

DoubleTab& Source_Transport_K_Eps_Realisable_aniso_therm_concen_VDF_Elem::calculer(DoubleTab& resu) const
{

  resu = 0.;
  return ajouter(resu);

}
