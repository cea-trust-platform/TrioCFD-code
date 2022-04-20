/****************************************************************************
* Copyright (c) 2022, CEA
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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Source_Transport_K_Realisable_aniso_therm_concen_VEF_Face.cpp
// Directory : $TURBULENCE_ROOT/src/Specializations/VEF/Sources/RANS
//
/////////////////////////////////////////////////////////////////////////////

#include <Source_Transport_K_Realisable_aniso_therm_concen_VEF_Face.h>
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
#include <TRUSTTrav.h>
#include <Fluide_Quasi_Compressible.h>
#include <Modele_Shih_Zhu_Lumley_VEF.h>
#include <Pb_Hydraulique_Turbulent.h>
#include <Pb_Hydraulique_Concentration_Turbulent.h>
#include <Pb_Thermohydraulique_Turbulent_QC.h>
#include <Pb_Thermohydraulique_Turbulent.h>
#include <Pb_Thermohydraulique_Concentration_Turbulent.h>
#include <Constituant.h>
#include <Param.h>

Implemente_instanciable(Source_Transport_K_Realisable_aniso_therm_concen_VEF_Face,"Source_Transport_K_Realisable_aniso_therm_concen_VEF_P1NC",Source_Transport_K_Realisable_VEF_Face);

Sortie& Source_Transport_K_Realisable_aniso_therm_concen_VEF_Face::printOn(Sortie& s) const
{
  return s << que_suis_je() ;
}

Entree& Source_Transport_K_Realisable_aniso_therm_concen_VEF_Face::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.lire_avec_accolades(is);
  return is ;
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



DoubleTab& Source_Transport_K_Realisable_aniso_therm_concen_VEF_Face::calculer(DoubleTab& resu) const
{
  resu = 0;
  return ajouter(resu);
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


















