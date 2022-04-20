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
// File      : Source_Transport_Eps_aniso_concen_VEF_Face.cpp
// Directory : $TURBULENCE_ROOT/src/Specializations/VEF/Sources/RANS
//
/////////////////////////////////////////////////////////////////////////////

#include <Source_Transport_Eps_aniso_concen_VEF_Face.h>
#include <Paroi_negligeable_VEF.h>
#include <Mod_turb_hyd_base.h>
#include <Transport_K_ou_Eps.h>
#include <Convection_Diffusion_Temperature.h>
#include <Convection_Diffusion_Concentration.h>
#include <Modele_turbulence_scal_base.h>
#include <Fluide_base.h>
#include <Probleme_base.h>
#include <Champ_Uniforme.h>
#include <Zone_VEF.h>
#include <Champ_P1NC.h>
#include <TRUSTTrav.h>
#include <Fluide_Quasi_Compressible.h>
#include <Pb_Hydraulique_Turbulent.h>
#include <Pb_Hydraulique_Concentration_Turbulent.h>
#include <Pb_Thermohydraulique_Turbulent_QC.h>
#include <Pb_Thermohydraulique_Turbulent.h>
#include <Pb_Thermohydraulique_Concentration_Turbulent.h>
#include <Param.h>
#include <Constituant.h>

#include <Modele_turbulence_hyd_K_Eps_Bicephale.h>
Implemente_instanciable_sans_constructeur(Source_Transport_Eps_aniso_concen_VEF_Face,"Source_Transport_Eps_aniso_concen_VEF_P1NC",Source_Transport_Eps_VEF_Face);


Sortie& Source_Transport_Eps_aniso_concen_VEF_Face::printOn(Sortie& s) const
{
  return s << que_suis_je() ;
}


Entree& Source_Transport_Eps_aniso_concen_VEF_Face::readOn(Entree& is)
{
  const Probleme_base& problem = mon_equation->probleme();
  if (!sub_type(Pb_Hydraulique_Concentration_Turbulent,problem)) error(que_suis_je(),problem.que_suis_je());
  Param param(que_suis_je());
  param.ajouter("C1_eps", &C1);
  param.ajouter("C2_eps", &C2);
  param.ajouter("C3_eps", &C3);
  param.lire_avec_accolades(is);
  Cerr << "C1_eps = " << C1 << finl;
  Cerr << "C2_eps = " << C2 << finl;
  Cerr << "C3_eps = " << C3 << finl;
  return is;
}


void Source_Transport_Eps_aniso_concen_VEF_Face::associer_pb(const Probleme_base& pb)
{
  if (pb.nombre_d_equations()<2)
    {
      Cerr<<"The Eps source term "<<que_suis_je()<<" cannot be activated"<<finl;
      Cerr<<"for a "<<pb.que_suis_je()<<" problem."<<finl;
    }

  const Equation_base& eqn = pb.equation(1);
  const Milieu_base& milieu = pb.equation(0).milieu();
  const Fluide_base& fluide = ref_cast(Fluide_base,milieu);

  if (sub_type(Fluide_Quasi_Compressible,fluide))
    {
      Cerr<<"The Eps source term "<<que_suis_je()<<" cannot be activated"<<finl;
      Cerr<<"with a "<<milieu.que_suis_je()<<" medium."<<finl;
      exit();
    }

  Source_Transport_Eps_VEF_Face::associer_pb(pb);

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


DoubleTab& Source_Transport_Eps_aniso_concen_VEF_Face::calculer(DoubleTab& resu) const
{
  resu = 0;
  return ajouter(resu);
}



DoubleTab& Source_Transport_Eps_aniso_concen_VEF_Face::ajouter(DoubleTab& resu) const
{
  Source_Transport_Eps_VEF_Face::ajouter(resu);
  //
  // Modifs VB : plutot que de calculer P, on appelle Source_Transport_Eps_VEF_Face::ajouter(resu)
  // et on ajoute directement G
  // On en profite pour faire des tests sur LeK_MIN
  //
  const Zone_VEF& zone_VEF = la_zone_VEF.valeur();
  const Zone_Cl_VEF& zcl_VEF_co = ref_cast(Zone_Cl_VEF,eq_concentration->zone_Cl_dis().valeur());
  const DoubleTab& K   = mon_eq_transport_K->inconnue().valeurs();
  const DoubleTab& Eps = mon_eq_transport_Eps->inconnue().valeurs();
  const DoubleTab& concen = eq_concentration->inconnue().valeurs();
  const Modele_turbulence_scal_base& le_modele_scalaire =
    ref_cast(Modele_turbulence_scal_base,eq_concentration->get_modele(TURBULENCE).valeur());
  const DoubleTab& lambda_turb = le_modele_scalaire.conductivite_turbulente().valeurs();
  const DoubleVect& g = gravite->valeurs();
  const Champ_Don& ch_beta_concen = beta_c.valeur();

  const DoubleVect& volumes_entrelaces = zone_VEF.volumes_entrelaces();
  int nb_face = zone_VEF.nb_faces();

  DoubleTrav G(nb_face);

  int nb_consti = eq_concentration->constituant().nb_constituants();

  calculer_terme_destruction_K_gen(zone_VEF,zcl_VEF_co,G,concen,lambda_turb,ch_beta_concen,g,nb_consti);

  double C3_loc;
  double LeK_MIN = mon_eq_transport_Eps->modele_turbulence().get_LeK_MIN() ;
  for (int face=0; face<nb_face; face++)
    {
      if (K(face) >= LeK_MIN)
        {
          C3_loc = C3 ;
          if ( G(face) > 0.) C3_loc = 0 ;
          resu(face) += C1*(1-C3_loc)*G(face)*volumes_entrelaces(face)
                        *Eps(face)/K(face);

        }
    }
  return resu;
}










