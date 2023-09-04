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
// File      : Source_Transport_proto.cpp
// Directory : $TURBULENCE_ROOT/src/Specializations/VDF/Sources/new
//
/////////////////////////////////////////////////////////////////////////////

#include <Pb_Thermohydraulique_Concentration_Turbulent.h>
#include <Pb_Hydraulique_Concentration_Turbulent.h>
#include <Pb_Thermohydraulique_Turbulent_QC.h>
#include <Pb_Thermohydraulique_Turbulent.h>
#include <Pb_Hydraulique_Turbulent.h>
#include <Source_Transport_proto.h>
#include <Fluide_Dilatable_base.h>
#include <Param.h>

Entree& Source_Transport_proto::readOn_proto( Entree& is, const Nom& nom)
{
  Param param(nom);
  param.ajouter("C1_eps", &C1);
  param.ajouter("C2_eps", &C2);
  param.lire_avec_accolades(is);
  Cerr << "C1_eps = " << C1 << finl;
  Cerr << "C2_eps = " << C2 << finl;
  return is;
}

Entree& Source_Transport_proto::readOn_nothing(Entree& is, const Nom& nom)
{
  Param param(nom);
  param.lire_avec_accolades(is);
  return is ;
}

Entree& Source_Transport_proto::readOn_anisotherme(Entree& is, const Nom& nom)
{
  Param param(nom);
  param.ajouter("C1_eps", &C1);
  param.ajouter("C2_eps", &C2);
  param.ajouter("C3_eps", &C3);
  param.lire_avec_accolades(is);
  Cerr << "C1_eps = " << C1 << finl;
  Cerr << "C2_eps = " << C2 << finl;
  Cerr << "C3_eps = " << C3 << finl;
  return is ;
}

Entree& Source_Transport_proto::readOn_concen(Entree& is, const Nom& nom) { return readOn_anisotherme(is,nom); }
Entree& Source_Transport_proto::readOn_anisotherme_concen(Entree& is, const Nom& nom) { return readOn_anisotherme(is,nom); }

Entree& Source_Transport_proto::readOn_real(Entree& is, const Nom& nom)
{
  Param param(nom);
  param.ajouter("C2_eps", &C2);
  param.ajouter("interpolation_viscosite_turbulente", &_interpolation_viscosite_turbulente);
  param.ajouter("coefficient_limiteur", &_coefficient_limiteur);
  param.lire_avec_accolades(is);
  Cerr << "C2_eps = " << C2 << finl;

  // Checking of the value given in the data deck for "interpolation_viscosite_turbulente"
  if ( _interpolation_viscosite_turbulente == 0 ) { Cerr << "Interpolation arithmetique de la viscosite turbulente aux faces (si VEF) = " << finl; }
  else if ( _interpolation_viscosite_turbulente == 1 ) { Cerr << "Interpolation harmonique de la viscosite turbulente aux faces (si VEF) = " << finl; }
  else if ( _interpolation_viscosite_turbulente == 2 ) { Cerr << "Interpolation harmonique ponderee par les volumes de maille de la viscosite turbulente aux faces (si VEF) = " << finl; }
  else if ( _interpolation_viscosite_turbulente == 3 ) { Cerr << "Interpolation harmonique ponderee par les tailles de maille de la viscosite turbulente aux faces (si VEF) = " << finl; }
  else
    {
      Cerr << "Error in 'interpolation_viscosite_turbulente' input value :" << _interpolation_viscosite_turbulente << finl;
      Process::exit();
    }

  return is;
}

Entree& Source_Transport_proto::readOn_anisotherme_real(Entree& is, const Nom& nom)
{
  Param param(nom);
  param.ajouter("C2_eps", &C2);
  param.ajouter("C3_eps", &C3);
  param.lire_avec_accolades(is);
  Cerr << "C2_eps = " << C2 << finl;
  Cerr << "C3_eps = " << C3 << finl;
  return is;
}

Entree& Source_Transport_proto::readOn_concen_real(Entree& is, const Nom& nom) { return readOn_anisotherme_real(is,nom); }
Entree& Source_Transport_proto::readOn_anisotherme_concen_real(Entree& is, const Nom& nom) { return readOn_anisotherme_real(is,nom); }

void Source_Transport_proto::verifier_pb_keps(const Probleme_base& pb, const Nom& nom)
{
  if (!sub_type(Pb_Hydraulique_Turbulent,pb) && !sub_type(Pb_Thermohydraulique_Turbulent_QC,pb)) error_keps(nom,pb.que_suis_je());
}

void Source_Transport_proto::verifier_pb_keps_anisotherme(const Probleme_base& pb, const Nom& nom)
{
  if (!sub_type(Pb_Thermohydraulique_Turbulent,pb)) error_keps(nom,pb.que_suis_je());
}

void Source_Transport_proto::verifier_pb_keps_concen(const Probleme_base& pb, const Nom& nom)
{
  if (!sub_type(Pb_Hydraulique_Concentration_Turbulent,pb)) error_keps(nom,pb.que_suis_je());
}

void Source_Transport_proto::verifier_pb_keps_anisotherme_concen(const Probleme_base& pb, const Nom& nom)
{
  if (!sub_type(Pb_Thermohydraulique_Concentration_Turbulent,pb)) error_keps(nom,pb.que_suis_je());
}

void Source_Transport_proto::verifier_milieu_anisotherme(const Probleme_base& pb, const Nom& nom)
{
  const Milieu_base& milieu = pb.equation(1).milieu(); // eq thermique
  if (pb.nombre_d_equations()<2) error_keps(nom,pb.que_suis_je());
  if (sub_type(Fluide_Dilatable_base,ref_cast(Fluide_base,milieu))) error_keps(nom,milieu.que_suis_je());
}

void Source_Transport_proto::verifier_milieu_concen(const Probleme_base& pb, const Nom& nom)
{
  const Milieu_base& milieu = pb.equation(0).milieu(); // XXX : Attention pas eq 1 car Constituant derive pas de Fluide_base ! donc eq hydro
  if (pb.nombre_d_equations()<2) error_keps(nom,pb.que_suis_je());
  if (sub_type(Fluide_Dilatable_base,ref_cast(Fluide_base,milieu))) error_keps(nom,milieu.que_suis_je());
}

void Source_Transport_proto::verifier_milieu_anisotherme_concen(const Probleme_base& pb, const Nom& nom)
{
  const Milieu_base& milieu = pb.equation(1).milieu(); // eq thermique
  if (pb.nombre_d_equations()<3) error_keps(nom,pb.que_suis_je());
  if (sub_type(Fluide_Dilatable_base,ref_cast(Fluide_base,milieu))) error_keps(nom,milieu.que_suis_je());
}

void Source_Transport_proto::verifier_beta_concen(const Fluide_base& fluide)
{
  if (!fluide.beta_c().non_nul())
    {
      Cerr << "You forgot to define beta_co field in the fluid. It is mandatory when using the K-Eps model (buoyancy effects)." << finl;
      Cerr << "If you don't want buoyancy effects, then specify: beta_co champ_uniforme 1 0." << finl;
      Process::exit();
    }
}

void Source_Transport_proto::associer_pb_proto(const Probleme_base& pb)
{
  eq_hydraulique = pb.equation(0);
}

void Source_Transport_proto::associer_pb_anisotherme(const Probleme_base& pb)
{
  const Fluide_base& fluide = ref_cast(Fluide_base,pb.equation(1).milieu());
  beta_t = fluide.beta_t();
  gravite = fluide.gravite();
  eq_thermique = ref_cast(Convection_Diffusion_Temperature,pb.equation(1));
}

void Source_Transport_proto::associer_pb_concen(const Probleme_base& pb)
{
  const Fluide_base& fluide = ref_cast(Fluide_base,pb.equation(0).milieu()); // XXX : Attention pas eq 1 car Constituant derive pas de Fluide_base !
  verifier_beta_concen(fluide);
  beta_c = fluide.beta_c();
  gravite = fluide.gravite();
  eq_concentration = ref_cast(Convection_Diffusion_Concentration,pb.equation(1));
}

void Source_Transport_proto::associer_pb_anisotherme_concen(const Probleme_base& pb)
{
  const Fluide_base& fluide = ref_cast(Fluide_base,pb.equation(1).milieu()); // a partir de l'eq thermique
  verifier_beta_concen(fluide);
  beta_t = fluide.beta_t();
  beta_c = fluide.beta_c();
  gravite = fluide.gravite();
  eq_thermique = ref_cast(Convection_Diffusion_Temperature,pb.equation(1));
  eq_concentration = ref_cast(Convection_Diffusion_Concentration,pb.equation(2));
}
