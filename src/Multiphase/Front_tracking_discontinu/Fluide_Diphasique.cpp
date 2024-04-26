/****************************************************************************
* Copyright (c) 2024, CEA
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
// File:        Fluide_Diphasique.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/17
//
//////////////////////////////////////////////////////////////////////////////

#include <Fluide_Diphasique.h>
#include <Champ_Uniforme.h>
#include <Interprete.h>
#include <TRUST_Ref.h>
#include <Motcle.h>
#include <Param.h>

// XD fluid_diph_lu objet_lecture nul 0 Single fluid to be read.
// XD  attr fluid_name chaine fluid_name 0 Name of the fluid which is part of the diphasic fluid.
// XD  attr single_fld fluide_incompressible single_fld 0 Definition of the single fluid part of a multiphasic fluid.

Implemente_instanciable_sans_constructeur(Fluide_Diphasique, "Fluide_Diphasique", Milieu_base);
// XD fluide_diphasique milieu_base fluide_diphasique -1 fluid_diph_lu 0 Two-phase fluid.

Sortie& Fluide_Diphasique::printOn(Sortie& os) const { return os; }

Entree& Fluide_Diphasique::readOn(Entree& is) { return Milieu_base::readOn(is); }

void Fluide_Diphasique::set_param(Param& param)
{
  param.ajouter("sigma", &sigma_, Param::REQUIRED); // XD_ADD_P champ_don_base surfacic tension (J/m2)
  param.ajouter("fluide0|phase0", &phase0_, Param::REQUIRED); // XD_ADD_P fluid_diph_lu first phase fluid
  param.ajouter("fluide1|phase1", &phase1_, Param::REQUIRED); // XD_ADD_P fluid_diph_lu second phase fluid
  param.ajouter("chaleur_latente", &chaleur_latente_); // XD_ADD_P champ_don_base phase changement enthalpy h(phase1_) - h(phase0_) (J/kg/K)
  param.ajouter("formule_mu", &formule_mu_); // XD_ADD_P chaine (into=[standard,arithmetic,harmonic]) formula used to calculate average
  Milieu_base::set_additional_params(param); // XD ref gravite field_base
}

void Fluide_Diphasique::verifier_coherence_champs(int& err, Nom& msg)
{
  msg = "";
  if (!sub_type(Fluide_Incompressible, phase0_.valeur()) && !sub_type(Fluide_Incompressible, phase1_.valeur()))
    {
      msg += "Both phases defined in the bloc Fluide_Diphasique must be of type Fluide_Incompressible. \n";
      err = 1;
    }

  if (phase0_->a_gravite() || phase1_->a_gravite())
    {
      msg += "The gravity field should be defined in the Fluide_Diphasique bloc and not inside the Incompressible fluids bloc. \n";
      err = 1;
    }

  if (!sub_type(Champ_Uniforme, sigma_.valeur()))
    {
      msg += "The surface tension sigma must be specify with a Champ_Uniforme type field. \n";
      err = 1;
    }
  else
    {
      if (sigma_(0, 0) < 0)
        {
          msg += "The surface tension sigma is not positive. \n";
          err = 1;
        }
    }

  if (chaleur_latente_.non_nul())
    if (!sub_type(Champ_Uniforme, chaleur_latente_.valeur()))
      {
        msg += "The latent heat chaleur_latente must be specify with a Champ_Uniforme type field. \n";
        err = 1;
      }

  Milieu_base::verifier_coherence_champs(err, msg);
}

const Fluide_Incompressible& Fluide_Diphasique::fluide_phase(int phase) const
{
  assert(phase == 0 || phase == 1);
  return (phase == 0) ? ref_cast(Fluide_Incompressible, phase0_.valeur()) : ref_cast(Fluide_Incompressible, phase1_.valeur());
}

double Fluide_Diphasique::sigma() const
{
  return sigma_(0, 0);
}

double Fluide_Diphasique::chaleur_latente() const
{
  if (!chaleur_latente_.non_nul())
    {
      Cerr << "Fluide_Diphasique::chaleur_latente() : The latent heat has not been specified." << finl;
      exit();
    }
  return chaleur_latente_(0, 0);
}

// These values are used in the switch of Navier_Stokes_FT_Disc::FT_disc_calculer_champs_rho_mu_nu_dipha
int Fluide_Diphasique::formule_mu() const
{
  if (formule_mu_ == "standard")
    return 0;
  else if ((formule_mu_ == "arithmetique") or (formule_mu_ == "arithmetic"))
    return 1;
  else if ((formule_mu_ == "harmonique") or (formule_mu_ == "harmonic"))
    return 2;
  else
    return -1;
}

int Fluide_Diphasique::initialiser(const double temps)
{
  phase0_->initialiser(temps);
  phase1_->initialiser(temps);
  initialiser_porosite(temps);
  return 1;
}

void Fluide_Diphasique::mettre_a_jour(double temps)
{
  Milieu_base::mettre_a_jour(temps);
  phase0_->mettre_a_jour(temps);
  phase1_->mettre_a_jour(temps);
}
void Fluide_Diphasique::discretiser(const Probleme_base& pb, const Discretisation_base& dis)
{
  phase0_->discretiser(pb, dis);
  phase1_->discretiser(pb, dis);
  discretiser_porosite(pb, dis);
  discretiser_diametre_hydro(pb, dis);
}
