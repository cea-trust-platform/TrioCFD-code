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
// File:        Fluide_Diphasique.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/17
//
//////////////////////////////////////////////////////////////////////////////
#include <Fluide_Diphasique.h>
#include <Fluide_Incompressible.h>
#include <Motcle.h>
#include <Interprete.h>
#include <Ref_Fluide_Diphasique.h>
#include <Param.h>
#include <Champ_Uniforme.h>

Implemente_instanciable_sans_constructeur(Fluide_Diphasique,"Fluide_Diphasique",Milieu_base);

Implemente_ref(Fluide_Diphasique);
//////////// XD fluide_diphasique milieu_v2_base fluide_diphasique -1 Two-phase fluid.
// XD fluide_diphasique milieu_base fluide_diphasique -1 Two-phase fluid.

Fluide_Diphasique::Fluide_Diphasique()
{
  indic_rayo_ = NONRAYO;
}

Sortie& Fluide_Diphasique::printOn(Sortie& os) const
{
  return os;
}

// cf Milieu_base::readOn
Entree& Fluide_Diphasique::readOn(Entree& is)
{
  // Default value for formule_mu_
  formule_mu_="standard";
  Milieu_base::readOn(is);
  return is;
}

void Fluide_Diphasique::set_param(Param& param)
{
  param.ajouter("sigma",&sigma_,Param::REQUIRED); // XD_ADD_P champ_don_base surfacic tension (J/m2)
  param.ajouter_non_std("fluide0",(this),Param::REQUIRED); // XD_ADD_P chaine first phase fluid
  param.ajouter_non_std("fluide1",(this),Param::REQUIRED); // XD_ADD_P chaine second phase fluid
  param.ajouter("chaleur_latente",&chaleur_latente_); // XD_ADD_P champ_don_base phase changement enthalpy h(phase1_) - h(phase0_) (J/kg/K)
  param.ajouter("formule_mu",&formule_mu_); // XD_ADD_P chaine (into=[standard,arithmetic,harmonic]) formula used to calculate average
  Milieu_base::set_additional_params(param);
}

int Fluide_Diphasique::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if ((mot=="fluide0") || (mot=="fluide1"))
    {
      Nom nom_objet;
      is >> nom_objet;
      const Objet_U& objet = Interprete::objet(nom_objet);
      const Fluide_Incompressible& fluide = ref_cast(Fluide_Incompressible,objet);
      if (mot=="fluide0")
        phase0_ = fluide;
      else if (mot=="fluide1")
        phase1_ = fluide;
      return 1;
    }
  else
    return Milieu_base::lire_motcle_non_standard(mot,is);
  return 1;
}

void Fluide_Diphasique::verifier_coherence_champs(int& err,Nom& msg)
{
  msg="";
  if (!sub_type(Champ_Uniforme,sigma_.valeur()))
    {
      msg += "The surface tension sigma must be specify with a Champ_Uniforme type field. \n";
      err = 1;
    }
  else
    {
      if (sigma_(0,0) < 0)
        {
          msg += "The surface tension sigma is not positive. \n";
          err = 1;
        }
    }

  if (chaleur_latente_.non_nul())
    {
      if (!sub_type(Champ_Uniforme,chaleur_latente_.valeur()))
        {
          msg += "The latent heat chaleur_latente must be specify with a Champ_Uniforme type field. \n";
          err = 1;
        }
    }
  Milieu_base::verifier_coherence_champs(err,msg);
}

const Fluide_Incompressible&
Fluide_Diphasique::fluide_phase(int phase) const
{
  assert(phase == 0 || phase == 1);
  if (phase == 0)
    return phase0_;
  else
    return phase1_;
}

double Fluide_Diphasique::sigma() const
{
  return sigma_(0,0);
}

double Fluide_Diphasique::chaleur_latente() const
{
  if (!chaleur_latente_.non_nul())
    {
      Cerr << "Fluide_Diphasique::chaleur_latente()\n";
      Cerr << " The latent heat has not been specified." << finl;
      exit();
    }
  return chaleur_latente_(0,0);
}

int Fluide_Diphasique::formule_mu() const
// These values are used in the switch of Navier_Stokes_FT_Disc::FT_disc_calculer_champs_rho_mu_nu_dipha
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
  phase0_.initialiser(temps);
  phase1_.initialiser(temps);
  initialiser_porosite(temps);
  return 1;
}

void Fluide_Diphasique::mettre_a_jour(double temps)
{
  // en particulier, mise a jour de g qui peut dependre de t...
  Milieu_base::mettre_a_jour(temps);
  phase0_.mettre_a_jour(temps);
  phase1_.mettre_a_jour(temps);
}
void Fluide_Diphasique::discretiser(const Probleme_base& pb, const  Discretisation_base& dis)
{
  // sigma chaleur latente phase_0 phase_1  diffusivite revoir
  phase0_.discretiser(pb,dis);
  phase1_.discretiser(pb,dis);
  discretiser_porosite(pb,dis);
  discretiser_diametre_hydro(pb, dis);
}



// ===========================================================================
// Appels invalides :

const Champ& Fluide_Diphasique::masse_volumique() const
{
  Cerr << "Error : Fluide_Diphasique::masse_volumique()" << finl;
  assert(0);
  exit();
  throw;
  return phase0_.masse_volumique();
}

Champ& Fluide_Diphasique::masse_volumique()
{
  Cerr << "Error : Fluide_Diphasique::masse_volumique()" << finl;
  assert(0);
  exit();
  throw;
  return masse_volumique();
}

const Champ_Don& Fluide_Diphasique::diffusivite() const
{
  Cerr << "Error : Fluide_Diphasique::diffusivite()" << finl;
  assert(0);
  exit();
  throw;
  return diffusivite();
}

Champ_Don& Fluide_Diphasique::diffusivite()
{
  Cerr << "Error : Fluide_Diphasique::diffusivite()" << finl;
  assert(0);
  exit();
  throw;
  return diffusivite();
}

const Champ_Don& Fluide_Diphasique::conductivite() const
{
  Cerr << "Error : Fluide_Diphasique::conductivite()" << finl;
  assert(0);
  exit();
  throw;
  return conductivite();
}

Champ_Don& Fluide_Diphasique::conductivite()
{
  Cerr << "Error : Fluide_Diphasique::conductivite()" << finl;
  assert(0);
  exit();
  throw;
  return conductivite();
}

const Champ_Don& Fluide_Diphasique::capacite_calorifique() const
{
  Cerr << "Error : Fluide_Diphasique::capacite_calorifique()" << finl;
  assert(0);
  exit();
  throw;
  return capacite_calorifique();
}

Champ_Don& Fluide_Diphasique::capacite_calorifique()
{
  Cerr << "Error : Fluide_Diphasique::capacite_calorifique()" << finl;
  assert(0);
  exit();
  throw;
  return capacite_calorifique();
}

const Champ_Don& Fluide_Diphasique::beta_t() const
{
  Cerr << "Error : Fluide_Diphasique::beta_t()" << finl;
  assert(0);
  exit();
  throw;
  return beta_t();
}

Champ_Don& Fluide_Diphasique::beta_t()
{
  Cerr << "Error : Fluide_Diphasique::beta_t()" << finl;
  assert(0);
  exit();
  throw;
  return beta_t();
}

