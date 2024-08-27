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
// File:        Echelle_temporelle_turbulente.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Equations
// Version:     /main/52
//
//////////////////////////////////////////////////////////////////////////////

#include <Echelle_temporelle_turbulente.h>
#include <Pb_Multiphase.h>
#include <Discret_Thyd.h>
#include <Domaine_VF.h>
#include <Domaine.h>
#include <Avanc.h>
#include <Debog.h>
#include <Frontiere_dis_base.h>
#include <EcritureLectureSpecial.h>
#include <Champ_Uniforme.h>
#include <Matrice_Morse.h>
#include <Navier_Stokes_std.h>
#include <TRUSTTrav.h>
#include <Neumann_sortie_libre.h>
#include <Op_Conv_negligeable.h>
#include <Param.h>
#include <Schema_Implicite_base.h>
#include <SETS.h>
#include <EChaine.h>
#include <Scalaire_impose_paroi.h>
#include <Echange_global_impose.h>

#define old_forme

Implemente_instanciable(Echelle_temporelle_turbulente,"Echelle_temporelle_turbulente",Convection_diffusion_turbulence_multiphase);


/*! @brief Simple appel a: Convection_diffusion_turbulence_multiphase::printOn(Sortie&)
 *
 * @param (Sortie& is) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Echelle_temporelle_turbulente::printOn(Sortie& is) const
{
  return Convection_diffusion_turbulence_multiphase::printOn(is);
}

/*! @brief Verifie si l'equation a une inconnue et un fluide associe et appelle Convection_Diffusion_std::readOn(Entree&).
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree& is) le flot d'entree modifie
 */
Entree& Echelle_temporelle_turbulente::readOn(Entree& is)
{
  Convection_diffusion_turbulence_multiphase::readOn(is);

  terme_convectif.set_fichier("Convection_echelle_temporelle_turbulente");
  terme_convectif.set_description((Nom)"Turbulent time scale transfer rate=Integral(-tau*ndS)");
  terme_diffusif.set_fichier("Diffusion_echelle_temporelle_turbulente");
  terme_diffusif.set_description((Nom)"Turbulent time scale transfer rate=Integral(nu*grad(tau)*ndS)");
  return is;
}

const Champ_Don& Echelle_temporelle_turbulente::diffusivite_pour_transport() const
{
  return ref_cast(Fluide_base,milieu()).viscosite_cinematique();
}

const Champ_base& Echelle_temporelle_turbulente::diffusivite_pour_pas_de_temps() const
{
  return ref_cast(Fluide_base,milieu()).viscosite_cinematique();
}

/*! @brief Discretise l'equation.
 *
 */
void Echelle_temporelle_turbulente::discretiser()
{
  int nb_valeurs_temp = schema_temps().nb_valeurs_temporelles();
  double temps = schema_temps().temps_courant();
  const Discret_Thyd& dis=ref_cast(Discret_Thyd, discretisation());
  Cerr << "Turbulent time scale discretization" << finl;
  //On utilise temperature pour la directive car discretisation identique
  dis.discretiser_champ("temperature",domaine_dis(),"tau","s", 1,nb_valeurs_temp,temps,l_inco_ch);//une seule compo, meme en multiphase
  l_inco_ch->fixer_nature_du_champ(scalaire);
  l_inco_ch->fixer_nom_compo(0, Nom("tau"));
  champs_compris_.ajoute_champ(l_inco_ch);
  Equation_base::discretiser();
  Cerr << "Echelle_temporelle_turbulente::discretiser() ok" << finl;
}

void Echelle_temporelle_turbulente::calculer_tau(const Objet_U& obj, DoubleTab& val, DoubleTab& bval, tabs_t& deriv)
{
  const Equation_base& eqn = ref_cast(Equation_base, obj);
  const DoubleTab& tau = eqn.inconnue()->valeurs();

  /* valeurs du champ */
  int i, n, N = val.line_size(), Nl = val.dimension_tot(0);
  for (i = 0; i < Nl; i++)
    for (n = 0; n < N; n++) val(i, n) = tau(i, n);

  /* on ne peut utiliser valeur_aux_bords que si ch_rho a un domaine_dis_base */
  const DoubleTab& b_tau = eqn.inconnue()->valeur_aux_bords();
  int Nb = b_tau.dimension_tot(0);
  for (i = 0; i < Nb; i++)
    for (n = 0; n < N; n++) bval(i, n) = b_tau(i, n);

  //derivee en tau : 1
  DoubleTab& d_tau = deriv["tau"];
  for (d_tau.resize(Nl, N), i = 0; i < Nl; i++)
    for (n = 0; n < N; n++) d_tau(i, n) = 1.;
}
