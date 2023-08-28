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
// File:        Energie_cinetique_turbulente.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Equations
// Version:     /main/52
//
//////////////////////////////////////////////////////////////////////////////

#include <Energie_cinetique_turbulente.h>
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
#include <Neumann_paroi.h>
#include <Scalaire_impose_paroi.h>
#include <Echange_global_impose.h>

#define old_forme

Implemente_instanciable(Energie_cinetique_turbulente,"Energie_cinetique_turbulente",Convection_diffusion_turbulence_multiphase);

/*! @brief Simple appel a: Convection_diffusion_turbulence_multiphase::printOn(Sortie&)
 *
 * @param (Sortie& is) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Energie_cinetique_turbulente::printOn(Sortie& is) const
{
  return Convection_diffusion_turbulence_multiphase::printOn(is);
}

/*! @brief Appelle Convection_diffusion_turbulence_multiphase::readOn(Entree&).
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree& is) le flot d'entree modifie
 */
Entree& Energie_cinetique_turbulente::readOn(Entree& is)
{
  Convection_diffusion_turbulence_multiphase::readOn(is);
  terme_convectif.set_fichier("Convection_energie_cinetique_turbulente");
  terme_convectif.set_description((Nom)"Turbulent kinetic energy transfer rate=Integral(-rho*k*ndS) [W] if SI units used");
  terme_diffusif.set_fichier("Diffusion_energie_cinetique_turbulente");
  terme_diffusif.set_description((Nom)"Turbulent kinetic energy transfer rate=Integral(mu*grad(k)*ndS) [W] if SI units used");
  return is;
}

const Champ_Don& Energie_cinetique_turbulente::diffusivite_pour_transport() const
{
  return ref_cast(Fluide_base,milieu()).viscosite_dynamique();
}

const Champ_base& Energie_cinetique_turbulente::diffusivite_pour_pas_de_temps() const
{
  return ref_cast(Fluide_base,milieu()).viscosite_cinematique();
}

/*! @brief Discretise l'equation.
 *
 */
void Energie_cinetique_turbulente::discretiser()
{
  int nb_valeurs_temp = schema_temps().nb_valeurs_temporelles();
  double temps = schema_temps().temps_courant();
  const Discret_Thyd& dis=ref_cast(Discret_Thyd, discretisation());
  Cerr << "Turbulent kinetic energy discretization" << finl;
  //On utilise temperature pour la directive car discretisation identique
  dis.discretiser_champ("temperature",domaine_dis(),"k","J/kg", 1,nb_valeurs_temp,temps,l_inco_ch);//une seule compo, meme en multiphase
  l_inco_ch.valeur().fixer_nature_du_champ(scalaire);
  l_inco_ch.valeur().fixer_nom_compo(0, Nom("k"));
  champs_compris_.ajoute_champ(l_inco_ch);
  Equation_base::discretiser();
  Cerr << "Energie_cinetique_turbulente::discretiser() ok" << finl;
}

void Energie_cinetique_turbulente::calculer_alpha_rho_k(const Objet_U& obj, DoubleTab& val, DoubleTab& bval, tabs_t& deriv)
{
  const Equation_base& eqn = ref_cast(Equation_base, obj);
  const Fluide_base& fl = ref_cast(Fluide_base, eqn.milieu());
  const Champ_base& ch_rho = fl.masse_volumique();
  const Champ_Inc_base *ch_alpha = sub_type(Pb_Multiphase, eqn.probleme()) ? &ref_cast(Pb_Multiphase, eqn.probleme()).equation_masse().inconnue().valeur() : NULL,
                        *pch_rho = sub_type(Champ_Inc_base, ch_rho) ? &ref_cast(Champ_Inc_base, ch_rho) : NULL; //pas toujours un Champ_Inc
  const DoubleTab* alpha = ch_alpha ? &ch_alpha->valeurs() : NULL, &rho = ch_rho.valeurs(), &k = eqn.inconnue().valeurs();

  /* valeurs du champ */
  int i, n, N = val.line_size(), Nl = val.dimension_tot(0), cR = sub_type(Champ_Uniforme, ch_rho);
  for (i = 0; i < Nl; i++)
    for (n = 0; n < N; n++) val(i, n) = (alpha ? (*alpha)(i, n) : 1) * rho(!cR * i, n) * k(i, n);

  /* on ne peut utiliser valeur_aux_bords que si ch_rho a un domaine_dis_base */
  DoubleTab b_al = ch_alpha ? ch_alpha->valeur_aux_bords() : DoubleTab();
  DoubleTab b_rho, b_k = eqn.inconnue()->valeur_aux_bords();
  int Nb = b_k.dimension_tot(0);
  if (ch_rho.a_un_domaine_dis_base()) b_rho = ch_rho.valeur_aux_bords();
  else b_rho.resize(Nb, rho.line_size()), ch_rho.valeur_aux(ref_cast(Domaine_VF, eqn.domaine_dis().valeur()).xv_bord(), b_rho);
  for (i = 0; i < Nb; i++)
    for (n = 0; n < N; n++) bval(i, n) = (alpha ? b_al(i, n) : 1) * b_rho(i, n) * b_k(i, n);

  if (alpha)//derivee en alpha : rho * k
    {
      DoubleTab& d_a = deriv["alpha"];
      for (d_a.resize(Nl, N), i = 0; i < Nl; i++)
        for (n = 0; n < N; n++) d_a(i, n) = rho(!cR * i, n) * k(i, n);
    }
  //derivee en k : alpha * rho
  DoubleTab& d_k = deriv["k"];
  for (d_k.resize(Nl, N), i = 0; i < Nl; i++)
    for (n = 0; n < N; n++) d_k(i, n) = (alpha ? (*alpha)(i, n) : 1) * rho(!cR * i, n);

  /* derivees a travers rho */
  if (pch_rho)
    for (auto && n_d :pch_rho->derivees())
      {
        DoubleTab& d_v = deriv[n_d.first];
        for (d_v.resize(Nl, N), i = 0; i < Nl; i++)
          for (n = 0; n < N; n++)
            d_v(i, n) = (alpha ? (*alpha)(i, n) : 1) * k(i, n) * n_d.second(i, n);
      }
}
