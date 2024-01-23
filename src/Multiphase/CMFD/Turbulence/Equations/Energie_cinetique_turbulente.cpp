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

#include <Energie_cinetique_turbulente.h>
#include <EcritureLectureSpecial.h>
#include <Scalaire_impose_paroi.h>
#include <Champ_Face_PolyMAC_P0.h>
#include <Echange_global_impose.h>
#include <Schema_Implicite_base.h>
#include <Neumann_sortie_libre.h>
#include <Op_Conv_negligeable.h>
#include <Frontiere_dis_base.h>
#include <Navier_Stokes_std.h>
#include <Champ_Uniforme.h>
#include <Matrice_Morse.h>
#include <Neumann_paroi.h>
#include <Probleme_base.h>
#include <Discret_Thyd.h>
#include <Domaine_VF.h>
#include <TRUSTTrav.h>
#include <EChaine.h>
#include <Domaine.h>
#include <Avanc.h>
#include <Param.h>
#include <Debog.h>
#include <SETS.h>

Implemente_instanciable(Energie_cinetique_turbulente,"Energie_cinetique_turbulente",Convection_diffusion_turbulence_multiphase);

Sortie& Energie_cinetique_turbulente::printOn(Sortie& is) const
{
  return Convection_diffusion_turbulence_multiphase::printOn(is);
}

Entree& Energie_cinetique_turbulente::readOn(Entree& is)
{
  Convection_diffusion_turbulence_multiphase::readOn(is);
  terme_convectif.set_fichier("Convection_energie_cinetique_turbulente");
  terme_convectif.set_description((Nom)"Turbulent kinetic energy transfer rate=Integral(-rho*k*ndS) [W] if SI units used");
  terme_diffusif.set_fichier("Diffusion_energie_cinetique_turbulente");
  terme_diffusif.set_description((Nom)"Turbulent kinetic energy transfer rate=Integral(mu*grad(k)*ndS) [W] if SI units used");
  return is;
}
void Energie_cinetique_turbulente::set_param(Param& param)
{
  Convection_diffusion_turbulence_multiphase::set_param(param);
  param.ajouter("limit_coef",&coef_limit_); // X_D attr limit_coef flottant limit_coef 1 Coefficient of the limiter (min (K, coef * v^2)). Default value of coef is set to 0.1
  param.ajouter_flag("limit_K",&limit_k_); // X_D attr limit_K entier limit_K 1 Flag to activate the limiter on K. Default value is 0 (deactivated)
}

const Champ_Don& Energie_cinetique_turbulente::diffusivite_pour_transport() const
{
  return ref_cast(Fluide_base,milieu()).viscosite_cinematique();
}

const Champ_base& Energie_cinetique_turbulente::diffusivite_pour_pas_de_temps() const
{
  return ref_cast(Fluide_base,milieu()).viscosite_cinematique();
}

void Energie_cinetique_turbulente::discretiser()
{
  int nb_valeurs_temp = schema_temps().nb_valeurs_temporelles();
  double temps = schema_temps().temps_courant();
  const Discret_Thyd& dis=ref_cast(Discret_Thyd, discretisation());
  Cerr << "Turbulent kinetic energy discretization" << finl;
  //On utilise temperature pour la directive car discretisation identique
  dis.discretiser_champ("temperature",domaine_dis(),"k","J/kg", 1,nb_valeurs_temp,temps,l_inco_ch);//une seule compo, meme en multiphase
  l_inco_ch->fixer_nature_du_champ(scalaire);
  l_inco_ch->fixer_nom_compo(0, Nom("k"));
  champs_compris_.ajoute_champ(l_inco_ch);
  Equation_base::discretiser();
  Cerr << "Energie_cinetique_turbulente::discretiser() ok" << finl;
}

void Energie_cinetique_turbulente::mettre_a_jour(double temps)
{
  // XXX : appel a la classe mere
  Convection_diffusion_turbulence_multiphase::mettre_a_jour(temps);

  const Navier_Stokes_std& eqv = ref_cast(Navier_Stokes_std, probleme().equation(0));

  if (probleme().discretisation().is_polymac_p0() && limit_k_ == 1)
    if ( temps > schema_temps().temps_courant() && coef_limit_ > 0 )
      {
        Cerr << "Limiting the value of K : coeff used = " << coef_limit_ << finl;
        const Champ_Face_PolyMAC_P0& ch_vit = ref_cast(Champ_Face_PolyMAC_P0, eqv.vitesse().valeur());
        const Domaine_PolyMAC_P0& domaine = ref_cast(Domaine_PolyMAC_P0, domaine_dis().valeur());
        DoubleTab& k_val = inconnue()->valeurs();
        const int N = k_val.line_size(), D = dimension;

        for (int e = 0; e < domaine.nb_elem(); e++)
          for (int n = 0; n < N; n++)
            {
              double norm_v2 = 0;
              for (int d = 0 ; d<D ; d++) norm_v2 += ch_vit.passe()(e, N*d+n);
              k_val(e, n) = std::min(k_val(e, n), coef_limit_ * norm_v2 );
            }
        k_val.echange_espace_virtuel();
        inconnue()->passe() = k_val;
      }
}

void Energie_cinetique_turbulente::calculer_alpha_rho_k(const Objet_U& obj, DoubleTab& val, DoubleTab& bval, tabs_t& deriv)
{
  const Equation_base& eqn = ref_cast(Equation_base, obj);

  /*  const Fluide_base& fl = ref_cast(Fluide_base, eqn.milieu());
    const Champ_base& ch_rho = fl.masse_volumique();
    const Champ_Inc_base *ch_alpha = sub_type(Pb_Multiphase, eqn.probleme()) ? &ref_cast(Pb_Multiphase, eqn.probleme()).equation_masse().inconnue().valeur() : nullptr,
                          *pch_rho = sub_type(Champ_Inc_base, ch_rho) ? &ref_cast(Champ_Inc_base, ch_rho) : nullptr; //pas toujours un Champ_Inc
    const DoubleTab* alpha = ch_alpha ? &ch_alpha->valeurs() : nullptr, &rho = ch_rho.valeurs(), &k = eqn.inconnue()->valeurs();
  */
  const DoubleTab& k = eqn.inconnue()->valeurs();

  /* valeurs du champ */
  int i, n, N = val.line_size(), Nl = val.dimension_tot(0);
  for (i = 0; i < Nl; i++)
    for (n = 0; n < N; n++) val(i, n) = k(i, n);

  /* on ne peut utiliser valeur_aux_bords que si ch_rho a un domaine_dis_base */
  const DoubleTab& b_k = eqn.inconnue()->valeur_aux_bords();
  int Nb = b_k.dimension_tot(0);
  for (i = 0; i < Nb; i++)
    for (n = 0; n < N; n++) bval(i, n) = b_k(i, n);

  //derivee en k : 1.
  DoubleTab& d_k = deriv["k"];
  for (d_k.resize(Nl, N), i = 0; i < Nl; i++)
    for (n = 0; n < N; n++) d_k(i, n) = 1.;
}
