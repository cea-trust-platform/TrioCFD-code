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
// File:        Transport_turbulent_GGDH_WIT.cpp
// Directory:   $TRUST_ROOT/src/Turbulence/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Transport_turbulent_GGDH_WIT.h>
#include <Param.h>
#include <Probleme_base.h>

#include <Pb_Multiphase.h>
#include <TRUSTTrav.h>
#include <MD_Vector_tools.h>
#include <math.h>
#include <TRUSTTab_parts.h>
#include <cmath>
#include <Viscosite_turbulente_multiple.h>

Implemente_instanciable(Transport_turbulent_GGDH_WIT, "Transport_turbulent_GGDH_WIT", Transport_turbulent_base);

Sortie& Transport_turbulent_GGDH_WIT::printOn(Sortie& os) const
{
  return os;
}

Entree& Transport_turbulent_GGDH_WIT::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("Aspect_ratio", &gamma_); // rapport d'aspet des bulles
  param.ajouter("Influence_area", &delta_); // paramètre modèle d'Alméras 2014 (taille du sillage)
  param.ajouter("C_s", &C_s); // paramètre modèle d'Alméras 2014
  //param.ajouter("vitesse_rel_attendue", &ur_user, Param::REQUIRED); // valeur de ur à prendre si u_r(i,0)=0
  param.ajouter("Limiteur_alpha", &limiteur_alpha_, Param::REQUIRED); // valeur minimal de (1-alpha) pour utiliser le modèle d'Alméras
  param.lire_avec_accolades_depuis(is);

  if ((C_s <0) && (dimension == 2)) C_s = 1;
  if ((C_s <0) && (dimension == 3)) C_s = 1.5;

  return is;
}

void Transport_turbulent_GGDH_WIT::modifier_mu(const Convection_Diffusion_std& eq, const Viscosite_turbulente_base& visc_turb, DoubleTab& nu) const
{
  const DoubleTab& mu0 = eq.diffusivite_pour_transport().passe(), &nu0 = eq.diffusivite_pour_pas_de_temps().passe(), //viscosites moleculaires
                   alp = pb_->get_champ("alpha").passe(), diam = pb_->get_champ("diametre_bulles").valeurs(),
                   &tab_u = pb_->get_champ("vitesse").passe();
  int i, nl = nu.dimension(0), N = alp.dimension(1), d, db, D = dimension, i_part=-1;
  //if (N!= 1) Process::exit("Only the liquid phase can have WIT for now");
  DoubleTrav Rij(0, N, D, D);
  MD_Vector_tools::creer_tableau_distribue(nu.get_md_vector(), Rij);

  ConstDoubleTab_parts p_u(tab_u); //en PolyMAC_P0, tab_u contient (nf.u) aux faces, puis (u_i) aux elements
  for (i = 0; i < p_u.size(); i++)
    if (p_u[i].get_md_vector() == Rij.get_md_vector()) i_part = i; //on cherche une partie ayant le meme support
  if (i_part < 0) Process::exit("Viscosite_turbulente_WIF : inconsistency between velocity and Rij!");
  const DoubleTab& u = p_u[i_part]; //le bon tableau
  DoubleTrav u_r(u.dimension(0), 1);
  for (i = 0; i < u_r.dimension(0); i++)
    {
      for (d = 0; d < D; d++) u_r(i, 0) += (u(i, d, 1) - u(i, d, 0))*(u(i, d, 1) - u(i, d, 0)); // relative speed = gas speed - liquid speed
      u_r(i, 0) = std::sqrt(u_r(i, 0));
    }

  visc_turb.reynolds_stress(Rij);
  //formule pour passer de nu a mu : mu0 / nu0 * C_s * temps_carac * <u'i u'_j>
  for (i = 0; i < nl; i++)
    if (alp(i, 0) >= limiteur_alpha_)
      for (d = 0; d < D; d++)
        for (db = 0; db < D; db++)
          {
            //double ur = (u_r(i,0)!=0) ? u_r(i,0) : ur_user; // vitesse relative (on ne peut pas diviser par 0)
            double temps_carac = (u_r(i,0)!=0) ? 2./3. * 1./(delta_*delta_*delta_)*diam(i, 1) / (pow(gamma_, 2./3.)*alp(i, 1)*u_r(i,0)) : 0 ; // si u_r=0 alors pas de WIT donc pas de diffusion du WIT
            nu(i, 0, d, db) += mu0(i, 0) / nu0(i, 0) * C_s * std::min(std::max(temps_carac * Rij(i, 0, d, db), visc_turb.min_ev_ratio() * nu(i, 0, d, db)), visc_turb.max_ev_ratio() * nu(i, 0, d, db));
          }
}
