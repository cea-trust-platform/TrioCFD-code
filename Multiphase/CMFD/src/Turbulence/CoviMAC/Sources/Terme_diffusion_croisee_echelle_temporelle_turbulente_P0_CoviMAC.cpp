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
// File:        Travail_pression_CoviMAC.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/CoviMAC/Sources
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Terme_diffusion_croisee_echelle_temporelle_turbulente_P0_CoviMAC.h>
#include <Zone_CoviMAC.h>
#include <Zone_Cl_CoviMAC.h>
#include <Champ_P0_CoviMAC.h>
#include <Equation_base.h>
#include <Pb_Multiphase.h>
#include <Milieu_composite.h>
#include <grad_Champ_Face_CoviMAC.h>
#include <Matrix_tools.h>
#include <Array_tools.h>

Implemente_instanciable(Terme_diffusion_croisee_echelle_temporelle_turbulente_P0_CoviMAC,"Terme_diffusion_croisee_echelle_temporelle_turbulente_P0_CoviMAC", Source_base);
// XD Terme_diffusion_croisee_echelle_temporelle_turbulente_P0_CoviMAC source_base Terme_diffusion_croisee_echelle_temporelle_turbulente 0 Source term which corresponds to the cross-diffusion term that appears in the turbulent tau equation (k-tau turbulence model)

Sortie& Terme_diffusion_croisee_echelle_temporelle_turbulente_P0_CoviMAC::printOn(Sortie& os) const
{
  return os;
}

Entree& Terme_diffusion_croisee_echelle_temporelle_turbulente_P0_CoviMAC::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("sigma_d", &sigma_d);
  param.lire_avec_accolades_depuis(is);
  return is;
}

void Terme_diffusion_croisee_echelle_temporelle_turbulente_P0_CoviMAC::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  const Zone_CoviMAC& 		zone 		= ref_cast(Zone_CoviMAC, equation().zone_dis().valeur());
  const Champ_P0_CoviMAC& 	ch_k 		= ref_cast(Champ_P0_CoviMAC, equation().probleme().get_champ("k"));	// Champ k
  const DoubleTab& 			k 			= ch_k.valeurs();	// ou valeurs passees ?
  const Champ_P0_CoviMAC& 	ch_tau 		= ref_cast(Champ_P0_CoviMAC, equation().inconnue().valeur()); 	// Champ tau
  const DoubleTab& 			tau 		= ch_tau.valeurs();
  const int 				N 			= tau.line_size(), nb_elem = zone.nb_elem();
  int e, n;

  if (!matrices.count(k.le_nom().getString()) || !matrices.count(tau.le_nom().getString())) return;

  for (auto &&n_m : matrices) if (n_m.first == "alpha" || n_m.first == "tau" || n_m.first == "temperature" || n_m.first == "pression")
      {
        Matrice_Morse& mat = *n_m.second, mat2;
        const DoubleTab& dep = equation().probleme().get_champ(n_m.first.c_str()).valeurs();
        int m,
            nc = dep.dimension_tot(0),	// nombre d'elements total
            M  = dep.line_size();		// nombre de composantes
        IntTrav sten(0, 2);
        sten.set_smart_resize(1);
        if (n_m.first == "alpha" || n_m.first == "temperature")	// N <= M
          for (e = 0; e < nb_elem; e++) for (n = 0; n < N; n++) sten.append_line(N * e + n, M * e + n); // il faudrait en fait faire la boucle sur la liste des id_composites des phases avec une equation tau
        if (n_m.first == "pression" || n_m.first == "tau" )		// M <= N
          for (e = 0; e < nb_elem; e++) for (n = 0, m = 0; n < N; n++, m+=(M>1)) sten.append_line(N * e + n, M * e + m); // il faudrait en fait faire la boucle sur la liste des id_composites des phases avec une equation tau
        Matrix_tools::allocate_morse_matrix(N * zone.nb_elem_tot(), M * nc, sten, mat2);
        mat.nb_colonnes() ? mat += mat2 : mat = mat2;
      }
}

// la phase dont la turbulence est decrite avec le modele k-tau doit etre ecrite en premier dans le bloc phases { } du jeu de donnees
// si plusieurs phases sont turbulentes et sont decrites par le modele k-tau, alors elles doivent se suivre dans le bloc phases { } du jeu de donnees
void Terme_diffusion_croisee_echelle_temporelle_turbulente_P0_CoviMAC::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Champ_P0_CoviMAC& 	ch_k 				= ref_cast(Champ_P0_CoviMAC, equation().probleme().get_champ("k"));	// Champ k
  const DoubleTab& 			k_passe				= ch_k.passe();
  const Champ_P0_CoviMAC& 	ch_tau 				= ref_cast(Champ_P0_CoviMAC, equation().inconnue().valeur()); 		// Champ tau
  const DoubleTab& 			tau_passe			= ch_tau.passe();
  const Champ_Inc_base& 	ch_alpha_rho_tau 	= equation().champ_conserve();
  const DoubleTab& 			alpha_rho_tau		= ch_alpha_rho_tau.valeurs();
  const tabs_t& 			der_alpha_rho_tau 	= ch_alpha_rho_tau.derivees(); // dictionnaire des derivees
  const Zone_CoviMAC& 		zone 				= ref_cast(Zone_CoviMAC, equation().zone_dis().valeur());
  const Conds_lim& 		   	cls_tau				= ch_tau.zone_Cl_dis().les_conditions_limites(); 	// conditions aux limites du champ tau
  const Conds_lim&          cls_k 				= ch_k.zone_Cl_dis().les_conditions_limites(); 		// conditions aux limites du champ k
  const IntTab&				fcl_tau 			= ch_tau.fcl(); // tableaux utilitaires sur les CLs : fcl(f, .) = (type de la CL, no de la CL, indice dans la CL)
  const IntTab&             fcl_k 				= ch_k.fcl();	// tableaux utilitaires sur les CLs : fcl(f, .) = (type de la CL, no de la CL, indice dans la CL)

  const int Mt = tau_passe.dimension(1), nf = zone.nb_faces(), D = dimension, nb_elem = zone.nb_elem() ;
  const int Ntau = tau_passe.line_size(), Np = equation().probleme().get_champ("pression").valeurs().line_size(), Na = equation().probleme().get_champ("alpha").valeurs().line_size(), Nt = equation().probleme().get_champ("temperature").valeurs().line_size();

  /* Calcul de grad de tau aux faces */

  DoubleTrav grad_f_tau(0, Mt);

  MD_Vector_tools::creer_tableau_distribue(equation().probleme().get_champ("vitesse").valeurs().get_md_vector(), grad_f_tau);
  ch_tau.init_grad(0);
  IntTab& f_d_tau = ch_tau.fgrad_d, f_e_tau = ch_tau.fgrad_e; // Tables utilisees dans zone_CoviMAC::fgrad
  DoubleTab f_w_tau = ch_tau.fgrad_w;

  for (int n = 0; n < Ntau; n++) for (int f = 0; f < nf; f++)
      {
        grad_f_tau(f, n) = 0;
        for (int j = f_d_tau(f); j < f_d_tau(f+1) ; j++)
          {
            int e = f_e_tau(j);
            int f_bord;
            if (e < nb_elem) //contrib d'un element
              {
                double val_e = tau_passe(e, n);
                grad_f_tau(f, n) += f_w_tau(j) * val_e;
              }
            else if (fcl_tau(f_bord = e - nb_elem, 0) == 3) //contrib d'un bord : seul Dirichlet contribue
              {
                double val_f_bord = sqrt(ref_cast(Dirichlet, cls_tau[fcl_tau(f_bord, 1)].valeur()).val_imp(fcl_tau(f_bord, 2), n));
                grad_f_tau(f, n) += f_w_tau(j) * val_f_bord;
              }
          }
      }

  /* Calcul de grad de k aux faces */

  DoubleTrav grad_f_k(0, Mt);
  MD_Vector_tools::creer_tableau_distribue(equation().probleme().get_champ("vitesse").valeurs().get_md_vector(), grad_f_k); // on cree un tableau distribue de meme structure que la vitesse
  ch_k.init_grad(0);
  IntTab& f_d_k = ch_k.fgrad_d, f_e_k = ch_k.fgrad_e;  // Tables utilisees dans zone_CoviMAC::fgrad
  DoubleTab f_w_k = ch_k.fgrad_w;

  for (int n = 0; n < Ntau; n++) for (int f = 0; f < nf; f++)
      {
        grad_f_k(f, n) = 0;
        for (int j = f_d_k(f); j < f_d_k(f+1) ; j++)
          {
            int e = f_e_k(j);
            int f_bord;
            if (e < nb_elem) //contrib d'un element
              {
                double val_e = k_passe(e, n);
                grad_f_k(f, n) += f_w_k(j) * val_e;
              }
            else if (fcl_k(f_bord = e - nb_elem, 0) == 3) //contrib d'un bord : seul Dirichlet contribue
              {
                double val_f_bord = sqrt(ref_cast(Dirichlet, cls_k[fcl_k(f_bord, 1)].valeur()).val_imp(fcl_k(f_bord, 2), n));
                grad_f_k(f, n) += f_w_k(j) * val_f_bord;
              }
          }
      }

  /* Calcul de grad(tau).(grad k) */

  DoubleTrav grad_f_tau_dot_grad_f_k(0, Mt);
  MD_Vector_tools::creer_tableau_distribue(equation().probleme().get_champ("vitesse").valeurs().get_md_vector(), grad_f_tau_dot_grad_f_k); // on cree un tableau distribue de meme structure que la vitesse
  zone.init_ve();
  for (int n = 0; n < Ntau; n++) for (int e = 0; e < nb_elem; e++)
      {
        grad_f_tau_dot_grad_f_k(e, n) = 0;
        std::vector<double> grad_tau(D), grad_k(D);
        for (int j = zone.ved(e); j < zone.ved(e + 1); j++) for (int f = zone.vej(j), d = 0; d < D; d++)
            {
              grad_tau[d] += zone.vec(j, d) * grad_f_tau(f, n);
              grad_k[d] 	+= zone.vec(j, d) * grad_f_k(f, n);
            }
        for (int d = 0 ; d<D ; d++) grad_f_tau_dot_grad_f_k(e, n) += grad_tau[d] * grad_k[d]; // produit scalaire
      }

  /* remplissage des matrices et du second membre */

  Matrice_Morse *Mtau = matrices.count(ch_tau.le_nom().getString()) ? matrices.at(ch_tau.le_nom().getString()) : NULL;
  Matrice_Morse *Ma = matrices.count("alpha") ? matrices.at("alpha") : NULL;
  Matrice_Morse *Mp = matrices.count("pression") ? matrices.at("pression") : NULL;
  Matrice_Morse *Mtemp	= matrices.count("temperature") ? matrices.at("temperature") : NULL;

  for (int e = 0; e < nb_elem; e++) for (int ntau = 0, mp = 0; ntau < Ntau; ntau++, mp += (Np>1))
      {
        secmem(e, ntau) += sigma_d * alpha_rho_tau(e, ntau) * std::min(grad_f_tau_dot_grad_f_k(e, ntau), 0.);
        if (Ma)	  (*Ma)(Ntau * e + ntau, Na * e + ntau)   -= sigma_d * (der_alpha_rho_tau.count("alpha") ? der_alpha_rho_tau.at("alpha")(e,ntau) : NULL ) * std::min(grad_f_tau_dot_grad_f_k(e, ntau), 0.); // derivee en alpha
        if (Mtemp)(*Mtemp)(Ntau * e + ntau, Nt * e + ntau)-= sigma_d * (der_alpha_rho_tau.count("temperature") ? der_alpha_rho_tau.at("temperature")(e,ntau) : NULL ) * std::min(grad_f_tau_dot_grad_f_k(e, ntau), 0.); // derivee par rapport a la temperature
        if (Mp)	  (*Mp)(Ntau * e + ntau, Np * e + mp)   -= sigma_d * (der_alpha_rho_tau.count("pression") ? der_alpha_rho_tau.at("pression")(e,ntau) : NULL ) * std::min(grad_f_tau_dot_grad_f_k(e, ntau), 0.); // derivee par rapport a la pression
        if (Mtau) (*Mtau)(Ntau * e + ntau, Ntau * e + ntau) -= sigma_d * (der_alpha_rho_tau.count("tau") ? der_alpha_rho_tau.at("tau")(e,ntau) : NULL ) * std::min(grad_f_tau_dot_grad_f_k(e, ntau), 0.); // derivee en tau
      }
}


