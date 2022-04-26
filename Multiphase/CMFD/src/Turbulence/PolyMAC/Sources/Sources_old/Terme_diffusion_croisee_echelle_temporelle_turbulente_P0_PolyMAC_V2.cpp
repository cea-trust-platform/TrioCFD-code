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
// File:        Travail_pression_PolyMAC_V2.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/PolyMAC_V2/Sources
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Terme_diffusion_croisee_echelle_temporelle_turbulente_P0_PolyMAC_V2.h>
#include <Zone_PolyMAC_V2.h>
#include <Zone_Cl_PolyMAC.h>
#include <Champ_P0_PolyMAC_V2.h>
#include <Equation_base.h>
#include <Probleme_base.h>
#include <Milieu_composite.h>
#include <grad_Champ_Face_PolyMAC_V2.h>
#include <Matrix_tools.h>
#include <Array_tools.h>

Implemente_instanciable(Terme_diffusion_croisee_echelle_temporelle_turbulente_P0_PolyMAC_V2,"Terme_diffusion_croisee_echelle_temporelle_turbulente_P0_PolyMAC_V2", Source_base);
// XD Terme_diffusion_croisee_echelle_temporelle_turbulente_P0_PolyMAC_V2 source_base Terme_diffusion_croisee_echelle_temporelle_turbulente 0 Source term which corresponds to the cross-diffusion term that appears in the turbulent tau equation (k-tau turbulence model)

Sortie& Terme_diffusion_croisee_echelle_temporelle_turbulente_P0_PolyMAC_V2::printOn(Sortie& os) const
{
  return os;
}

Entree& Terme_diffusion_croisee_echelle_temporelle_turbulente_P0_PolyMAC_V2::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("sigma_d", &sigma_d);
  param.lire_avec_accolades_depuis(is);
  return is;
}

void Terme_diffusion_croisee_echelle_temporelle_turbulente_P0_PolyMAC_V2::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  const Zone_PolyMAC_V2& 		zone 		= ref_cast(Zone_PolyMAC_V2, equation().zone_dis().valeur());
  const Champ_P0_PolyMAC_V2& 	ch_tau 		= ref_cast(Champ_P0_PolyMAC_V2, equation().inconnue().valeur()); 		// Champ tau
  const DoubleTab& 			tau 		= ch_tau.valeurs();
  const int N = tau.line_size(), nb_elem = zone.nb_elem();
  int e, n;

  assert(N == 1); // si N > 1 il vaut mieux iterer sur les id_composites des phases turbulentes
  assert(ref_cast(Champ_P0_PolyMAC_V2, equation().probleme().get_champ("k")).valeurs().line_size() == 1);

  for (auto &&n_m : matrices) if (n_m.first == "alpha" || n_m.first == "temperature" || n_m.first == "pression")
      {
        Matrice_Morse& mat = *n_m.second, mat2;
        const DoubleTab& dep = equation().probleme().get_champ(n_m.first.c_str()).valeurs();
        int m,
            nc = dep.dimension_tot(0),	// nombre d'elements total
            M  = dep.line_size();		// nombre de composantes
        IntTrav sten(0, 2);
        sten.set_smart_resize(1);
        if (n_m.first == "alpha" || n_m.first == "temperature")	// N <= M
          for (e = 0; e < nb_elem; e++) for (n = 0; n < N; n++) sten.append_line(N * e + n, M * e + n);
        if (n_m.first == "pression")
          for (e = 0; e < nb_elem; e++) for (n = 0; n < N; n++) for (m = 0; m<M; m++) sten.append_line(N * e + n, M * e + m);
        Matrix_tools::allocate_morse_matrix(N * zone.nb_elem_tot(), M * nc, sten, mat2);
        mat.nb_colonnes() ? mat += mat2 : mat = mat2;
      }
}

void Terme_diffusion_croisee_echelle_temporelle_turbulente_P0_PolyMAC_V2::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Champ_P0_PolyMAC_V2& 	ch_k 				= ref_cast(Champ_P0_PolyMAC_V2, equation().probleme().get_champ("k"));	// Champ k
  const DoubleTab& 			k_passe				= ch_k.passe();
  const Champ_P0_PolyMAC_V2& 	ch_tau 				= ref_cast(Champ_P0_PolyMAC_V2, equation().inconnue().valeur()); 		// Champ tau
  const DoubleTab& 			tau_passe			= ch_tau.passe();
  const Champ_Inc_base& 	ch_alpha_rho_tau 	= equation().champ_conserve();
  const DoubleTab& 			alpha_rho_tau		= ch_alpha_rho_tau.valeurs();
  const tabs_t& 			der_alpha_rho_tau 	= ch_alpha_rho_tau.derivees(); // dictionnaire des derivees
  const Zone_PolyMAC_V2& 		zone 				= ref_cast(Zone_PolyMAC_V2, equation().zone_dis().valeur());
  const Conds_lim& 		   	cls_tau				= ch_tau.zone_Cl_dis().les_conditions_limites(); 	// conditions aux limites du champ tau
  const Conds_lim&          cls_k 				= ch_k.zone_Cl_dis().les_conditions_limites(); 		// conditions aux limites du champ k
  const IntTab&				fcl_tau 			= ch_tau.fcl(); // tableaux utilitaires sur les CLs : fcl(f, .) = (type de la CL, no de la CL, indice dans la CL)
  const IntTab&             fcl_k 				= ch_k.fcl(), &e_f = zone.elem_faces(), &f_e = zone.face_voisins();	// tableaux utilitaires sur les CLs : fcl(f, .) = (type de la CL, no de la CL, indice dans la CL)

  const int Mt = tau_passe.dimension(1), nf = zone.nb_faces(), D = dimension, nb_elem = zone.nb_elem(), nb_elem_tot = zone.nb_elem() ;
  const int Ntau = tau_passe.line_size(), Np = equation().probleme().get_champ("pression").valeurs().line_size(), Na = equation().probleme().get_champ("alpha").valeurs().line_size(), Nt = equation().probleme().get_champ("temperature").valeurs().line_size();

  assert(Ntau == 1 || k_passe.line_size() == 1); // si Ntau > 1 il vaut mieux iterer sur les id_composites des phases turbulentes decrites par un modele k-tau dans le calcul de grad_f_k et dans le remplissage des matrices

  /* Calcul de grad de tau aux faces */

  DoubleTrav grad_f_tau(nf, Mt);
  ch_tau.init_grad(0); // Initialisation des tables fgrad_d, fgrad_e, fgrad_w qui dependent de la discretisation et du type de conditions aux limites --> pas de mises a jour necessaires
  const IntTab& f_d_tau = ch_tau.fgrad_d, &f_e_tau = ch_tau.fgrad_e; // Tables utilisees dans zone_PolyMAC_V2::fgrad pour le calcul du gradient
  const DoubleTab& f_w_tau = ch_tau.fgrad_w, &xp = zone.xp(), &xv = zone.xv();
  const DoubleVect& fs = zone.face_surfaces(), &ve = zone.volumes();

  for (int n = 0; n < Ntau; n++) for (int f = 0; f < nf; f++)
      {
        grad_f_tau(f, n) = 0;
        for (int j = f_d_tau(f); j < f_d_tau(f+1) ; j++)
          {
            int e = f_e_tau(j);
            int f_bord;
            if (e < nb_elem_tot) //contribution d'un element
              {
                double val_e = tau_passe(e, n);
                grad_f_tau(f, n) += f_w_tau(j) * val_e;
              }
            else if (fcl_tau(f_bord = e - nb_elem_tot, 0) == 6) //contribution d'un bord : seul Dirichlet contribue
              {
                double val_f_bord = ref_cast(Dirichlet, cls_tau[fcl_tau(f_bord, 1)].valeur()).val_imp(fcl_tau(f_bord, 2), n);
                grad_f_tau(f, n) += f_w_tau(j) * val_f_bord;
              }
          }
      }

  /* Calcul de grad de k aux faces */

  DoubleTrav grad_f_k(nf, Mt);
  ch_k.init_grad(0); // Initialisation des tables fgrad_d, fgrad_e, fgrad_w qui dependent de la discretisation et du type de conditions aux limites --> pas de mises a jour necessaires
  const IntTab& f_d_k = ch_k.fgrad_d, &f_e_k = ch_k.fgrad_e;  // Tables utilisees dans zone_PolyMAC_V2::fgrad pour le calcul du gradient
  const DoubleTab& f_w_k = ch_k.fgrad_w;

  for (int n = 0; n < Ntau; n++) for (int f = 0; f < nf; f++)
      {
        grad_f_k(f, n) = 0;
        for (int j = f_d_k(f); j < f_d_k(f+1) ; j++)
          {
            int e = f_e_k(j);
            int f_bord;
            if (e < nb_elem_tot) //contribution d'un element
              {
                double val_e = k_passe(e, n);
                grad_f_k(f, n) += f_w_k(j) * val_e;
              }
            else if (fcl_k(f_bord = e - nb_elem_tot, 0) == 6) //contribution d'un bord : seul Dirichlet contribue
              {
                double val_f_bord = ref_cast(Dirichlet, cls_k[fcl_k(f_bord, 1)].valeur()).val_imp(fcl_k(f_bord, 2), n);
                grad_f_k(f, n) += f_w_k(j) * val_f_bord;
              }
          }
      }

  /* Calcul de grad(tau).(grad k) */

  DoubleTrav grad_f_tau_dot_grad_f_k(0, Mt);
  MD_Vector_tools::creer_tableau_distribue(equation().probleme().get_champ("vitesse").valeurs().get_md_vector(), grad_f_tau_dot_grad_f_k); // on cree un tableau distribue de meme structure que la vitesse
  for (int n = 0; n < Ntau; n++) for (int e = 0; e < nb_elem; e++)
      {
        grad_f_tau_dot_grad_f_k(e, n) = 0;
        std::vector<double> grad_tau(D), grad_k(D);
        for (int j = 0, f; j < e_f.dimension(1) && (f = e_f(e, j)) >= 0; j++) for (int d = 0; d < D; d++)
            {
              grad_tau[d] += (e == f_e(f, 0) ? 1 : -1) * fs(f) * (xv(f, d) - xp(e, d)) / ve(e) * grad_f_tau(f, n);
              grad_k[d]   += (e == f_e(f, 0) ? 1 : -1) * fs(f) * (xv(f, d) - xp(e, d)) / ve(e) * grad_f_k(f, n);
            }
        for (int d = 0 ; d < D ; d++) grad_f_tau_dot_grad_f_k(e, n) += grad_tau[d] * grad_k[d]; // produit scalaire
      }

  /* remplissage des matrices et du second membre */

  Matrice_Morse *Mtau = matrices.count(ch_tau.le_nom().getString()) ? matrices.at(ch_tau.le_nom().getString()) : nullptr;
  Matrice_Morse *Ma = matrices.count("alpha") ? matrices.at("alpha") : nullptr;
  Matrice_Morse *Mp = matrices.count("pression") ? matrices.at("pression") : nullptr;
  Matrice_Morse *Mtemp	= matrices.count("temperature") ? matrices.at("temperature") : nullptr;

  for (int e = 0; e < nb_elem; e++) for (int ntau = 0; ntau < Ntau; ntau++)
      {
        secmem(e, ntau) += sigma_d * alpha_rho_tau(e, ntau) * std::min(grad_f_tau_dot_grad_f_k(e, ntau), 0.);
        if (!(Ma==nullptr))    (*Ma)(Ntau * e + ntau, Na * e + ntau)   	-= sigma_d * (der_alpha_rho_tau.count("alpha") ? der_alpha_rho_tau.at("alpha")(e,ntau) : 0 ) * std::min(grad_f_tau_dot_grad_f_k(e, ntau), 0.); // derivee en alpha
        if (!(Mtemp==nullptr)) (*Mtemp)(Ntau * e + ntau, Nt * e + ntau)	-= sigma_d * (der_alpha_rho_tau.count("temperature") ? der_alpha_rho_tau.at("temperature")(e,ntau) : 0 ) * std::min(grad_f_tau_dot_grad_f_k(e, ntau), 0.); // derivee par rapport a la temperature
        if (!(Mtau==nullptr))  (*Mtau)(Ntau * e + ntau, Ntau * e + ntau)-= sigma_d * (der_alpha_rho_tau.count("tau") ? der_alpha_rho_tau.at("tau")(e,ntau) : 0 ) * std::min(grad_f_tau_dot_grad_f_k(e, ntau), 0.); // derivee en tau
        for (int mp = 0; mp<Np; mp++) if (!(Mp==nullptr))
            (*Mp)(Ntau * e + ntau, Np * e + mp)     	-= sigma_d * (der_alpha_rho_tau.count("pression") ? der_alpha_rho_tau.at("pression")(e,ntau) : 0 ) * std::min(grad_f_tau_dot_grad_f_k(e, ntau), 0.); // derivee par rapport a la pression

      }
}


