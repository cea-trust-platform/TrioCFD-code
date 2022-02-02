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
// File:        Diffusion_croisee_echelle_temp_taux_diss_turb_CoviMAC.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/CoviMAC/Sources
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Diffusion_croisee_echelle_temp_taux_diss_turb_CoviMAC.h>
#include <Zone_CoviMAC.h>
#include <Zone_Cl_CoviMAC.h>
#include <Champ_P0_CoviMAC.h>
#include <Equation_base.h>
#include <Probleme_base.h>
#include <Pb_Multiphase.h>
#include <Milieu_composite.h>
#include <grad_Champ_Face_CoviMAC.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <Echelle_temporelle_turbulente.h>
#include <Taux_dissipation_turbulent.h>

Implemente_instanciable(Diffusion_croisee_echelle_temp_taux_diss_turb_CoviMAC,"Diffusion_croisee_echelle_temp_taux_diss_turb_P0_CoviMAC", Source_base);

Sortie& Diffusion_croisee_echelle_temp_taux_diss_turb_CoviMAC::printOn(Sortie& os) const
{
  return os;
}

Entree& Diffusion_croisee_echelle_temp_taux_diss_turb_CoviMAC::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("sigma_d", &sigma_d);
  param.lire_avec_accolades_depuis(is);
  return is;
}

void Diffusion_croisee_echelle_temp_taux_diss_turb_CoviMAC::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  const Zone_CoviMAC& 		zone 		= ref_cast(Zone_CoviMAC, equation().zone_dis().valeur());
  const Champ_P0_CoviMAC& 	ch_diss 		= ref_cast(Champ_P0_CoviMAC, equation().inconnue().valeur()); 		// Champ tau
  const int N = ch_diss.valeurs().line_size(), nb_elem = zone.nb_elem();
  int e, n;

  assert(N == 1); // si N > 1 il vaut mieux iterer sur les id_composites des phases turbulentes
  assert(ref_cast(Champ_P0_CoviMAC, equation().probleme().get_champ("k")).valeurs().line_size() == 1);

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

void Diffusion_croisee_echelle_temp_taux_diss_turb_CoviMAC::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Zone_CoviMAC& 		zone 				= ref_cast(Zone_CoviMAC, equation().zone_dis().valeur());
  const Champ_P0_CoviMAC& 	ch_k 				= ref_cast(Champ_P0_CoviMAC, equation().probleme().get_champ("k"));	// Champ k
  const DoubleTab& 	  		k_passe				= ch_k.passe();
  const Conds_lim&          cls_k 				= ch_k.zone_Cl_dis().les_conditions_limites(); 		// conditions aux limites du champ k
  const IntTab&             fcl_k 				= ch_k.fcl();	// tableaux utilitaires sur les CLs : fcl(f, .) = (type de la CL, no de la CL, indice dans la CL)

  const Champ_P0_CoviMAC& 	ch_diss 		= ref_cast(Champ_P0_CoviMAC, equation().inconnue().valeur()); 		// Champ tau ou omega
  const DoubleTab& 			diss_passe			= ch_diss.passe();
  const DoubleTab& 			      diss			= ch_diss.valeurs();
  const Conds_lim& 		   	cls_diss			= ch_diss.zone_Cl_dis().les_conditions_limites(); 	// conditions aux limites du champ tau ou omega
  const IntTab&				   fcl_diss 			= ch_diss.fcl(); // tableaux utilitaires sur les CLs : fcl(f, .) = (type de la CL, no de la CL, indice dans la CL)

  const int nf = zone.nb_faces(), D = dimension, nb_elem = zone.nb_elem(), nb_elem_tot = zone.nb_elem() ;
  const int N = diss_passe.line_size(), Np = equation().probleme().get_champ("pression").valeurs().line_size(), Na = equation().probleme().get_champ("alpha").valeurs().line_size(), Nt = equation().probleme().get_champ("temperature").valeurs().line_size();

  std::string Type_diss = ""; // omega or tau dissipation
  if sub_type(Echelle_temporelle_turbulente, equation()) Type_diss = "tau";
  else if sub_type(Taux_dissipation_turbulent, equation()) Type_diss = "omega";
  if (Type_diss == "") abort();


  assert(N == 1 || k_passe.line_size() == 1); // si Ntau > 1 il vaut mieux iterer sur les id_composites des phases turbulentes decrites par un modele k-tau dans le calcul de grad_f_k et dans le remplissage des matrices

  /* Calcul de grad de tau ou de omega aux faces */

  DoubleTrav grad_f_diss(nf, N);
  ch_diss.init_grad(0); // Initialisation des tables fgrad_d, fgrad_e, fgrad_w qui dependent de la discretisation et du type de conditions aux limites --> pas de mises a jour necessaires
  IntTab& f_d_tau = ch_diss.fgrad_d, f_e_tau = ch_diss.fgrad_e; // Tables utilisees dans zone_CoviMAC::fgrad pour le calcul du gradient
  DoubleTab f_w_tau = ch_diss.fgrad_w;

  for (int n = 0; n < N; n++) for (int f = 0; f < nf; f++)
      {
        grad_f_diss(f, n) = 0;
        for (int j = f_d_tau(f); j < f_d_tau(f+1) ; j++)
          {
            int e = f_e_tau(j);
            int f_bord;
            if (e < nb_elem_tot) //contribution d'un element
              {
                double val_e = diss_passe(e, n);
                grad_f_diss(f, n) += f_w_tau(j) * val_e;
              }
            else if (fcl_diss(f_bord = e - nb_elem_tot, 0) == 3) //contribution d'un bord : seul Dirichlet contribue
              {
                double val_f_bord = ref_cast(Dirichlet, cls_diss[fcl_diss(f_bord, 1)].valeur()).val_imp(fcl_diss(f_bord, 2), n);
                grad_f_diss(f, n) += f_w_tau(j) * val_f_bord;
              }
          }
      }

  /* Calcul de grad de k aux faces */

  DoubleTrav grad_f_k(nf, N);
  ch_k.init_grad(0); // Initialisation des tables fgrad_d, fgrad_e, fgrad_w qui dependent de la discretisation et du type de conditions aux limites --> pas de mises a jour necessaires
  IntTab& f_d_k = ch_k.fgrad_d, f_e_k = ch_k.fgrad_e;  // Tables utilisees dans zone_CoviMAC::fgrad pour le calcul du gradient
  DoubleTab f_w_k = ch_k.fgrad_w;

  for (int n = 0; n < N; n++) for (int f = 0; f < nf; f++)
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
            else if (fcl_k(f_bord = e - nb_elem_tot, 0) == 3) //contribution d'un bord : seul Dirichlet contribue
              {
                double val_f_bord = ref_cast(Dirichlet, cls_k[fcl_k(f_bord, 1)].valeur()).val_imp(fcl_k(f_bord, 2), n);
                grad_f_k(f, n) += f_w_k(j) * val_f_bord;
              }
          }
      }

  /* Calcul de grad(tau/omega).(grad k) */

  DoubleTrav grad_f_diss_dot_grad_f_k(nb_elem_tot, N);
  zone.init_ve();
  for (int n = 0; n < N; n++) for (int e = 0; e < nb_elem; e++)
      {
        grad_f_diss_dot_grad_f_k(e, n) = 0;
        std::vector<double> grad_diss(D), grad_k(D);
        for (int j = zone.ved(e); j < zone.ved(e + 1); j++) for (int f = zone.vej(j), d = 0; d < D; d++)
            {
              grad_diss[d] += zone.vec(j, d) * grad_f_diss(f, n);
              grad_k[d]   += zone.vec(j, d) * grad_f_k(f, n);
            }
        for (int d = 0 ; d < D ; d++) grad_f_diss_dot_grad_f_k(e, n) += grad_diss[d] * grad_k[d]; // produit scalaire
      }

  /* remplissage des matrices et du second membre */

  Matrice_Morse *M = matrices.count(ch_diss.le_nom().getString()) ? matrices.at(ch_diss.le_nom().getString()) : nullptr;
  Matrice_Morse *Ma = matrices.count("alpha") ? matrices.at("alpha") : nullptr;
  Matrice_Morse *Mp = matrices.count("pression") ? matrices.at("pression") : nullptr;
  Matrice_Morse *Mtemp	= matrices.count("temperature") ? matrices.at("temperature") : nullptr;

  for (int e = 0; e < nb_elem; e++) for (int n = 0; n < N; n++)
      {
        if (Type_diss == "tau")
          {
            const Champ_Inc_base& 	ch_alpha_rho_tau 	= equation().champ_conserve();
            const DoubleTab& 			alpha_rho_tau		= ch_alpha_rho_tau.valeurs();
            const tabs_t& 			der_alpha_rho_tau 	= ch_alpha_rho_tau.derivees(); // dictionnaire des derivees
            double secmem_en = sigma_d * alpha_rho_tau(e, n) * std::min(grad_f_diss_dot_grad_f_k(e, n), 0.);
            secmem(e, n) += secmem_en;
            if (!(Ma==nullptr))    (*Ma)(N * e + n, Na * e + n)   	-= sigma_d * (der_alpha_rho_tau.count("alpha") ? der_alpha_rho_tau.at("alpha")(e,n) : 0 ) * std::min(grad_f_diss_dot_grad_f_k(e, n), 0.); // derivee en alpha
            if (!(Mtemp==nullptr)) (*Mtemp)(N * e + n, Nt * e + n)	-= sigma_d * (der_alpha_rho_tau.count("temperature") ? der_alpha_rho_tau.at("temperature")(e,n) : 0 ) * std::min(grad_f_diss_dot_grad_f_k(e, n), 0.); // derivee par rapport a la temperature
            if (!(M==nullptr))     (*M)(N * e + n, N * e + n)       -= sigma_d * (der_alpha_rho_tau.count("tau") ? der_alpha_rho_tau.at("tau")(e,n) : 0 ) * std::min(grad_f_diss_dot_grad_f_k(e, n), 0.); // derivee en tau
            for (int mp = 0; mp<Np; mp++) if (!(Mp==nullptr))
                (*Mp)(N * e + n, Np * e + mp)     	-= sigma_d * (der_alpha_rho_tau.count("pression") ? der_alpha_rho_tau.at("pression")(e,n) : 0 ) * std::min(grad_f_diss_dot_grad_f_k(e, n), 0.); // derivee par rapport a la pression
          }
        else if (Type_diss == "omega")
          {
            const Champ_base& 		ch_alpha_rho 	= sub_type(Pb_Multiphase,equation().probleme()) ? ref_cast(Pb_Multiphase,equation().probleme()).eq_masse.champ_conserve() : equation().milieu().masse_volumique();
            const DoubleTab& 			alpha_rho		= ch_alpha_rho.valeurs();
            const tabs_t& 			der_alpha_rho 	= ref_cast(Champ_Inc_base, ch_alpha_rho).derivees(); // dictionnaire des derivees

            if (diss(e,n)>0) // Else everything = 0
              {
                secmem(e, n) += sigma_d * alpha_rho(e, n) / diss(e, n)* std::max(grad_f_diss_dot_grad_f_k(e, n), 0.) ;
                if (!(Ma==nullptr))    (*Ma)(N * e + n, Na * e + n)   	-= sigma_d * (der_alpha_rho.count("alpha") ? der_alpha_rho.at("alpha")(e,n) : 0 ) * std::max(grad_f_diss_dot_grad_f_k(e, n), 0.); // derivee en alpha
                if (!(Mtemp==nullptr)) (*Mtemp)(N * e + n, Nt * e + n)	-= sigma_d * (der_alpha_rho.count("temperature") ? der_alpha_rho.at("temperature")(e,n) : 0 ) * std::max(grad_f_diss_dot_grad_f_k(e, n), 0.); // derivee par rapport a la temperature
                for (int mp = 0; mp<Np; mp++) if (!(Mp==nullptr))
                    (*Mp)(N * e + n, Np * e + mp)     	-= sigma_d * (der_alpha_rho.count("pression") ? der_alpha_rho.at("pression")(e,n) : 0 ) * std::max(grad_f_diss_dot_grad_f_k(e, n), 0.); // derivee par rapport a la pression
                if (!(M==nullptr))     (*M)(N * e + n, N * e + n) -= sigma_d * alpha_rho(e, n) * (-1/(diss(e,n)*diss(e,n))) * std::max(grad_f_diss_dot_grad_f_k(e, n), 0.); // derivee en omega
              }
          }
      }
}


