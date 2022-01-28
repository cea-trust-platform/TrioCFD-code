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
// File:        Diffusion_supplementaire_echelle_temp_turb_CoviMAC.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/CoviMAC
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Diffusion_supplementaire_echelle_temp_turb_CoviMAC.h>

#include <Zone_CoviMAC.h>
#include <Zone_Cl_CoviMAC.h>
#include <Champ_P0_CoviMAC.h>
#include <Pb_Multiphase.h>
#include <grad_Champ_Face_CoviMAC.h>
#include <Champ_Uniforme.h>
#include <Flux_interfacial_base.h>
#include <Milieu_composite.h>
#include <Operateur_Diff.h>
#include <Op_Diff_CoviMAC_base.h>
#include <Op_Diff_Turbulent_CoviMAC_Face.h>
#include <Op_Diff_Turbulent_CoviMAC_Elem.h>
#include <QDM_Multiphase.h>
#include <Viscosite_turbulente_k_tau.h>
#include <Energie_cinetique_turbulente.h>
#include <Echelle_temporelle_turbulente.h>
#include <Transport_turbulent_SGDH.h>
#include <Array_tools.h>
#include <Matrix_tools.h>
#include <Dirichlet.h>

#include <cmath>
#include <vector>

Implemente_instanciable(Diffusion_supplementaire_echelle_temp_turb_CoviMAC,"Diffusion_supplementaire_echelle_temp_turb_P0_CoviMAC", Source_base);

Sortie& Diffusion_supplementaire_echelle_temp_turb_CoviMAC::printOn(Sortie& os) const
{
  return os;
}

Entree& Diffusion_supplementaire_echelle_temp_turb_CoviMAC::readOn(Entree& is)
{
  return is;
}

void Diffusion_supplementaire_echelle_temp_turb_CoviMAC::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  const Zone_CoviMAC&       zone = ref_cast(Zone_CoviMAC, equation().zone_dis().valeur());
  const Champ_P0_CoviMAC&                   tau = ref_cast(Champ_P0_CoviMAC, equation().inconnue().valeur());

  if (!matrices.count("tau") || semi_impl.count("tau")) return;

  int N = equation().inconnue().valeurs().line_size(), ne = zone.nb_elem(), ne_tot = zone.nb_elem_tot() ;

  tau.init_grad(0);
  IntTab& f_d = tau.fgrad_d, f_e = tau.fgrad_e;             // Tables used in zone_CoviMAC::fgrad


  for (auto &&i_m : matrices) if (i_m.first == "tau")
      {
        Matrice_Morse mat;
        IntTrav stencil(0, 2);
        stencil.set_smart_resize(1);

        for(int e = 0 ; e < ne ; e++) for (int n=0 ; n<N ; n++) for (int i = zone.ved(e); i < zone.ved(e + 1); i++)
              {
                int f = zone.vej(i);
                for (int j = f_d(f); j < f_d(f+1) ; j++)
                  {
                    int ed = f_e(j);
                    if ( ed < ne_tot) //contrib d'un element ; contrib d'un bord : pas de derivee
                      {
                        stencil.append_line(N * e + n, N * ed + n);
                      }
                  }
              }
        tableau_trier_retirer_doublons(stencil);
        Matrix_tools::allocate_morse_matrix(ne_tot, ne_tot, stencil, mat);
        i_m.second->nb_colonnes() ? *i_m.second += mat : *i_m.second = mat;

      }
}

void Diffusion_supplementaire_echelle_temp_turb_CoviMAC::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Zone_CoviMAC&                      zone = ref_cast(Zone_CoviMAC, equation().zone_dis().valeur());
  const Zone_Cl_CoviMAC&                    zcl = ref_cast(Zone_Cl_CoviMAC, equation().zone_Cl_dis().valeur());
  const Echelle_temporelle_turbulente&       eq = ref_cast(Echelle_temporelle_turbulente, equation());
  const Navier_Stokes_std&               eq_qdm = ref_cast(Navier_Stokes_std, equation().probleme().equation(0));
  const Champ_P0_CoviMAC&                   tau = ref_cast(Champ_P0_CoviMAC, equation().inconnue().valeur());
  const DoubleTab&                      tab_tau = semi_impl.count("tau") ? tau.passe() : tau.valeurs();
  const DoubleTab&                  tab_k_passe = equation().probleme().get_champ("k").passe();
  const IntTab&                             fcl = tau.fcl();
  const Conds_lim&                          cls = zcl.les_conditions_limites();
  const Op_Diff_Turbulent_CoviMAC_Elem& Op_diff_loc = ref_cast(Op_Diff_Turbulent_CoviMAC_Elem, eq.operateur(0).l_op_base());
  const DoubleTab&                       mu_tot = Op_diff_loc.nu();
  const Op_Diff_Turbulent_CoviMAC_Face& op_diff = ref_cast(Op_Diff_Turbulent_CoviMAC_Face, eq_qdm.operateur(0).l_op_base());
  const Viscosite_turbulente_k_tau&   visc_turb = ref_cast(Viscosite_turbulente_k_tau, op_diff.corr.valeur());
  const DoubleTab&                      nu_visc	= equation().probleme().get_champ("viscosite_cinematique").passe();

  int N = tab_tau.dimension(1), nf = zone.nb_faces(), ne = zone.nb_elem(), ne_tot = zone.nb_elem_tot(), nf_tot = zone.nb_faces_tot(), D = dimension ;

  DoubleTrav grad_sqrt_tau(eq_qdm.vitesse()->valeurs().dimension_tot(0), N);
  tau.init_grad(0);
  IntTab& f_d = tau.fgrad_d, f_e = tau.fgrad_e;             // Tables used in zone_CoviMAC::fgrad
  DoubleTab f_w = tau.fgrad_w;

  // Calculation of grad of root of tau at faces
  Matrice_Morse *Mtau = (matrices.count("tau") && !semi_impl.count("tau")) ? matrices.at("tau") : nullptr;

  for (int f = 0; f < nf; f++) for (int n=0 ; n<N ; n++)
      {
        grad_sqrt_tau(f, n) = 0;
        for (int j = f_d(f); j < f_d(f+1) ; j++)
          {
            int e = f_e(j);
            int f_bord;
            if ( e < ne_tot) //contrib d'un element
              {
                double val_e = sqrt(std::max(0., tab_tau(e, n))) ; // Si implicite ou pas
                grad_sqrt_tau(f, n) += f_w(j) * val_e;
              }
            else if (fcl(f_bord = e - ne_tot, 0) == 3) //contrib d'un bord : seul Dirichlet contribue
              {
                double val_f_bord = sqrt(ref_cast(Dirichlet, cls[fcl(f_bord, 1)].valeur()).val_imp(fcl(f_bord, 2), n));
                grad_sqrt_tau(f, n) += f_w(j) * val_f_bord;
              }
          }
      }

  // Calculation of grad of root of tau at elements ; same if implicit or explicit
  zone.init_ve();
  for (int e = 0; e < ne_tot; e++) for (int n=0 ; n<N ; n++) for (int d = 0; d < D; d++)
        {
          grad_sqrt_tau(nf_tot + D*e+d, n) = 0;
          for (int j = zone.ved(e); j < zone.ved(e + 1); j++)
            {
              int f = zone.vej(j);
              grad_sqrt_tau(nf_tot + D*e+d, n) += zone.vec(j, d) * grad_sqrt_tau(f, n);
            }
        }


  // Creation of a limiter so that (|additionnal diffusion| < |5*diffusion|) ; on regarde seulement le second membre

  DoubleTrav sec_m(tab_tau); //residus
  matrices_t mat_m; //derivees
  eq.operateur(0).l_op_base().ajouter_blocs(mat_m, sec_m, semi_impl);

  for(int e = 0 ; e < ne ; e++) for (int n=0 ; n<N ; n++)
      {
        // Second membre
        double secmem_en = 0 ;
        for (int d = 0; d<D; d++) secmem_en += grad_sqrt_tau(nf_tot + D*e+d, n)*grad_sqrt_tau(nf_tot + D*e+d, n) ;
        secmem_en *= -8  * mu_tot(e, n);

        // limiteur
        if (limiter_*std::abs(sec_m(e,n)) < std::abs(secmem_en)) secmem_en *= limiter_*std::abs(sec_m(e,n))/std::abs(secmem_en);
        // Jacobienne ; si on est dans le limiteur, pas de derivee
        else if (!(Mtau == nullptr))
          {
            for (int i = zone.ved(e); i < zone.ved(e + 1); i++)
              {
                int f = zone.vej(i);
                for (int j = f_d(f); j < f_d(f+1) ; j++)
                  {
                    int ed = f_e(j);
                    if ( ed < ne_tot) for (int d = 0 ; d<D ; d++) //contrib d'un element ; contrib d'un bord : pas de derivee
                        {
                          double inv_sqrt_tau = sqrt(tab_k_passe(ed, n) / std::max(tab_k_passe(ed, n) * tab_tau(ed, n), visc_turb.limiteur() * nu_visc(ed, n)));
                          (*Mtau)(N * e + n, N * ed + n) += 8  * mu_tot(e, n) * sqrt(2) * grad_sqrt_tau(nf_tot + D*ed+d, n) * zone.vec(i, d) * f_w(j) * inv_sqrt_tau;
                        }
                  }
              }
          }

        secmem(e, n) += secmem_en;

      }
}


//        double inv_tau = tab_k_passe(e, n) / std::max(tab_k_passe(e, n) * tab_tau_passe(e, n), visc_turb.limiteur() * nu_visc(e, n));
//       double inv_tau = tab_k_passe(e, n) / std::max(tab_k_passe(e, n) * tab_tau(e, n), visc_turb.limiteur() * nu_visc(e, n));
//       double secmem_en = -2  * inv_tau *nu_tot(e, n) * sq_grad_sqrt_tau(e, n) ;

/*        if (!(Mtau==nullptr))           	 // derivee en tau
          {
            if (tab_k_passe(e, n) * tab_tau(e, n) > visc_turb.limiteur() * nu_visc(e, n))
              (*Mtau)(N * e + n, N * e + n)   += (-2) * nu_tot(e, n) * sq_grad_sqrt_tau(e, n) / (tab_tau(e, n)*tab_tau(e, n));
            else
              (*Mtau)(N * e + n, N * e + n) += 0;
          }*/
