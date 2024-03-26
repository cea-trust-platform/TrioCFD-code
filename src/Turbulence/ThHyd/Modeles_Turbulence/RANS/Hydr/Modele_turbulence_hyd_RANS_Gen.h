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
// File:        Modele_turbulence_hyd_RANS_Gen.h
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/RANS/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modele_turbulence_hyd_RANS_Gen_included
#define Modele_turbulence_hyd_RANS_Gen_included

#include <Modifier_pour_fluide_dilatable.h>
#include <Schema_Temps_base.h>
#include <Champ_Inc_P0_base.h>
#include <Equation_base.h>
#include <Champ_Inc.h>
#include <type_traits>

enum class MODELE_TYPE { K_EPS , K_EPS_2_COUCHES, K_EPS_BAS_REYNOLDS , K_EPS_REALISABLE, K_OMEGA };

template <typename MODELE>
class Modele_turbulence_hyd_RANS_Gen
{
private:
  template <MODELE_TYPE M_TYPE>
  void print_evolution(const Champ_Inc& , const Schema_Temps_base& , double , int );

public:
  DoubleTab& complete_viscosity_field(const int , const Domaine_dis_base& , Champ_Inc& );

  template <MODELE_TYPE M_TYPE> std::enable_if_t<(M_TYPE !=  MODELE_TYPE::K_EPS_2_COUCHES && M_TYPE !=  MODELE_TYPE::K_EPS_BAS_REYNOLDS), void>
  calculate_limit_viscosity(Champ_Inc& , double );

  template <MODELE_TYPE M_TYPE> std::enable_if_t<M_TYPE ==  MODELE_TYPE::K_EPS_2_COUCHES, void>
  calculate_limit_viscosity(Champ_Inc& , double );

  template <MODELE_TYPE M_TYPE> std::enable_if_t<M_TYPE ==  MODELE_TYPE::K_EPS_BAS_REYNOLDS, void>
  calculate_limit_viscosity(Champ_Inc& , double );
};

template <typename MODELE>
DoubleTab& Modele_turbulence_hyd_RANS_Gen<MODELE>::complete_viscosity_field(const int n, const Domaine_dis_base& dom_dis, Champ_Inc& visc)
{
  visc->associer_domaine_dis_base(dom_dis);
  visc->nommer("diffusivite_turbulente");
  visc->fixer_nb_comp(1);
  visc->fixer_nb_valeurs_nodales(n);
  visc->fixer_unite("inconnue");
  visc->changer_temps(0.);
  return visc->valeurs();
}

template<typename MODELE> template<MODELE_TYPE M_TYPE>
std::enable_if_t<(M_TYPE !=  MODELE_TYPE::K_EPS_2_COUCHES && M_TYPE !=  MODELE_TYPE::K_EPS_BAS_REYNOLDS), void>
Modele_turbulence_hyd_RANS_Gen<MODELE>::calculate_limit_viscosity(Champ_Inc& ch_K_Eps_ou_Omega, double LeCmu)
{
  static constexpr bool IS_K_OMEGA = (M_TYPE == MODELE_TYPE::K_OMEGA);
  auto *z_class = static_cast<MODELE*>(this); // CRTP --> I love you :*
  const Milieu_base& mil = z_class->equation().probleme().milieu();

  // on divise par rho en QC pour revenir a K et Eps/Omega
  if (z_class->equation().probleme().is_dilatable())
    diviser_par_rho_si_dilatable(ch_K_Eps_ou_Omega.valeurs(), mil);

  print_evolution<M_TYPE>(ch_K_Eps_ou_Omega, z_class->equation().schema_temps(), LeCmu, 1);

  z_class->loi_paroi().calculer_hyd(ch_K_Eps_ou_Omega);
  z_class->controler();
  z_class->calculer_viscosite_turbulente(ch_K_Eps_ou_Omega.temps());
  z_class->limiter_viscosite_turbulente();

  // on remultiplie par rho
  if (z_class->equation().probleme().is_dilatable())
    {
      multiplier_par_rho_si_dilatable(ch_K_Eps_ou_Omega.valeurs(), mil);
      if (!IS_K_OMEGA)
        correction_nut_et_cisaillement_paroi_si_qc(*z_class);
    }

  z_class->viscosite_turbulente().valeurs().echange_espace_virtuel();

  print_evolution<M_TYPE>(ch_K_Eps_ou_Omega, z_class->equation().schema_temps(), LeCmu, 0);
}

template<typename MODELE> template<MODELE_TYPE M_TYPE>
std::enable_if_t<M_TYPE == MODELE_TYPE::K_EPS_2_COUCHES, void>
Modele_turbulence_hyd_RANS_Gen<MODELE>::calculate_limit_viscosity(Champ_Inc& ch_K_Eps_ou_Omega, double LeCmu)
{
  auto *z_class = static_cast<MODELE*>(this); // CRTP --> I love you :*
  const Milieu_base& mil = z_class->equation().probleme().milieu();

  z_class->controler();

  // on divise par rho en QC pour revenir a K et Eps/Omega
  if (z_class->equation().probleme().is_dilatable())
    diviser_par_rho_si_dilatable(ch_K_Eps_ou_Omega.valeurs(), mil);

  z_class->calculer_viscosite_turbulente(ch_K_Eps_ou_Omega.temps());
  z_class->loi_paroi().calculer_hyd(ch_K_Eps_ou_Omega);
  z_class->limiter_viscosite_turbulente();

  // on remultiplie par rho
  if (z_class->equation().probleme().is_dilatable())
    {
      multiplier_par_rho_si_dilatable(ch_K_Eps_ou_Omega.valeurs(), mil);
      correction_nut_et_cisaillement_paroi_si_qc(*z_class);
    }

  z_class->viscosite_turbulente().valeurs().echange_espace_virtuel();
}

template<typename MODELE> template<MODELE_TYPE M_TYPE>
std::enable_if_t<M_TYPE == MODELE_TYPE::K_EPS_BAS_REYNOLDS, void>
Modele_turbulence_hyd_RANS_Gen<MODELE>::calculate_limit_viscosity(Champ_Inc& ch_K_Eps_ou_Omega, double LeCmu)
{
  auto *z_class = static_cast<MODELE*>(this); // CRTP --> I love you :*
  z_class->controler();
  z_class->calculer_viscosite_turbulente(ch_K_Eps_ou_Omega.temps());
  z_class->limiter_viscosite_turbulente();
  z_class->viscosite_turbulente().valeurs().echange_espace_virtuel();
}

template <typename MODELE> template <MODELE_TYPE M_TYPE>
void Modele_turbulence_hyd_RANS_Gen<MODELE>::print_evolution(const Champ_Inc& le_champ_K_Eps_ou_Omega, const Schema_Temps_base& sch, const double LeCmu, const int avant)
{
  static constexpr bool IS_K_OMEGA = (M_TYPE == MODELE_TYPE::K_OMEGA), IS_2_COUCHES = (M_TYPE == MODELE_TYPE::K_EPS_2_COUCHES),
                        IS_BAS_REYNOLDS = (M_TYPE == MODELE_TYPE::K_EPS_BAS_REYNOLDS);

  if (IS_2_COUCHES || IS_BAS_REYNOLDS) return; /* Do nothing */

  if (sch.nb_pas_dt() == 0 || sch.limpr())
    {
      const DoubleTab& K_Eps_ou_Omega = le_champ_K_Eps_ou_Omega.valeurs();
      double k_min = DMAXFLOAT, eps_ou_omega_min = DMAXFLOAT, nut_min = DMAXFLOAT;
      double k_max = 0, eps_ou_omega_max = 0, nut_max = 0;
      int loc_k_min = -1, loc_eps_ou_omega_min = -1, loc_nut_min = -1;
      int loc_k_max = -1, loc_eps_ou_omega_max = -1, loc_nut_max = -1;
      int size = K_Eps_ou_Omega.dimension(0);

      if (size < 0)
        {
          if (!sub_type(Champ_Inc_P0_base, le_champ_K_Eps_ou_Omega.valeur()))
            {
              Cerr << "Unsupported field in Modele_turbulence_hyd_RANS_Gen::imprimer_evolution()" << finl;
              Process::exit(-1);
            }
          size = le_champ_K_Eps_ou_Omega.valeur().equation().domaine_dis().domaine().nb_elem();
        }

      for (int n = 0; n < size; n++)
        {
          const double k = K_Eps_ou_Omega(n, 0);
          const double epsOuomega = K_Eps_ou_Omega(n, 1);
          double nut = 0;

          const double num = IS_K_OMEGA ? k : LeCmu * k * k;

          if (epsOuomega > 0)
            nut = num / epsOuomega;

          if (k < k_min)
            {
              k_min = k;
              loc_k_min = n;
            }
          else if (k > k_max)
            {
              k_max = k;
              loc_k_max = n;
            }

          if (epsOuomega < eps_ou_omega_min)
            {
              eps_ou_omega_min = epsOuomega;
              loc_eps_ou_omega_min = n;
            }
          else if (epsOuomega > eps_ou_omega_max)
            {
              eps_ou_omega_max = epsOuomega;
              loc_eps_ou_omega_max = n;
            }

          if (nut < nut_min)
            {
              nut_min = nut;
              loc_nut_min = n;
            }
          else if (nut > nut_max)
            {
              nut_max = nut;
              loc_nut_max = n;
            }
        }

      ArrOfDouble values(3);

      // min values
      values[0] = k_min;
      values[1] = eps_ou_omega_min;
      values[2] = nut_min;
      mp_min_for_each_item(values);
      k_min = values[0];
      eps_ou_omega_min = values[1];
      nut_min = values[2];

      // max values
      values[0] = k_max;
      values[1] = eps_ou_omega_max;
      values[2] = nut_max;
      mp_max_for_each_item(values);
      k_max = values[0];
      eps_ou_omega_max = values[1];
      nut_max = values[2];

      // ecriture
      if (Process::je_suis_maitre())
        {
          Cout << finl << "K_Eps/Omega evolution (" << (avant ? "before" : "after") << " law of the wall applies) at time " << le_champ_K_Eps_ou_Omega.temps() << ":" << finl;
          Cout << "std::min(k)=" << k_min;
          if (Process::is_sequential()) Cout << " located at node " << loc_k_min;
          Cout << finl;
          Cout << "std::min(eps/omega)=" << eps_ou_omega_min;
          if (Process::is_sequential()) Cout << " located at node " << loc_eps_ou_omega_min;
          Cout << finl;
          Cout << "std::min(nut)=" << nut_min;
          if (Process::is_sequential()) Cout << " located at node " << loc_nut_min;
          Cout << finl;
          Cout << "std::max(k)=" << k_max;
          if (Process::is_sequential()) Cout << " located at node " << loc_k_max;
          Cout << finl;
          Cout << "std::max(eps/omega)=" << eps_ou_omega_max;
          if (Process::is_sequential()) Cout << " located at node " << loc_eps_ou_omega_max;
          Cout << finl;
          Cout << "std::max(nut)=" << nut_max;
          if (Process::is_sequential()) Cout << " located at node " << loc_nut_max;
          Cout << finl;
        }
    }
}

#endif /* Modele_turbulence_hyd_RANS_Gen_included */
