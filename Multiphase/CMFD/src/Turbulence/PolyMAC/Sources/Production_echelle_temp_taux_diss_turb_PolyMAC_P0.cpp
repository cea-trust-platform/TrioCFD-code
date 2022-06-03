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
// File:        Production_echelle_temp_taux_diss_turb_PolyMAC_P0.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/PolyMAC_P0
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Production_echelle_temp_taux_diss_turb_PolyMAC_P0.h>

#include <Zone_PolyMAC_P0.h>
#include <Champ_Elem_PolyMAC_P0.h>
#include <Probleme_base.h>
#include <grad_Champ_Face_PolyMAC_P0.h>
#include <Champ_Uniforme.h>
#include <Flux_interfacial_base.h>
#include <Milieu_composite.h>
#include <Operateur_Diff.h>
#include <Op_Diff_Turbulent_PolyMAC_P0_Face.h>
#include <Navier_Stokes_std.h>
#include <Viscosite_turbulente_k_tau.h>
#include <Energie_cinetique_turbulente.h>
#include <Array_tools.h>
#include <Matrix_tools.h>
#include <Echelle_temporelle_turbulente.h>
#include <Taux_dissipation_turbulent.h>

Implemente_instanciable(Production_echelle_temp_taux_diss_turb_PolyMAC_P0,"Production_echelle_temp_taux_diss_turb_Elem_PolyMAC_P0", Source_base);

Sortie& Production_echelle_temp_taux_diss_turb_PolyMAC_P0::printOn(Sortie& os) const
{
  return os;
}

Entree& Production_echelle_temp_taux_diss_turb_PolyMAC_P0::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("alpha_omega", &alpha_omega_, Param::REQUIRED);
  param.lire_avec_accolades_depuis(is);
  Cout << "alpha_omega = " << alpha_omega_ << finl ;

  return is;
}

void Production_echelle_temp_taux_diss_turb_PolyMAC_P0::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  const Zone_PolyMAC_P0&       zone = ref_cast(Zone_PolyMAC_P0, equation().zone_dis().valeur());
  int ne = zone.nb_elem(), ne_tot = zone.nb_elem_tot(), N = equation().inconnue().valeurs().line_size();

  for (auto &&i_m : matrices) if (i_m.first == "k")
      {
        Matrice_Morse mat;
        IntTrav stencil(0, 2);
        stencil.set_smart_resize(1);
        for (int e = 0; e < ne; e++) for(int n = 0; n<N ; n++) stencil.append_line(N * e + n, N * e + n);
        tableau_trier_retirer_doublons(stencil);
        Matrix_tools::allocate_morse_matrix(ne_tot, ne_tot, stencil, mat);
        i_m.second->nb_colonnes() ? *i_m.second += mat : *i_m.second = mat;
      }
}

void Production_echelle_temp_taux_diss_turb_PolyMAC_P0::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Zone_PolyMAC_P0&                      zone = ref_cast(Zone_PolyMAC_P0, equation().zone_dis().valeur());
  const Probleme_base&                       pb = ref_cast(Probleme_base, equation().probleme());
  const Navier_Stokes_std&                  eq_qdm = ref_cast(Navier_Stokes_std, pb.equation(0));
  const grad_Champ_Face_PolyMAC_P0&           grad = ref_cast(grad_Champ_Face_PolyMAC_P0, eq_qdm.get_champ("gradient_vitesse"));
  const DoubleTab&                     tab_grad = grad.passe();
  const Op_Diff_Turbulent_PolyMAC_P0_Face& Op_diff = ref_cast(Op_Diff_Turbulent_PolyMAC_P0_Face, eq_qdm.operateur(0).l_op_base());
  const Viscosite_turbulente_base&    visc_turb = ref_cast(Viscosite_turbulente_base, Op_diff.correlation().valeur());
  const DoubleTab&                     tab_diss = ref_cast(Champ_Elem_PolyMAC_P0, equation().inconnue().valeur()).valeurs(); // tau ou omega selon l'equation
  const DoubleTab&                           nu = pb.get_champ("viscosite_cinematique").passe();
  const DoubleTab&                        tab_k = ref_cast(Champ_Elem_PolyMAC_P0, pb.get_champ("k")).valeurs();
  const DoubleTab&                      tab_rho = equation().probleme().get_champ("masse_volumique").passe();
  const DoubleTab&                      tab_alp = equation().probleme().get_champ("alpha").passe();
  const DoubleVect& pe = zone.porosite_elem(), &ve = zone.volumes();

  int Nph = pb.get_champ("vitesse").valeurs().line_size(), nf_tot = zone.nb_faces_tot(), ne = zone.nb_elem(), D = dimension ;
  int N = tab_diss.line_size();

  DoubleTrav Rij(0, Nph, D, D);
  MD_Vector_tools::creer_tableau_distribue(eq_qdm.pression()->valeurs().get_md_vector(), Rij); //Necessary to compare size in reynolds_stress()
  visc_turb.reynolds_stress(Rij);
  assert((ne == Rij.dimension(0)));

  std::string Type_diss = ""; // omega or tau dissipation
  if sub_type(Echelle_temporelle_turbulente, equation()) Type_diss = "tau";
  else if sub_type(Taux_dissipation_turbulent, equation()) Type_diss = "omega";
  if (Type_diss == "") abort();

  DoubleTrav prod_scal(Rij.dimension_tot(0), Nph);
  for (int e = 0; e < ne; e++) for(int n = 0; n<N ; n++)
      {
        prod_scal(e, n) = 0;
        for (int d_U = 0; d_U < D; d_U++) for (int d_X = 0; d_X < D; d_X++)
            prod_scal(e, n) += Rij(e, n, d_U, d_X) * tab_grad(nf_tot + d_X + e * D , D * n + d_U) ;
      }

  // Derivees
  for (auto &&i_m : matrices)
    {
      if (i_m.first == "tau")
        {
          Matrice_Morse& mat = *i_m.second;
          for (int e = 0; e < ne; e++) for(int n = 0; n<N ; n++)
              {
                double deriv = std::min(prod_scal(e, n),0.) ; // So the production is always negative for tau and positive for omega
                if (tab_k(e, n) * tab_diss(e, n) > visc_turb.limiteur() * nu(e, n))
                  deriv *=    pe(e) * ve(e) * (-1) * alpha_omega_* tab_alp(e, n) * tab_rho(e, n)/tab_k(e, n) ;
                else deriv *= pe(e) * ve(e) * (-2) * alpha_omega_* tab_alp(e, n) * tab_rho(e, n) * tab_diss(e, n)/(visc_turb.limiteur() * nu(e, n));
                mat(N * e + n, N * e + n) += deriv;
              }
        }

      else if (i_m.first == "omega")
        {
          Matrice_Morse& mat = *i_m.second;
          for (int e = 0; e < ne; e++) for(int n = 0; n<N ; n++)
              {
                double deriv = std::min(prod_scal(e, n),0.) ; // So the production is always negative for tau and positive for omega
                if (tab_diss(e, n) <= 0) deriv *= 0 ;
                else if (tab_k(e, n)/tab_diss(e, n)<visc_turb.limiteur() * nu(e, n)) deriv *= 0 ;
                else deriv *= pe(e) * ve(e) * alpha_omega_* tab_alp(e, n) * tab_rho(e, n) / tab_k(e, n) ;
                mat(N * e + n, N * e + n) += deriv;
              }
        }

      else if (i_m.first == "k")
        {
          Matrice_Morse& mat = *i_m.second;
          if (Type_diss == "tau")
            {
              for (int e = 0; e < ne; e++) for(int n = 0; n<N ; n++)
                  {
                    double deriv = std::min(prod_scal(e, n),0.) ; // So the production is always negative for tau and positive for omega
                    if (tab_k(e, n) * tab_diss(e, n) > visc_turb.limiteur() * nu(e, n))
                      deriv *= pe(e) * ve(e) * alpha_omega_* tab_alp(e, n) * tab_rho(e, n) * tab_diss(e, n)/(tab_k(e, n)*tab_k(e, n)) ;
                    else deriv *= 0;
                    mat(N * e + n, N * e + n) += deriv;
                  }
            }
          else if (Type_diss == "omega")
            for (int e = 0; e < ne; e++) for(int n = 0; n<N ; n++)
                {
                  double deriv = std::min(prod_scal(e, n),0.) ; // So the production is always negative for tau and positive for omega
                  if (tab_diss(e, n) <= 0) deriv *= 0 ;
                  else if (tab_k(e, n)/tab_diss(e, n)<visc_turb.limiteur() * nu(e, n)) deriv *= 0 ;
                  else deriv *= pe(e) * ve(e) * (-1) * alpha_omega_* tab_alp(e, n) * tab_rho(e, n) * tab_diss(e, n)/ (tab_k(e, n)*tab_k(e, n)) ;
                  mat(N * e + n, N * e + n) += deriv;
                }


        }
    }

  // Second membre
  for(int e = 0 ; e < ne ; e++) for(int n = 0; n<N ; n++)
      {
        double secmem_en = std::min(prod_scal(e, n),0.) ; // So the production is always negative for tau and positive for omega
        if (Type_diss == "tau")
          secmem_en *= pe(e) * ve(e) *      alpha_omega_* tab_alp(e, n) * tab_rho(e, n)*tab_diss(e, n)*tab_diss(e, n)/std::max(tab_k(e, n) * tab_diss(e, n), visc_turb.limiteur() * nu(e, n)) ;
        else if (Type_diss == "omega")
          secmem_en *= pe(e) * ve(e) * (-1)*alpha_omega_* tab_alp(e, n) * tab_rho(e, n)* ((tab_diss(e, n) <= 0) ? 0 : 1/(std::max(tab_k(e, n)/tab_diss(e, n), visc_turb.limiteur() * nu(e, n)))) ;
        secmem(e, n) += secmem_en;
      }
}
