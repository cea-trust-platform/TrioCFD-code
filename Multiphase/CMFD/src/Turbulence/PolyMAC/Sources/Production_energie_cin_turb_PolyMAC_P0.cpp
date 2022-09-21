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
// File:        Production_energie_cin_turb_PolyMAC_P0.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/PolyMAC_P0
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Production_energie_cin_turb_PolyMAC_P0.h>

#include <Zone_PolyMAC_P0.h>
#include <Champ_Elem_PolyMAC_P0.h>
#include <Matrix_tools.h>
#include <Probleme_base.h>
#include <grad_Champ_Face_PolyMAC_P0.h>
#include <Champ_Uniforme.h>
#include <Flux_interfacial_base.h>
#include <Milieu_composite.h>
#include <Operateur_Diff.h>
#include <Op_Diff_Turbulent_PolyMAC_P0_Face.h>
#include <Navier_Stokes_std.h>
#include <Viscosite_turbulente_base.h>
#include <TRUSTTab_parts.h>
#include <Loi_paroi_adaptative.h>
#include <Pb_Multiphase.h>

Implemente_instanciable(Production_energie_cin_turb_PolyMAC_P0,"Production_energie_cin_turb_Elem_PolyMAC_P0", Source_base);
// XD Production_energie_cin_turb source_base Production_energie_cin_turb 0 Production source term for the TKE equation


Sortie& Production_energie_cin_turb_PolyMAC_P0::printOn(Sortie& os) const
{
  return os;
}

Entree& Production_energie_cin_turb_PolyMAC_P0::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("limiter_production", &limiter_prod_);
  param.lire_avec_accolades_depuis(is);

  return is;
}

void Production_energie_cin_turb_PolyMAC_P0::completer()
{
  const Zone_PolyMAC_P0&                      zone = ref_cast(Zone_PolyMAC_P0, equation().zone_dis().valeur());
  fac_.resize(zone.nb_elem_tot(), equation().inconnue().valeur().valeurs().line_size(), 2); // 3rd column : derivative
  fac_ = 1. ;

  for (int e = 0 ; e < fac_.dimension_tot(0) ; e++)
    for (int n = 0 ; n < fac_.dimension_tot(1) ; n++) fac_(e, n, 1) = 0. ; // 0 derivative by default

  if (ref_cast(Pb_Multiphase, equation().probleme()).has_correlation("Loi_paroi"))
    correlation_loi_paroi_ = ref_cast(Pb_Multiphase, equation().probleme()).get_correlation("Loi_paroi");
}

void Production_energie_cin_turb_PolyMAC_P0::mettre_a_jour(double tps)
{
  if (tps != mon_temps_) {calculer_fac() ; mon_temps_ = tps;}
}
void Production_energie_cin_turb_PolyMAC_P0::calculer_fac()
{
  if (correlation_loi_paroi_.non_nul())
    {
      const Zone_PolyMAC_P0&                      zone = ref_cast(Zone_PolyMAC_P0, equation().zone_dis().valeur());
      const IntTab& f_e = zone.face_voisins();

      Loi_paroi_adaptative& corr_loi_paroi = ref_cast(Loi_paroi_adaptative, correlation_loi_paroi_.valeur().valeur());
      DoubleTab& u_tau = corr_loi_paroi.get_tab("u_tau");
      const DoubleTab& tab_k = equation().probleme().get_champ("k").passe();

      int nf = zone.nb_faces(), N = tab_k.line_size() ;

      fac_ = 1.;
      for (int e = 0 ; e < fac_.dimension_tot(0) ; e++)
        for (int n = 0 ; n < fac_.dimension_tot(1) ; n++) fac_(e, n, 1) = 0. ; // 0 derivative by default

      for(int f = 0 ; f < nf ; f++)
        for (int n = 0 ; n < N ; n++)
          if (corr_loi_paroi.a_calculer(f))
            if (u_tau(f, n) > 1.e-6)
              {
                int e = f_e(f, 0) ;

                double kp = tab_k(e, n) / (u_tau(f, n)*u_tau(f, n));
                if ( (kp-3.33) > 1. ) fac_(e, n, 0) = 1. - limiter_prod_ ;
                else if ( (kp-3.33) < 0. ) fac_(e, n, 0) = 1. ;
                else fac_(e, n, 0) = 1. - limiter_prod_*(kp-3.33), fac_(e, n, 1) = - limiter_prod_  / (u_tau(f, n)*u_tau(f, n)) ;
              }
    }
}

void Production_energie_cin_turb_PolyMAC_P0::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
// empty : no derivative for the turbulent kinetic energy production to add in the blocks
}

void Production_energie_cin_turb_PolyMAC_P0::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Zone_PolyMAC_P0&                   zone = ref_cast(Zone_PolyMAC_P0, equation().zone_dis().valeur());
  const Probleme_base&                       pb = ref_cast(Probleme_base, equation().probleme());
  const Navier_Stokes_std&               eq_qdm = ref_cast(Navier_Stokes_std, pb.equation(0));
  const DoubleTab&                     tab_grad = pb.get_champ("gradient_vitesse").passe();
  const Op_Diff_Turbulent_PolyMAC_P0_Face& Op_diff = ref_cast(Op_Diff_Turbulent_PolyMAC_P0_Face, eq_qdm.operateur(0).l_op_base());
  const Viscosite_turbulente_base&    visc_turb = ref_cast(Viscosite_turbulente_base, Op_diff.correlation().valeur());
  const DoubleTab&                      tab_rho = equation().probleme().get_champ("masse_volumique").passe();
  const DoubleTab&                      tab_alp = equation().probleme().get_champ("alpha").passe();
  const DoubleVect& pe = equation().milieu().porosite_elem(), &ve = zone.volumes();

  int Nph = pb.get_champ("vitesse").valeurs().dimension(1), nb_elem = zone.nb_elem(), D = dimension, nf_tot = zone.nb_faces_tot() ;
  int N = equation().inconnue()->valeurs().line_size();

//  DoubleTrav Rij(0, Nph, D, D);
// MD_Vector_tools::creer_tableau_distribue(eq_qdm.pression()->valeurs().get_md_vector(), Rij); //Necessary to compare size in reynolds_stress()
//  visc_turb.reynolds_stress(Rij);

  DoubleTrav nut(0, Nph);
  MD_Vector_tools::creer_tableau_distribue(eq_qdm.pression()->valeurs().get_md_vector(), nut); //Necessary to compare size in reynolds_stress()
  visc_turb.eddy_viscosity(nut);

  Matrice_Morse *mat = matrices.count("k") ? matrices.at("k") : nullptr;

  for(int e = 0 ; e < nb_elem ; e++)
    for(int n = 0; n<N ; n++)
      {
        double secmem_en = 0.;
        for (int d_U = 0; d_U < D; d_U++)
          for (int d_X = 0; d_X < D; d_X++)
            secmem_en += ( tab_grad(nf_tot + d_U + e * D , D * n + d_X) + tab_grad(nf_tot + d_X + e * D , D * n + d_U) ) * tab_grad(nf_tot + d_X + e * D , D * n + d_U) ;
//            secmem_en += Rij(e, n, d_X, d_U) * tab_grad(nf_tot + d_X + e * D , D * n + d_U) ;
        secmem_en *= pe(e) * ve(e) * tab_alp(e, n) * tab_rho(e, n) * nut(e, n) ;
//        secmem_en *= (-1) * pe(e) * ve(e) * tab_alp(e, n) * tab_rho(e, n) ;

//        secmem(e, n) += fac_(e, n, 0) * std::max(secmem_en, 0.);
        secmem(e, n) += std::max(secmem_en, 0.);

        if (mat) (*mat)(N * e + n, N * e + n) -= 0.;//fac_(e, n, 1) * std::max(secmem_en, 0.);
      }
}
