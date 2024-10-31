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
// File:        Production_energie_cin_turb_PolyVEF.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/PolyVEF
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Production_energie_cin_turb_PolyVEF.h>

#include <Op_Diff_Turbulent_PolyVEF_Face.h>
#include <Echelle_temporelle_turbulente.h>
#include <Taux_dissipation_turbulent.h>
#include <Viscosite_turbulente_base.h>
#include <Domaine_PolyVEF.h>
#include <Navier_Stokes_std.h>
#include <Loi_paroi_base.h>
#include <Pb_Multiphase.h>
#include <Synonyme_info.h>

Implemente_instanciable(Production_energie_cin_turb_PolyVEF,"Production_energie_cin_turb_Elem_PolyVEF_P0", Source_Production_energie_cin_turb);
Add_synonym(Production_energie_cin_turb_PolyVEF, "Production_energie_cin_turb_Elem_PolyVEF_P0P1");
Add_synonym(Production_energie_cin_turb_PolyVEF, "Production_energie_cin_turb_Elem_PolyVEF_P0P1NC");

Sortie& Production_energie_cin_turb_PolyVEF::printOn(Sortie& os) const {return Source_Production_energie_cin_turb::printOn(os);}
Entree& Production_energie_cin_turb_PolyVEF::readOn(Entree& is) { return Source_Production_energie_cin_turb::readOn(is);}

void Production_energie_cin_turb_PolyVEF::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Domaine_PolyVEF&             domaine = ref_cast(Domaine_PolyVEF, equation().domaine_dis());
  const Probleme_base&                       pb = ref_cast(Probleme_base, equation().probleme());
  const Navier_Stokes_std&               eq_qdm = ref_cast(Navier_Stokes_std, pb.equation(0));
  const Viscosite_turbulente_base&    visc_turb = ref_cast(Viscosite_turbulente_base, ref_cast(Op_Diff_Turbulent_PolyVEF_Face, eq_qdm.operateur(0).l_op_base()).correlation());
  const DoubleVect& pe = equation().milieu().porosite_elem(), &pf = equation().milieu().porosite_face(), &ve = domaine.volumes();
  const IntTab& e_f = domaine.elem_faces();

  std::string Type_diss = ""; // omega or tau dissipation
  for (int i = 0 ; i < equation().probleme().nombre_d_equations() ; i++)
    {
      if      sub_type(Echelle_temporelle_turbulente, equation().probleme().equation(i)) Type_diss = "tau";
      else if sub_type(Taux_dissipation_turbulent, equation().probleme().equation(i)) Type_diss = "omega";
    }

  int D = dimension, Nph = pb.get_champ("vitesse").valeurs().dimension(1)/D, nb_elem = domaine.nb_elem() ;
  int N = equation().inconnue().valeurs().line_size();
  int e, f, n, i, ok;

  const DoubleTab& tab_rho = equation().probleme().get_champ("masse_volumique").passe(),
                   *palp = sub_type(Pb_Multiphase, pb) ? &equation().probleme().get_champ("alpha").passe() : nullptr,
                    &nu =  equation().probleme().get_champ("viscosite_cinematique").passe(),
                     // &k = equation().probleme().get_champ("k").valeurs(),
                     // &pk = semi_impl.count("k")  ? semi_impl.at("k") : equation().probleme().get_champ("k").passe(),
                     &tab_grad = pb.get_champ("gradient_vitesse").passe();
  // *diss = equation().probleme().has_champ(Type_diss) ? &equation().probleme().get_champ(Type_diss).valeurs() : nullptr,
  //  *pdiss = semi_impl.count(Type_diss) ? &semi_impl.at(Type_diss) : (equation().probleme().has_champ(Type_diss) ? &equation().probleme().get_champ(Type_diss).passe() : nullptr );

  const Loi_paroi_base& lp = ref_cast(Loi_paroi_base, pb.get_correlation("Loi_paroi"));
  // double min_ev_ratio = visc_turb.min_ev_ratio(), max_ev_ratio = visc_turb.max_ev_ratio();
  int c_nu = nu.dimension_tot(0) == 1;

  DoubleTrav nut(0, Nph);
  domaine.domaine().creer_tableau_elements(nut);
  visc_turb.eddy_viscosity(nut);

  if (Type_diss == "")
    {
      for( e = 0 ; e < nb_elem ; e++)
        {
          for (i = 0, ok = 1; ok && i < e_f.dimension(1) && (f = e_f(e, i)) >= 0; i++)
            ok &= std::abs(pf(f) - pe(e)) < 1e-8;
          if (!ok) break;
          for( n = 0; n < N ; n++)
            {
              double secmem_en = 0.;
              for (int d_U = 0; d_U < D; d_U++)
                for (int d_X = 0; d_X < D; d_X++)
                  secmem_en += ( tab_grad( e, Nph * ( D*d_U+d_X ) + n) + tab_grad( e,  Nph * ( D*d_X+d_U ) + n) ) * tab_grad( e,  Nph * ( D*d_U+d_X ) + n) ;
              secmem_en *= pe(e) * ve(e) * (palp ? (*palp)(e, n) : 1) * tab_rho(e, n) * nut(e, n) ;

              secmem(e, n) += std::max(secmem_en, 0.);
            }
        }
    }
  else
    {
      for( e = 0 ; e < nb_elem ; e++)
        {
          for (i = 0, ok = 1; ok && i < e_f.dimension(1) && (f = e_f(e, i)) >= 0; i++)
            ok &= std::abs(pf(f) - pe(e)) < 1e-8;
          if (!ok) break;
          for(n = 0; n<N ; n++)
            {
              double grad_grad = 0.;
              for (int d_U = 0; d_U < D; d_U++)
                for (int d_X = 0; d_X < D; d_X++)
                  grad_grad += ( tab_grad( e, Nph * ( D*d_U+d_X ) + n) + tab_grad( e,  Nph * ( D*d_X+d_U ) + n) ) * tab_grad( e,  Nph * ( D*d_U+d_X ) + n) ;
              double fac = std::max(grad_grad, 0.) * pe(e) * ve(e) ;

              double utau = 0., yloc = 0., nut_p = nut(e, n);
              int floc, cloc;
              for (cloc = 0 ; cloc < e_f.dimension(1) && (floc = e_f(e, cloc)) >= 0 ; cloc++)
                if (lp.get_utau(floc) > utau) utau = lp.get_utau(floc), yloc = domaine.dist_face_elem0(floc,e);

              fac += 0.5 * std::pow(utau,4)/(nut_p*nut_p) * pe(e) * ve(e) * std::tanh( std::pow(1./80. * utau*yloc/nu(!c_nu * e,n), 3) ) ;
              secmem(e, n) += fac * nut_p;

              // if (Type_diss == "tau")
              //   {
              //     secmem(e, n) += fac * nut_p;
              //     for (auto &&i_m : matrices)
              //       {
              //         Matrice_Morse& mat = *i_m.second;
              //         if (i_m.first == "k")         mat(N * e + n,  N * e + n) -= fac * (*diss)(e,n) ;
              //         if (i_m.first == "tau")       mat(N * e + n,  N * e + n) -= fac * k(e,n);
              //       }
              //   }
              // else if (Type_diss == "omega")
              //   {
              //     secmem(e, n) += fac * std::max(k(e, n) / std::max((*pdiss)(e, n), omega_min_), limiter_ * nu(!c_nu * e,n));
              //     if (k(e, n) / std::max((*pdiss)(e, n), omega_min_) >= limiter_ * nu(!c_nu * e,n))
              //     for (auto &&i_m : matrices)
              //       {
              //         Matrice_Morse& mat = *i_m.second;
              //         if (i_m.first == "k") mat(N * e + n,  N * e + n) -= fac / std::max((*pdiss)(e, n), omega_min_);
              //       }
              //   }
            }
        }
    }
}
