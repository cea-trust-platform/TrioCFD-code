/****************************************************************************
* Copyright (c) 2022, CEA
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
// File:        Frottement_interfacial_CoviMAC.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/CoviMAC
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Portance_interfaciale_CoviMAC.h>
#include <Portance_interfaciale_base.h>
#include <Zone_CoviMAC.h>
#include <Champ_Face_CoviMAC.h>
#include <Op_Grad_CoviMAC_Face.h>
#include <Zone_Cl_CoviMAC.h>
#include <Array_tools.h>
#include <Matrix_tools.h>
#include <Pb_Multiphase.h>
#include <Champ_Uniforme.h>
#include <Milieu_composite.h>
#include <Interface_base.h>

Implemente_instanciable(Portance_interfaciale_CoviMAC,"Portance_interfaciale_Face_CoviMAC", Source_base);

Sortie& Portance_interfaciale_CoviMAC::printOn(Sortie& os) const
{
  return os;
}

Entree& Portance_interfaciale_CoviMAC::readOn(Entree& is)
{
  const Pb_Multiphase& pbm = ref_cast(Pb_Multiphase, equation().probleme());
  if (pbm.has_correlation("Portance_interfaciale")) correlation_ = pbm.get_correlation("Portance_interfaciale"); //correlation fournie par le bloc correlation
  else correlation_.typer_lire(pbm, "Portance_interfaciale", is); //sinon -> on la lit
  return is;
}

void Portance_interfaciale_CoviMAC::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const // 100% explicit
{
}

void Portance_interfaciale_CoviMAC::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Champ_Face_CoviMAC& ch = ref_cast(Champ_Face_CoviMAC, equation().inconnue().valeur());
  const Zone_CoviMAC& zone = ref_cast(Zone_CoviMAC, equation().zone_dis().valeur());
  const IntTab& f_e = zone.face_voisins(), &fcl = ch.fcl(), &e_f = zone.elem_faces();
  const DoubleTab& n_f = zone.face_normales();
  const DoubleVect& pe = zone.porosite_elem(), &pf = zone.porosite_face(), &ve = zone.volumes(), &vf = zone.volumes_entrelaces(), &fs = zone.face_surfaces(), &dh_e = zone.diametre_hydraulique_elem();
  const DoubleTab& pvit = ch.passe(), &mu_f = ref_cast(Op_Grad_CoviMAC_Face, ref_cast(Navier_Stokes_std, equation()).operateur_gradient().valeur()).mu_f(),
                   &alpha = ref_cast(Pb_Multiphase, equation().probleme()).eq_masse.inconnue().passe(),
                    &press = ref_cast(QDM_Multiphase, equation()).pression().passe(),
                     &temp  = ref_cast(Pb_Multiphase, equation().probleme()).eq_energie.inconnue().passe(),
                      &rho   = equation().milieu().masse_volumique().passe(),
                       &mu    = ref_cast(Fluide_base, equation().milieu()).viscosite_dynamique().passe(),
                        &vort  = equation().probleme().get_champ("vorticite").passe();

  const Milieu_composite& milc = ref_cast(Milieu_composite, equation().milieu());

  int e, f, b, c, i, k, l, n, N = ch.valeurs().line_size(), Np = press.line_size(), D = dimension,
                              cR = (rho.dimension_tot(0) == 1), cM = (mu.dimension_tot(0) == 1);
  DoubleTrav a_l(N), p_l(N), T_l(N), rho_l(N), mu_l(N), sigma_l(N,N), dv(N, N), ddv(N, N, 4), ddv_c(4), coeff(N, N); //arguments pour coeff
  double dv_min = 0.1, fac_e, fac_f;
  const Portance_interfaciale_base& correlation_pi = ref_cast(Portance_interfaciale_base, correlation_.valeur());

  /* elements */
  for (e = 0; e < zone.nb_elem_tot(); e++)   for (b = 0; b < e_f.dimension(1) && (f = e_f(e, b)) >= 0; b++)
      if ( (e < zone.nb_elem()) || (f<zone.nb_faces()))
        {
          /* arguments de coeff */
          for (n = 0; n < N; n++) a_l(n)   = alpha(e, n);
          for (n = 0; n < N; n++) p_l(n)   = press(e, n * (Np > 1));
          for (n = 0; n < N; n++) T_l(n)   =  temp(e, n);
          for (n = 0; n < N; n++) rho_l(n) =   rho(!cR * e, n);
          for (n = 0; n < N; n++) mu_l(n)  =    mu(!cM * e, n);
          for (n = 0; n < N; n++)
            {
              for (k = 0; k < N; k++) if(milc.has_interface(n, k))
                  {
                    Interface_base& sat = milc.get_interface(n, k);
                    sigma_l(n,k) = sat.sigma_(temp(e,n), press(e,n * (Np > 1)));
                  }
            }

          for (k = 0; k < N; k++) for (l = 0; l < N; l++) dv(k, l) = std::max(ch.v_norm(pvit, pvit, e, -1, k, l, NULL, &ddv(k, l, 0)), dv_min);
          correlation_pi.coefficient(a_l, p_l, T_l, rho_l, mu_l, sigma_l, dh_e(e), dv, e, coeff);

          fac_e = pe(e) * ve(e);
          k=0;

          if (D==2)
            {
              if (e < zone.nb_elem())
                for (l = 1; l < N; l++)
                  {
                    i = zone.nb_faces_tot() + D * e;
                    secmem(i, k) += fac_e * coeff(k, l) * (pvit(i+1, l) -pvit(i+1, k)) * vort(e, 0) ;
                    secmem(i, l) -= fac_e * coeff(k, l) * (pvit(i+1, l) -pvit(i+1, k)) * vort(e, 0) ;
                    secmem(i+1,k)-= fac_e * coeff(k, l) * (pvit(i  , l) -pvit(i  , k)) * vort(e, 0) ;
                    secmem(i+1,l)+= fac_e * coeff(k, l) * (pvit(i  , l) -pvit(i  , k)) * vort(e, 0) ;
                  } // 100% explicit

              if (f<zone.nb_faces()) if (fcl(f, 0) < 2)
                  {
                    c = (e == f_e(f,0)) ? 0 : 1 ;
                    fac_f = pf(f) * vf(f);
                    for (l = 1; l < N; l++)
                      {
                        secmem(f, k) += fac_f * mu_f(f, k, c) * n_f(f, 0)/fs(f) * coeff(k, l) * (pvit(i+1, l) -pvit(i+1, k)) * vort(e, 0) ;
                        secmem(f, l) -= fac_f * mu_f(f, l, c) * n_f(f, 0)/fs(f) * coeff(k, l) * (pvit(i+1, l) -pvit(i+1, k)) * vort(e, 0) ;
                        secmem(f, k) -= fac_f * mu_f(f, k, c) * n_f(f, 1)/fs(f) * coeff(k, l) * (pvit(i  , l) -pvit(i  , k)) * vort(e, 0) ;
                        secmem(f, l) += fac_f * mu_f(f, l, c) * n_f(f, 1)/fs(f) * coeff(k, l) * (pvit(i  , l) -pvit(i  , k)) * vort(e, 0) ;
                      } // 100% explicit

                  }
            }

          if (D==3)
            {
              if (e < zone.nb_elem())
                for (l = 1; l < N; l++)
                  {
                    i = zone.nb_faces_tot() + D * e;
                    secmem(i, k) += fac_e * coeff(k, l) * ((pvit(i+1, l) -pvit(i+1, k)) * vort(D*e+ 2, k) - (pvit(i+2, l) -pvit(i+2, k)) * vort(D*e+ 1, k)) ;
                    secmem(i, l) -= fac_e * coeff(k, l) * ((pvit(i+1, l) -pvit(i+1, k)) * vort(D*e+ 2, k) - (pvit(i+2, l) -pvit(i+2, k)) * vort(D*e+ 1, k)) ;
                    secmem(i+1,k)+= fac_e * coeff(k, l) * ((pvit(i+2, l) -pvit(i+2, k)) * vort(D*e+ 0, k) - (pvit(i+0, l) -pvit(i+0, k)) * vort(D*e+ 2, k)) ;
                    secmem(i+1,l)-= fac_e * coeff(k, l) * ((pvit(i+2, l) -pvit(i+2, k)) * vort(D*e+ 0, k) - (pvit(i+0, l) -pvit(i+0, k)) * vort(D*e+ 2, k)) ;
                    secmem(i+2,k)+= fac_e * coeff(k, l) * ((pvit(i+0, l) -pvit(i+0, k)) * vort(D*e+ 1, k) - (pvit(i+1, l) -pvit(i+1, k)) * vort(D*e+ 0, k)) ;
                    secmem(i+2,l)-= fac_e * coeff(k, l) * ((pvit(i+0, l) -pvit(i+0, k)) * vort(D*e+ 1, k) - (pvit(i+1, l) -pvit(i+1, k)) * vort(D*e+ 0, k)) ;
                  } // 100% explicit

              if (f<zone.nb_faces()) if (fcl(f, 0) < 2)
                  {
                    c = (e == f_e(f,0)) ? 0 : 1 ;
                    fac_f = pf(f) * vf(f);
                    for (l = 1; l < N; l++)
                      {
                        secmem(f, k) += fac_f * mu_f(f, k, c) * n_f(f, 0)/fs(f) * coeff(k, l) * ((pvit(i+1, l) -pvit(i+1, k)) * vort(D*e+ 2, k) - (pvit(i+2, l) -pvit(i+2, k)) * vort(D*e+ 1, k)) ;
                        secmem(f, l) -= fac_f * mu_f(f, l, c) * n_f(f, 0)/fs(f) * coeff(k, l) * ((pvit(i+1, l) -pvit(i+1, k)) * vort(D*e+ 2, k) - (pvit(i+2, l) -pvit(i+2, k)) * vort(D*e+ 1, k)) ;
                        secmem(f, k) += fac_f * mu_f(f, k, c) * n_f(f, 0)/fs(f) * coeff(k, l) * ((pvit(i+2, l) -pvit(i+2, k)) * vort(D*e+ 0, k) - (pvit(i+0, l) -pvit(i+0, k)) * vort(D*e+ 2, k)) ;
                        secmem(f, l) -= fac_f * mu_f(f, l, c) * n_f(f, 0)/fs(f) * coeff(k, l) * ((pvit(i+2, l) -pvit(i+2, k)) * vort(D*e+ 0, k) - (pvit(i+0, l) -pvit(i+0, k)) * vort(D*e+ 2, k)) ;
                        secmem(f, k) += fac_f * mu_f(f, k, c) * n_f(f, 0)/fs(f) * coeff(k, l) * ((pvit(i+0, l) -pvit(i+0, k)) * vort(D*e+ 1, k) - (pvit(i+1, l) -pvit(i+1, k)) * vort(D*e+ 0, k)) ;
                        secmem(f, l) -= fac_f * mu_f(f, l, c) * n_f(f, 0)/fs(f) * coeff(k, l) * ((pvit(i+0, l) -pvit(i+0, k)) * vort(D*e+ 1, k) - (pvit(i+1, l) -pvit(i+1, k)) * vort(D*e+ 0, k)) ;
                      } // 100% explicit

                  }
            }

        }
}
