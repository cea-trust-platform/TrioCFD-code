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
// File:        Correction_Lubchenko_PolyMAC_P0.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Correction_Lubchenko_PolyMAC_P0.h>
#include <Pb_Multiphase.h>
#include <Portance_interfaciale_PolyMAC_P0.h>
#include <Portance_interfaciale_base.h>
#include <Op_Diff_Turbulent_PolyMAC_P0_Face.h>
#include <Dispersion_bulles_PolyMAC_P0.h>
#include <Dispersion_bulles_base.h>
#include <Champ_Face_PolyMAC_P0.h>
#include <Milieu_composite.h>
#include <Navier_Stokes_std.h>
#include <Viscosite_turbulente_base.h>
#include <Viscosite_turbulente_multiple.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <math.h>
#include <Dispersion_bulles_turbulente_Burns.h>
#include <Domaine_PolyMAC_P0.h>
#include <Probleme_base.h>
#include <Champ_Elem_PolyMAC_P0.h>

Implemente_instanciable(Correction_Lubchenko_PolyMAC_P0, "Correction_Lubchenko_Face_PolyMAC_P0", Source_base);
// XD Correction_Lubchenko source_base Correction_Lubchenko 1 not_set
// XD attr beta_lift floattant beta_lift 1 not_set
// XD attr beta_disp floattant beta_disp 1 not_set

Sortie& Correction_Lubchenko_PolyMAC_P0::printOn(Sortie& os) const
{
  return os;
}

Entree& Correction_Lubchenko_PolyMAC_P0::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("beta_lift", &beta_lift_);
  param.ajouter("beta_disp", &beta_disp_);
  param.ajouter("portee_disp", &portee_disp_);
  param.ajouter("portee_lift", &portee_lift_);
  param.ajouter("use_bif", &use_bif_);
  param.lire_avec_accolades_depuis(is);

  //identification des phases
  Pb_Multiphase *pbm = sub_type(Pb_Multiphase, equation().probleme()) ? &ref_cast(Pb_Multiphase, equation().probleme()) : nullptr;

  if (!pbm || pbm->nb_phases() == 1) Process::exit(que_suis_je() + " : not needed for single-phase flow!");
  for (int n = 0; n < pbm->nb_phases(); n++) //recherche de n_l, n_g : phase {liquide,gaz}_continu en priorite
    if (pbm->nom_phase(n).debute_par("liquide") && (n_l < 0 || pbm->nom_phase(n).finit_par("continu")))  n_l = n;

  if (n_l < 0) Process::exit(que_suis_je() + " : liquid phase not found!");

  pbm->creer_champ("distance_paroi_globale"); // Besoin de distance a la paroi

  return is;
}

void Correction_Lubchenko_PolyMAC_P0::completer() // We must wait for all readOn's to be sure that the bubble dispersion and lift correlations are created
{
  const Pb_Multiphase& pbm = ref_cast(Pb_Multiphase, equation().probleme());

  for (int i = 0 ; i < equation().sources().size() ; i++)
    if sub_type(Portance_interfaciale_PolyMAC_P0, equation().sources()(i).valeur())
      correlation_lift_ = ref_cast(Portance_interfaciale_PolyMAC_P0, equation().sources()(i).valeur()).correlation();

  for (int i = 0 ; i < equation().sources().size() ; i++)
    if sub_type(Dispersion_bulles_PolyMAC_P0, equation().sources()(i).valeur())
      correlation_dispersion_ = ref_cast(Dispersion_bulles_PolyMAC_P0, equation().sources()(i).valeur()).correlation();

  if ( (!correlation_lift_.non_nul()) || (!correlation_dispersion_.non_nul()) )
    Process::exit("Correction_Lubchenko_PolyMAC_P0::completer() : a dispersion_bulles and a portance_interfaciale force must be defined !");

  if sub_type(Op_Diff_Turbulent_PolyMAC_P0_Face, equation().operateur(0).l_op_base()) is_turb = 1;

  if (!pbm.has_champ("diametre_bulles")) Process::exit("Correction_Lubchenko_PolyMAC_P0::completer() : a bubble diameter must be defined !");
}


void Correction_Lubchenko_PolyMAC_P0::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const // The necessary dimensionner_bloc is taken care of in the dispersion_bulles_PolyMAC_P0 and Portance_interfaciale_PolyMAC_P0 functions
{
}

void Correction_Lubchenko_PolyMAC_P0::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  ajouter_blocs_disp(matrices, secmem, semi_impl);
  ajouter_blocs_lift(matrices, secmem, semi_impl);
  // Le terme source_BIF peut entrainer une mauvaise prédiction des forces normales à la paroi : une correction est mise en place
  if ((use_bif_) && (sub_type(Viscosite_turbulente_multiple, ref_cast(Op_Diff_Turbulent_PolyMAC_P0_Face, ref_cast(Navier_Stokes_std, equation().probleme().equation(0)).operateur(0).l_op_base()).correlation())))
    ajouter_blocs_BIF(matrices, secmem, semi_impl);
}

void Correction_Lubchenko_PolyMAC_P0::ajouter_blocs_disp(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl ) const
{
  const Pb_Multiphase& pbm = ref_cast(Pb_Multiphase, equation().probleme());
  const Champ_Face_PolyMAC_P0& ch = ref_cast(Champ_Face_PolyMAC_P0, equation().inconnue());
  const Domaine_PolyMAC_P0& domaine = ref_cast(Domaine_PolyMAC_P0, equation().domaine_dis());
  const IntTab& f_e = domaine.face_voisins(), &fcl = ch.fcl();
  const DoubleVect& pe = equation().milieu().porosite_elem(), &pf = equation().milieu().porosite_face(), &ve = domaine.volumes(), &vf = domaine.volumes_entrelaces(), &fs = domaine.face_surfaces();
  const DoubleTab& vf_dir = domaine.volumes_entrelaces_dir(), &n_f = domaine.face_normales();
  const DoubleTab& pvit = ch.passe(),
                   &alpha = pbm.equation_masse().inconnue().passe(),
                    &press = ref_cast(QDM_Multiphase, pbm.equation_qdm()).pression().passe(),
                     &temp  = pbm.equation_energie().inconnue().passe(),
                      &rho   = equation().milieu().masse_volumique().passe(),
                       &mu    = ref_cast(Fluide_base, equation().milieu()).viscosite_dynamique().passe(),
                        &y_elem = domaine.y_elem(),
                         &y_faces = domaine.y_faces(),
                          &n_y_elem = domaine.normale_paroi_elem(),
                           &n_y_faces = domaine.normale_paroi_faces(),
                            &d_bulles = equation().probleme().get_champ("diametre_bulles").valeurs(),
                             *k_turb = (equation().probleme().has_champ("k")) ? &equation().probleme().get_champ("k").passe() : nullptr,
                              * k_WIT = (equation().probleme().has_champ("k_WIT")) ? &equation().probleme().get_champ("k_WIT").passe() : nullptr ;
  const Milieu_composite& milc = ref_cast(Milieu_composite, equation().milieu());

  int N = pvit.line_size() , Np = press.line_size(), Nk = (k_turb) ? (*k_turb).dimension(1) : 1, D = dimension,
      nf_tot = domaine.nb_faces_tot(), nf = domaine.nb_faces(), ne_tot = domaine.nb_elem_tot(),
      cR = (rho.dimension_tot(0) == 1), cM = (mu.dimension_tot(0) == 1);

  DoubleTrav nut(domaine.nb_elem_tot(), N); //viscosite turbulente
  if (is_turb)
    ref_cast(Viscosite_turbulente_base, ref_cast(Op_Diff_Turbulent_PolyMAC_P0_Face, equation().operateur(0).l_op_base()).correlation()).eddy_viscosity(nut); //remplissage par la correlation

  // Input-output
  const Dispersion_bulles_base& correlation_db = ref_cast(Dispersion_bulles_base, correlation_dispersion_.valeur());
  Dispersion_bulles_base::input_t in;
  Dispersion_bulles_base::output_t out;
  in.alpha.resize(N);
  in.T.resize(N);
  in.p.resize(N);
  in.rho.resize(N);
  in.mu.resize(N);
  in.sigma.resize(N*(N-1)/2);
  in.k_turb.resize(N);
  in.nut.resize(N);
  in.d_bulles.resize(N);
  in.nv.resize(N, N);
  out.Ctd.resize(N, N);

  // There is no need to calculate the gradient of alpha here

  // Et pour les methodes span de la classe Interface pour choper la tension de surface
  const int nb_max_sat =  N * (N-1) /2; // oui !! suite arithmetique !!
  DoubleTrav Sigma_tab(ne_tot,nb_max_sat);

  // remplir les tabs ...
  for (int k = 0; k < N; k++)
    for (int l = k + 1; l < N; l++)
      {
        if (milc.has_saturation(k, l))
          {
            Saturation_base& z_sat = milc.get_saturation(k, l);
            const int ind_trav = (k*(N-1)-(k-1)*(k)/2) + (l-k-1); // Et oui ! matrice triang sup !
            // recuperer sigma ...
            const DoubleTab& sig = z_sat.get_sigma_tab();
            // fill in the good case
            for (int ii = 0; ii < ne_tot; ii++) Sigma_tab(ii, ind_trav) = sig(ii);
          }
        else if (milc.has_interface(k, l))
          {
            Interface_base& sat = milc.get_interface(k,l);
            const int ind_trav = (k*(N-1)-(k-1)*(k)/2) + (l-k-1); // Et oui ! matrice triang sup !
            for (int i = 0 ; i<ne_tot ; i++) Sigma_tab(i,ind_trav) = sat.sigma(temp(i,k),press(i,k * (Np > 1))) ;
          }
      }

  /* faces */
  for (int f = 0; f < nf; f++)
    if (fcl(f, 0) < 2)
      {
        in.alpha=0., in.T=0., in.p=0., in.rho=0., in.mu=0., in.sigma=0., in.k_turb=0., in.nut=0., in.d_bulles=0., in.nv=0.;
        int e;
        for (int c = 0; c < 2 && (e = f_e(f, c)) >= 0; c++)
          {
            for (int n = 0; n < N; n++)
              {
                in.alpha[n] += vf_dir(f, c)/vf(f) * alpha(e, n);
                in.p[n]     += vf_dir(f, c)/vf(f) * press(e, n * (Np > 1));
                in.T[n]     += vf_dir(f, c)/vf(f) * temp(e, n);
                in.rho[n]   += vf_dir(f, c)/vf(f) * rho(!cR * e, n);
                in.mu[n]    += vf_dir(f, c)/vf(f) * mu(!cM * e, n);
                in.nut[n]   += is_turb    ? vf_dir(f, c)/vf(f) * nut(e,n) : 0;
                in.d_bulles[n] += vf_dir(f, c)/vf(f) * d_bulles(e,n);
                for (int k = n+1; k < N; k++)
                  if (milc.has_interface(n,k))
                    {
                      const int ind_trav = (n*(N-1)-(n-1)*(n)/2) + (k-n-1);
                      in.sigma[ind_trav] += vf_dir(f, c) / vf(f) * Sigma_tab(e, ind_trav);
                    }
                for (int k = 0; k < N; k++)
                  in.nv(k, n) += vf_dir(f, c)/vf(f) * ch.v_norm(pvit, pvit, e, f, k, n, nullptr, nullptr);
              }
            for (int n = 0; n <Nk; n++) in.k_turb[n]   += (k_turb)   ? vf_dir(f, c)/vf(f) * (*k_turb)(e,0) : 0;
            in.k_WIT   += (k_WIT)   ? vf_dir(f, c)/vf(f) * (*k_WIT)(e,0) : 0.;
            in.e = e;
          }

        correlation_db.coefficient(in, out);

        double sum_alphag_wall = 0 ;
        for (int k = 0; k<N ; k++)
          if (k!=n_l) sum_alphag_wall += (y_faces(f)<portee_disp_*in.d_bulles[k]/2.) ? in.alpha[k] * (portee_disp_*in.d_bulles[k]-2*y_faces(f))/(portee_disp_*in.d_bulles[k]-y_faces(f)) :0 ;

        for (int k = 0; k < N; k++)
          if (k != n_l)
            if (y_faces(f)<portee_disp_*in.d_bulles[k]/2.)
              {
                double fac = 0 ;
                for (int d = 0 ; d<D ; d++) fac += n_y_faces(f, d) * n_f(f, d)/fs(f);

                fac *= beta_disp_*pf(f) * vf(f) ;
                secmem(f, k)   += fac * out.Ctd(k, n_l) * 1/y_faces(f) * in.alpha[k] * (portee_disp_*in.d_bulles[k]-2*y_faces(f))/(portee_disp_*in.d_bulles[k]-y_faces(f));
                secmem(f, k)   += fac * out.Ctd(n_l, k) * 1/y_faces(f) * sum_alphag_wall;
                secmem(f, n_l) -= fac * out.Ctd(k, n_l) * 1/y_faces(f) * in.alpha[k] * (portee_disp_*in.d_bulles[k]-2*y_faces(f))/(portee_disp_*in.d_bulles[k]-y_faces(f));
                secmem(f, n_l) -= fac * out.Ctd(n_l, k) * 1/y_faces(f) * sum_alphag_wall;
              }
      }

  /* elements */
  for (int e = 0; e < ne_tot; e++)
    {
      /* arguments de coeff */
      for (int n = 0; n < N; n++)
        {
          in.alpha[n] = alpha(e, n);
          in.p[n]     = press(e, n * (Np > 1));
          in.T[n]     = temp(e, n);
          in.rho[n]   = rho(!cR * e, n);
          in.mu[n]    = mu(!cM * e, n);
          in.nut[n]   = is_turb    ? nut(e,n) : 0;
          in.d_bulles[n] = d_bulles(e,n) ;
          for (int k = n+1; k < N; k++)
            if (milc.has_interface(n,k))
              {
                const int ind_trav = (n*(N-1)-(n-1)*(n)/2) + (k-n-1);
                in.sigma[ind_trav] = Sigma_tab(e, ind_trav);
              }
          for (int k = 0; k < N; k++)
            in.nv(k, n) = ch.v_norm(pvit, pvit, e, -1, k, n, nullptr, nullptr);
        }
      for (int n = 0; n <Nk; n++) in.k_turb[n]   = (k_turb) ? (*k_turb)(e,0) : 0;
      in.k_WIT   = (k_WIT) ? (*k_WIT)(e,0) : 0;
      in.e = e;

      correlation_db.coefficient(in, out);

      double sum_alphag_wall = 0 ;
      for (int k = 0; k<N ; k++)
        if (k!=n_l) sum_alphag_wall += (y_elem(e)<portee_disp_*d_bulles(e,k)/2.) ? alpha(e,k) *(portee_disp_*d_bulles(e,k)-2*y_elem(e))/(portee_disp_*d_bulles(e,k)-y_elem(e)) :0 ;
      for (int d = 0, i = nf_tot + D * e; d < D; d++, i++)
        for (int k = 0; k < N; k++)
          if (k != n_l)
            if (y_elem(e)<portee_disp_*d_bulles(e,k)/2)
              {
                double fac = beta_disp_*pe(e) * ve(e);

                secmem(i, k)   += fac * out.Ctd(k, n_l) * 1/y_elem(e) * alpha(e,k) * (portee_disp_*d_bulles(e,k)-2*y_elem(e))/(portee_disp_*d_bulles(e,k)-y_elem(e)) * n_y_elem(e, d);
                secmem(i, k)   += fac * out.Ctd(n_l, k) * 1/y_elem(e) * sum_alphag_wall                                      * n_y_elem(e, d);
                secmem(i, n_l) -= fac * out.Ctd(k, n_l) * 1/y_elem(e) * alpha(e,k) * (portee_disp_*d_bulles(e,k)-2*y_elem(e))/(portee_disp_*d_bulles(e,k)-y_elem(e)) * n_y_elem(e, d);
                secmem(i, n_l) -= fac * out.Ctd(n_l, k) * 1/y_elem(e) * sum_alphag_wall                                      * n_y_elem(e, d);
              }
    }
}


void Correction_Lubchenko_PolyMAC_P0::ajouter_blocs_lift(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl ) const
{
  const Pb_Multiphase& pbm = ref_cast(Pb_Multiphase, equation().probleme());
  const Champ_Face_PolyMAC_P0& ch = ref_cast(Champ_Face_PolyMAC_P0, equation().inconnue());
  const Domaine_PolyMAC_P0& domaine = ref_cast(Domaine_PolyMAC_P0, equation().domaine_dis());
  const IntTab& f_e = domaine.face_voisins(), &fcl = ch.fcl(), &e_f = domaine.elem_faces();
  const DoubleVect& pe = equation().milieu().porosite_elem(), &pf = equation().milieu().porosite_face(), &ve = domaine.volumes(), &vf = domaine.volumes_entrelaces(), &fs = domaine.face_surfaces();
  const DoubleTab& vf_dir = domaine.volumes_entrelaces_dir(), &n_f = domaine.face_normales();
  const DoubleTab& pvit = ch.passe(),
                   &alpha = pbm.equation_masse().inconnue().passe(),
                    &press = ref_cast(QDM_Multiphase, pbm.equation_qdm()).pression().passe(),
                     &temp  = pbm.equation_energie().inconnue().passe(),
                      &rho   = equation().milieu().masse_volumique().passe(),
                       &mu    = ref_cast(Fluide_base, equation().milieu()).viscosite_dynamique().passe(),
                        &vort  = equation().probleme().get_champ("vorticite").valeurs(),
                         &y_elem = domaine.y_elem(),
                          &y_faces = domaine.y_faces(),
                           &d_bulles = equation().probleme().get_champ("diametre_bulles").valeurs(),
                            &grad_v = equation().probleme().get_champ("gradient_vitesse").valeurs(),
                             * k_turb = (equation().probleme().has_champ("k")) ? &equation().probleme().get_champ("k").passe() : nullptr ;

  const Milieu_composite& milc = ref_cast(Milieu_composite, equation().milieu());

  int N = pvit.line_size() , Np = press.line_size(), Nk = (k_turb) ? (*k_turb).dimension(1) : 1, D = dimension,
      nf_tot = domaine.nb_faces_tot(), cR = (rho.dimension_tot(0) == 1), cM = (mu.dimension_tot(0) == 1);

  DoubleTrav vr_l(N,D), scal_ur(N), scal_u(N), pvit_l(N, D), vort_l( D==2 ? 1 :D), grad_l(D,D), scal_grad(D); // Requis pour corrections vort et u_l-u-g

  const Portance_interfaciale_base& correlation_pi = ref_cast(Portance_interfaciale_base, correlation_lift_.valeur());
  double vl_norm ;

  Portance_interfaciale_base::input_t in;
  Portance_interfaciale_base::output_t out;
  in.alpha.resize(N), in.T.resize(N), in.p.resize(N), in.rho.resize(N), in.mu.resize(N), in.sigma.resize(N*(N-1)/2), in.k_turb.resize(N), in.d_bulles.resize(N), in.nv.resize(N, N);
  out.Cl.resize(N, N);

  // Et pour les methodes span de la classe Interface pour choper la tension de surface
  const int ne_tot = domaine.nb_elem_tot(), nb_max_sat =  N * (N-1) /2; // oui !! suite arithmetique !!
  DoubleTrav Sigma_tab(ne_tot,nb_max_sat);

  // remplir les tabs ...
  for (int k = 0; k < N; k++)
    for (int l = k + 1; l < N; l++)
      {
        if (milc.has_saturation(k, l))
          {
            Saturation_base& z_sat = milc.get_saturation(k, l);
            const int ind_trav = (k*(N-1)-(k-1)*(k)/2) + (l-k-1); // Et oui ! matrice triang sup !
            // XXX XXX XXX
            // Attention c'est dangereux ! on suppose pour le moment que le champ de pression a 1 comp. Par contre la taille de res est nb_max_sat*nbelem !!
            // Aussi, on passe le Span le nbelem pour le champ de pression et pas nbelem_tot ....
            assert(press.line_size() == 1);
            assert(temp.line_size() == N);
            z_sat.get_sigma(temp.get_span_tot(), press.get_span_tot(), Sigma_tab.get_span_tot(), nb_max_sat, ind_trav);
          }
        else if (milc.has_interface(k, l))
          {
            Interface_base& sat = milc.get_interface(k,l);
            const int ind_trav = (k*(N-1)-(k-1)*(k)/2) + (l-k-1); // Et oui ! matrice triang sup !
            for (int i = 0 ; i<ne_tot ; i++) Sigma_tab(i,ind_trav) = sat.sigma(temp(i,k),press(i,k * (Np > 1))) ;
          }
      }

  /* elements */
  int f;
  for (int e = 0; e < ne_tot; e++)
    {
      /* arguments de coeff */
      for (int n=0; n<N; n++)
        {
          in.alpha[n] = alpha(e, n);
          in.p[n]     = press(e, n * (Np > 1));
          in.T[n]     = temp(e, n);
          in.rho[n]   = rho(!cR * e, n);
          in.mu[n]    = mu(!cM * e, n);
          in.d_bulles[n] = d_bulles(e,n) ;
          for (int k = n+1; k < N; k++)
            if (milc.has_interface(n,k))
              {
                const int ind_trav = (n*(N-1)-(n-1)*(n)/2) + (k-n-1);
                in.sigma[ind_trav] = Sigma_tab(e, ind_trav);
              }
          for (int k = 0; k < N; k++)
            in.nv(k, n) = ch.v_norm(pvit, pvit, e, -1, k, n, nullptr, nullptr);
        }
      for (int n = 0; n <Nk; n++) in.k_turb[n]   = (k_turb) ? (*k_turb)(e,0) : 0;

      correlation_pi.coefficient(in, out);

      double fac_e = beta_lift_*pe(e) * ve(e);
      int i = nf_tot + D * e;

      // Experimentation sur la portance
      vl_norm = 0;
      scal_ur = 0 ;
      for (int d = 0 ; d < D ; d++) vl_norm += pvit(i+d, n_l)*pvit(i+d, n_l);
      vl_norm = std::sqrt(vl_norm);
      if (vl_norm > 1.e-6)
        {
          for (int k = 0; k < N; k++)
            for (int d = 0 ; d < D ; d++) scal_ur(k) += pvit(i+d, n_l)/vl_norm * (pvit(i+d, k) -pvit(i+d, n_l));
          for (int k = 0; k < N; k++)
            for (int d = 0 ; d < D ; d++) vr_l(k, d)  = pvit(i+d, n_l)/vl_norm * scal_ur(k) ;
        }
      else for (int k=0 ; k<N ; k++)
          for (int d=0 ; d<D ; d++) vr_l(k, d) = pvit(i+d, k) -pvit(i+d, n_l) ;

      if (D==2)
        {
          for (int k = 0; k < N; k++)
            if (k!= n_l) // gas phase
              {
                // Damping of the lift force close to the wall;
                if      (y_elem(e) < .5*d_bulles(e,k)*portee_lift_) fac_e *= -1 ; // suppresses lift
                else if (y_elem(e) >    d_bulles(e,k)*portee_lift_) fac_e *=  0 ; // no effect
                else   fac_e *= (3*std::pow(2*y_elem(e)/(d_bulles(e,k)*portee_lift_)-1, 2) - 2*std::pow(2*y_elem(e)/(d_bulles(e,k)*portee_lift_)-1, 3)) - 1; // partial damping

                secmem(i, n_l) += fac_e * out.Cl(n_l, k) * vr_l(k, 1) * vort(e, 0) ;
                secmem(i,  k ) -= fac_e * out.Cl(n_l, k) * vr_l(k, 1) * vort(e, 0) ;
                secmem(i+1,n_l)-= fac_e * out.Cl(n_l, k) * vr_l(k, 0) * vort(e, 0) ;
                secmem(i+1, k )+= fac_e * out.Cl(n_l, k) * vr_l(k, 0) * vort(e, 0) ;
              } // 100% explicit
          for (int b = 0; b < e_f.dimension(1) && (f = e_f(e, b)) >= 0; b++)
            if (f<domaine.nb_faces())
              if (fcl(f, 0) < 2)
                for (int k = 0; k < N; k++)
                  if (k!= n_l) // gas phase
                    {
                      int c = (e == f_e(f, 0)) ? 0 : 1 ;
                      double fac_f = beta_lift_*pf(f) * vf_dir(f, c);  // Coherence with portance_interfaciale that calculates the correlation at the element

                      if   (y_elem(e) < .5*d_bulles(e,k)*portee_lift_) fac_f *= -1 ; // suppresses lift
                      else if (y_elem(e) > d_bulles(e,k)*portee_lift_) fac_f *=  0 ; // no effect
                      else fac_f *= (3*std::pow(2*y_elem(e)/(d_bulles(e,k)*portee_lift_)-1, 2) - 2*std::pow(2*y_elem(e)/(d_bulles(e,k)*portee_lift_)-1, 3)) - 1; // partial damping

                      secmem(f, n_l) += fac_f * n_f(f, 0)/fs(f) * out.Cl(n_l, k) * vr_l(k, 1) * vort(e, 0) ;
                      secmem(f,  k ) -= fac_f * n_f(f, 0)/fs(f) * out.Cl(n_l, k) * vr_l(k, 1) * vort(e, 0) ;
                      secmem(f, n_l) -= fac_f * n_f(f, 1)/fs(f) * out.Cl(n_l, k) * vr_l(k, 0) * vort(e, 0) ;
                      secmem(f,  k ) += fac_f * n_f(f, 1)/fs(f) * out.Cl(n_l, k) * vr_l(k, 0) * vort(e, 0) ;
                    } // 100% explicit

        }

      if (D==3)
        {
          for (int k = 0; k < N; k++)
            if (k!= n_l) // gas phase
              {
                // Damping of the lift force close to the wall;
                if (y_elem(e) < .5*d_bulles(e,k)*portee_lift_) fac_e *= -1 ; // suppresses lift
                else if (y_elem(e) >    d_bulles(e,k)*portee_lift_) fac_e *=  0 ; // no effect
                else fac_e *= (3*std::pow(2*y_elem(e)/(d_bulles(e,k)*portee_lift_)-1, 2) - 2*std::pow(2*y_elem(e)/(d_bulles(e,k)*portee_lift_)-1, 3)) - 1; // partial damping

                secmem(i, n_l) += fac_e * out.Cl(n_l, k) * (vr_l(k, 1) * vort(e, 2*N+n_l) - vr_l(k, 2) * vort(e, 1*N+n_l)) ;
                secmem(i,  k ) -= fac_e * out.Cl(n_l, k) * (vr_l(k, 1) * vort(e, 2*N+n_l) - vr_l(k, 2) * vort(e, 1*N+n_l)) ;
                secmem(i+1,n_l)+= fac_e * out.Cl(n_l, k) * (vr_l(k, 2) * vort(e, 0*N+n_l) - vr_l(k, 0) * vort(e, 2*N+n_l)) ;
                secmem(i+1, k )-= fac_e * out.Cl(n_l, k) * (vr_l(k, 2) * vort(e, 0*N+n_l) - vr_l(k, 0) * vort(e, 2*N+n_l)) ;
                secmem(i+2,n_l)+= fac_e * out.Cl(n_l, k) * (vr_l(k, 0) * vort(e, 1*N+n_l) - vr_l(k, 1) * vort(e, 0*N+n_l)) ;
                secmem(i+2, k )-= fac_e * out.Cl(n_l, k) * (vr_l(k, 0) * vort(e, 1*N+n_l) - vr_l(k, 1) * vort(e, 0*N+n_l)) ;
              } // 100% explicit
        }

    }

  int c, e, n, k, d, d2;
  double fac_f ;

  if (D==3)
    for (f = 0 ; f<domaine.nb_faces() ; f++)
      if (fcl(f, 0) < 2)
        {
          in.alpha=0., in.T=0., in.p=0., in.rho=0., in.mu=0., in.sigma=0., in.k_turb=0., in.d_bulles=0., in.nv=0.;
          for ( c = 0; c < 2 && (e = f_e(f, c)) >= 0; c++)
            {
              for (n = 0; n < N; n++)
                {
                  in.alpha[n] += vf_dir(f, c)/vf(f) * alpha(e, n);
                  in.p[n]     += vf_dir(f, c)/vf(f) * press(e, n * (Np > 1));
                  in.T[n]     += vf_dir(f, c)/vf(f) * temp(e, n);
                  in.rho[n]   += vf_dir(f, c)/vf(f) * rho(!cR * e, n);
                  in.mu[n]    += vf_dir(f, c)/vf(f) * mu(!cM * e, n);
                  in.d_bulles[n] += vf_dir(f, c)/vf(f) *d_bulles(e,n) ;
                  for (k = n+1; k < N; k++)
                    if (milc.has_interface(n,k))
                      {
                        const int ind_trav = (n*(N-1)-(n-1)*(n)/2) + (k-n-1);
                        in.sigma[ind_trav] += vf_dir(f, c) / vf(f) * Sigma_tab(e, ind_trav);
                      }
                  for (k = 0; k < N; k++)
                    in.nv(k, n) += vf_dir(f, c)/vf(f) * ch.v_norm(pvit, pvit, e, f, k, n, nullptr, nullptr);
                }
              for (n = 0; n <Nk; n++) in.k_turb[n]   += (k_turb)   ? vf_dir(f, c)/vf(f) * (*k_turb)(e,0) : 0;
            }

          correlation_pi.coefficient(in, out);

          grad_l = 0; // we fill grad_l so that grad_l(d, d2) = du_d/dx_d2 by averaging between both elements
          for (d = 0 ; d<D ; d++)
            for (d2 = 0 ; d2<D ; d2++)
              for (c=0 ; c<2  && (e = f_e(f, c)) >= 0; c++)
                grad_l(d, d2) += vf_dir(f, c)/vf(f)*grad_v(nf_tot + D*e + d2 , n_l * D + d) ;
          //We replace the n_l components by the one calculated without interpolation to elements
          scal_grad = 0 ; // scal_grad(d) = grad(u_d).n_f
          for (d = 0 ; d<D ; d++)
            for (d2 = 0 ; d2<D ; d2++)
              scal_grad(d) += grad_l(d, d2)*n_f(f, d2)/fs(f);
          for (d = 0 ; d<D ; d++)
            for (d2 = 0 ; d2<D ; d2++)
              grad_l(d, d2) += (grad_v(f ,n_l*D+d) - scal_grad(d)) * n_f(f, d2)/fs(f);
          // We calculate the local vorticity using this local gradient
          vort_l(0) = grad_l(2, 1) - grad_l(1, 2); // dUz/dy - dUy/dz
          vort_l(1) = grad_l(0, 2) - grad_l(2, 0); // dUx/dz - dUz/dx
          vort_l(2) = grad_l(1, 0) - grad_l(0, 1); // dUy/dx - dUx/dy

          // We also need to calculate relative velocity at the face
          pvit_l = 0 ;
          for (d = 0 ; d<D ; d++)
            for (k = 0 ; k<N ; k++)
              for (c=0 ; c<2 && (e = f_e(f, c)) >= 0; c++)
                pvit_l(k, d) += vf_dir(f, c)/vf(f)*pvit(domaine.nb_faces_tot()+D*e+d, k) ;
          scal_u = 0;
          for (k = 0 ; k<N ; k++)
            for (d = 0 ; d<D ; d++)
              scal_u(k) += pvit_l(k, d)*n_f(f, d)/fs(f);
          for (k = 0 ; k<N ; k++)
            for (d = 0 ; d<D ; d++)
              pvit_l(k, d) += (pvit(f, k) - scal_u(k)) * n_f(f, d)/fs(f) ; // Corect velocity at the face
          vl_norm = 0;
          scal_ur = 0;
          for (d = 0 ; d < D ; d++) vl_norm += pvit_l(n_l, d)*pvit_l(n_l, d);
          vl_norm = std::sqrt(vl_norm);
          if (vl_norm > 1.e-6)
            {
              for (k = 0; k < N; k++)
                for (d = 0 ; d < D ; d++) scal_ur(k) += pvit_l(n_l, d)/vl_norm * (pvit_l(k, d) -pvit_l(n_l, d));
              for (k = 0; k < N; k++)
                for (d = 0 ; d < D ; d++) vr_l(k, d)  = pvit_l(n_l, d)/vl_norm * scal_ur(k) ;
            }
          else for (k=0 ; k<N ; k++)
              for (d=0 ; d<D ; d++) vr_l(k, d) = pvit_l(k, d)-pvit_l(n_l, d) ;


          // Use local vairables for the calculation of secmem

          for (k = 0; k < N; k++)
            if (k!= n_l) // gas phase
              {
                fac_f = beta_lift_*pf(f) * vf(f);

                if   (y_faces(f) < .5*in.d_bulles(k)*portee_lift_) fac_f *= -1 ; // suppresses lift
                else if (y_faces(f) > in.d_bulles(k)*portee_lift_) fac_f *=  0 ; // no effect
                else  fac_f *= (3*std::pow(2*y_faces(f)/(in.d_bulles(k)*portee_lift_)-1, 2) - 2*std::pow(2*y_faces(f)/(in.d_bulles(k)*portee_lift_)-1, 3)) - 1; // partial damping

                secmem(f, n_l) += fac_f * n_f(f, 0)/fs(f) * out.Cl(n_l, k) * (vr_l(k, 1) * vort_l(2) - vr_l(k, 2) * vort_l(1)) ;
                secmem(f,  k ) -= fac_f * n_f(f, 0)/fs(f) * out.Cl(n_l, k) * (vr_l(k, 1) * vort_l(2) - vr_l(k, 2) * vort_l(1)) ;
                secmem(f, n_l) += fac_f * n_f(f, 1)/fs(f) * out.Cl(n_l, k) * (vr_l(k, 2) * vort_l(0) - vr_l(k, 0) * vort_l(2)) ;
                secmem(f,  k ) -= fac_f * n_f(f, 1)/fs(f) * out.Cl(n_l, k) * (vr_l(k, 2) * vort_l(0) - vr_l(k, 0) * vort_l(2)) ;
                secmem(f, n_l) += fac_f * n_f(f, 2)/fs(f) * out.Cl(n_l, k) * (vr_l(k, 0) * vort_l(1) - vr_l(k, 1) * vort_l(0)) ;
                secmem(f,  k ) -= fac_f * n_f(f, 2)/fs(f) * out.Cl(n_l, k) * (vr_l(k, 0) * vort_l(1) - vr_l(k, 1) * vort_l(0)) ;
              } // 100% explicit
        }


}


void Correction_Lubchenko_PolyMAC_P0::ajouter_blocs_BIF(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl ) const
{
  const Domaine_PolyMAC_P0&                      domaine = ref_cast(Domaine_PolyMAC_P0, equation().domaine_dis());
  const Probleme_base&                       pb = ref_cast(Probleme_base, equation().probleme());
  const Navier_Stokes_std&               eq_qdm = ref_cast(Navier_Stokes_std, pb.equation(0));
  const Op_Diff_Turbulent_PolyMAC_P0_Face& Op_diff = ref_cast(Op_Diff_Turbulent_PolyMAC_P0_Face, eq_qdm.operateur(0).l_op_base());
  const DoubleTab&                      tab_rho = equation().probleme().get_champ("masse_volumique").passe();
  const DoubleTab&                      tab_alp = equation().probleme().get_champ("alpha").passe();
  const DoubleTab& d_bulles = equation().probleme().get_champ("diametre_bulles").valeurs();
  const DoubleTab& y_elem = domaine.y_elem(),
                   &y_faces = domaine.y_faces();

  const DoubleTab&                       vf_dir = domaine.volumes_entrelaces_dir(), &xp = domaine.xp(), &xv = domaine.xv();
  const DoubleVect& pe = equation().milieu().porosite_elem(), &ve = domaine.volumes(), &fs = domaine.face_surfaces(), &vf = domaine.volumes_entrelaces();
  const DoubleTab& normales_f = domaine.face_normales();
  const IntTab& voisins_f = domaine.face_voisins(), &e_f = domaine.elem_faces(), &f_e = domaine.face_voisins();

  const Viscosite_turbulente_multiple&    visc_turb = ref_cast(Viscosite_turbulente_multiple, Op_diff.correlation());

  int N = pb.get_champ("vitesse").valeurs().dimension(1), D = dimension, nf_tot = domaine.nb_faces_tot(), nf = domaine.nb_faces(), ne_tot = domaine.nb_elem_tot() ;

  // On recupere les tensions de reynolds des termes de BIF
  DoubleTrav Rij(0, N, D, D);
  MD_Vector_tools::creer_tableau_distribue(eq_qdm.pression().valeurs().get_md_vector(), Rij); //Necessary to compare size in reynolds_stress()
  visc_turb.reynolds_stress_BIF(Rij);

  DoubleTrav grad_Rij(0, N, D, D);
  MD_Vector_tools::creer_tableau_distribue(eq_qdm.vitesse().valeurs().get_md_vector(), grad_Rij); //Necessary to exchange virtual elements after calculation of the gradient at the faces()

  const Champ_Elem_PolyMAC_P0& ch_alpha = ref_cast(Champ_Elem_PolyMAC_P0, equation().probleme().get_champ("alpha"));	// Champ alpha qui servira à obtenir les coeffs du gradient ; normalement toujours des CAL de Neumann ; terme source qui n'apparait qu'en multiphase
  ch_alpha.init_grad(0); // Initialisation des tables fgrad_d, fgrad_e, fgrad_w qui dependent de la discretisation et du type de conditions aux limites --> pas de mises a jour necessaires
  IntTab& fg_d = ch_alpha.fgrad_d, fg_e = ch_alpha.fgrad_e; // Tables utilisees dans domaine_PolyMAC_P0::fgrad pour le calcul du gradient
  DoubleTab fg_w = ch_alpha.fgrad_w;

  // On calcule le gradient de Rij aux faces
  for (int n = 0; n < N; n++)
    for (int f = 0; f < nf_tot; f++)
      for (int d_i = 0; d_i <D ; d_i++)
        for (int d_j = 0; d_j <D ; d_j++)
          {
            grad_Rij(f, n, d_i, d_j) = 0;
            // grad_Rij(face, phase, x-coord, y-coord)   or
            // grad_Rij(nb_face_tot + dimension*element*gradient_component, phase, x-coord, y-coord)
            for (int j = fg_d(f); j < fg_d(f+1) ; j++)
              {
                int e = fg_e(j);
                int f_bord;
                if ( (f_bord = e - ne_tot) < 0) //contribution d'un element
                  grad_Rij(f, n, d_i, d_j) += fg_w(j) * Rij(e, n, d_i, d_j);
                else if ( (ch_alpha.fcl()(f_bord, 0) == 1) || (ch_alpha.fcl()(f_bord, 0) == 2)
                          || (ch_alpha.fcl()(f_bord, 0) == 3) || (ch_alpha.fcl()(f_bord, 0) == 6))
                  {
                    Process::exit("You must have a neumann limit condition on alpha for RIJ_BIF to work !");
                  }
              }
          }

  // On interpole le gradient de Rij aux elements
  for (int n = 0; n < N; n++)
    for (int e = 0; e < ne_tot; e++)
      for (int d_d=0 ; d_d<D ; d_d++) // on derive / d_d
        for (int d_i = 0; d_i <D ; d_i++)
          for (int d_j = 0; d_j <D ; d_j++)
            {
              grad_Rij(nf_tot + D *e + d_d, n, d_i,  d_j) = 0;
              for (int j = 0, f; j < e_f.dimension(1) && (f = e_f(e, j)) >= 0; j++)
                grad_Rij(nf_tot + D *e + d_d, n, d_i,  d_j) += (e == f_e(f, 0) ? 1 : -1) * fs(f) * (xv(f, d_d) - xp(e, d_d)) / ve(e) * grad_Rij(f, n, d_i,  d_j);
            }

  // On calcule le second membre aux elements

  for(int e = 0 ; e < ne_tot ; e++)
    for(int n = 0; n<N ; n++)
      for (int d_i = 0; d_i < D; d_i++)
        {
          // Find max value of d_bulles
          double d_b_elem = 0;
          for(int k = 0; k<N ; k++)
            {
              if (k!=n_l && d_bulles(e,k)>d_b_elem) d_b_elem = d_bulles(e,k);
            }
          // BIF correction in the near wall region
          if (y_elem(e) < 0.5 * d_b_elem)
            {
              double secmem_en = 0;
              for (int d_j = 0; d_j < D; d_j++)
                secmem_en -= grad_Rij(nf_tot + D *e + d_j, n, d_i,  d_j) ; // On annule la contribution de BIF sur la qdm lorsque la distance à la paroi y<0.5*d_bulles
              secmem_en *= (-1) * pe(e) * ve(e) * tab_alp(e, n) * tab_rho(e, n) ; // For us, Rij = < u_i u_j >, therefore *(-1)
              secmem(nf_tot + e*D + d_i, n) += secmem_en;
            }
        }


  // On calcule le second membre aux faces

  int e;

  for(int f = 0 ; f < nf ; f++)
    for(int n = 0; n<N ; n++)
      for (int i = 0; i < 2 && (e = voisins_f(f, i)) >= 0; i++)
        {
          // Find max value of d_bulles
          double d_b_face = 0;
          for(int k = 0; k<N ; k++)
            {
              if (k!=n_l && d_bulles(e,k)>d_b_face) d_b_face = d_bulles(e,k);
            }
          // BIF correction in the near wall region
          double d_b = vf_dir(f, i)/vf(f) * d_b_face ;
          if (y_faces(f) < 0.5 * d_b)
            {
              DoubleTrav secmem_en(3); // Contains the vector of the divergence of R_ij at the face
              secmem_en = 0;
              for (int d_i = 0; d_i < D; d_i++)
                for (int d_j = 0; d_j < D; d_j++)
                  {
                    secmem_en(d_i) -= grad_Rij(nf_tot + D *e + d_j, n, d_i,  d_j) ;
                  }

              for (int d_i = 0; d_i < D; d_i++)
                secmem_en(d_i) *= (-1) * pe(e) * vf_dir(f, i) * tab_alp(e, n) * tab_rho(e, n);// For us, Rij = < u_i u_j >, therefore *(-1)
              double flux_face = domaine.dot(&normales_f(f, 0), &secmem_en(0));
              secmem(f, n) += flux_face;
            }

        }

}
