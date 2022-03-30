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
// File:        Op_Diff_Tau_CoviMAC_Elem.cpp
// Directory:   $TRUST_ROOT/src/CoviMAC/Operateurs
// Version:     1
//
//////////////////////////////////////////////////////////////////////////////

#include <Op_Diff_Tau_CoviMAC_Elem.h>
#include <Op_Diff_Turbulent_CoviMAC_Face.h>
#include <Pb_Multiphase.h>
#include <Viscosite_turbulente_base.h>
#include <Transport_turbulent_base.h>
#include <Champ_P0_CoviMAC.h>
#include <Zone_CoviMAC.h>
#include <Neumann_paroi.h>
#include <cmath>
#include <functional>

Implemente_instanciable( Op_Diff_Tau_CoviMAC_Elem, "Op_Diff_Tau_CoviMAC_Elem", Op_Diff_Turbulent_CoviMAC_Elem ) ;

Sortie& Op_Diff_Tau_CoviMAC_Elem::printOn( Sortie& os ) const
{
  Op_Diff_Turbulent_CoviMAC_Elem::printOn( os );
  return os;
}

Entree& Op_Diff_Tau_CoviMAC_Elem::readOn( Entree& is )
{
  //lecture de la correlation de viscosite turbulente
  Op_Diff_Turbulent_CoviMAC_Elem::readOn( is );
  return is;
}

void Op_Diff_Tau_CoviMAC_Elem::completer()
{
  Op_Diff_Turbulent_CoviMAC_Elem::completer();
}

void Op_Diff_Tau_CoviMAC_Elem::modifier_nu(DoubleTab& mu) const // Multiplication par 1/tau^2 pour la diffusion
{
  Op_Diff_Turbulent_CoviMAC_Elem::modifier_nu(mu);
  const DoubleTab& tau = equation().probleme().get_champ("tau").passe();
  int i, nl = mu.dimension_tot(0), n, N = mu.dimension(1), D = dimension;
  // Ici on divise mu par 1/tau^2
  if (mu.nb_dim() == 2) for (i = 0; i < nl; i++) for (n = 0; n < N; n++) //isotrope
        {
          mu(i, n) *= ((tau(i,n) > limiter_tau_) ? limiter_tau_/(tau(i,n)*tau(i,n)) : 1/limiter_tau_);
        }
  else if (mu.nb_dim() == 3) for (i = 0; i < nl; i++) for (n = 0; n < N; n++) for (int d = 0; d < D; d++) //anisotrope diagonal
          mu(i, n, d) *= ((tau(i,n) > limiter_tau_) ? limiter_tau_/(tau(i,n)*tau(i,n)) : 1/limiter_tau_);
  else for (i = 0; i < nl; i++) for (n = 0; n < N; n++) for (int d1 = 0; d1 < D; d1++) for (int d2 = 0; d2 < D; d2++) //anisotrope complet
            mu(i, n, d1, d2) *= ((tau(i,n) > limiter_tau_) ? limiter_tau_/(tau(i,n)*tau(i,n)) : 1/limiter_tau_);
  mu.echange_espace_virtuel();
}

double Op_Diff_Tau_CoviMAC_Elem::calculer_dt_stab() const
{
  const Zone_CoviMAC& zone = la_zone_poly_.valeur();
  const IntTab& e_f = zone.elem_faces();
  const DoubleTab& nf = zone.face_normales(),
                   *alp = sub_type(Pb_Multiphase, equation().probleme()) ? &ref_cast(Pb_Multiphase, equation().probleme()).eq_masse.inconnue().passe() : NULL,
                    &diffu = diffusivite_pour_pas_de_temps().valeurs(), &lambda = diffusivite().valeurs();
  const DoubleTab& tau = equation().probleme().get_champ("tau").passe();
  const DoubleVect& pe = zone.porosite_elem(), &vf = zone.volumes_entrelaces(), &ve = zone.volumes();
  update_nu();
  int i, d, e, f, n, N = equation().inconnue().valeurs().dimension(1), cD = diffu.dimension(0) == 1, cL = lambda.dimension(0) == 1, D = dimension;

  DoubleTrav nu_loc(nu_); //nu local sans le 1/tau^2 pour stabilite
  if (nu_.nb_dim() == 2) for (e = 0; e < zone.nb_elem(); e++) for (n = 0; n < N; n++) //isotrope
        nu_loc(e, n) = nu_(e,n) / ((tau(e,n) > limiter_tau_) ? limiter_tau_/(tau(e,n)*tau(e,n)) : 1/limiter_tau_);
  else if (nu_.nb_dim() == 3) for (e = 0; e < zone.nb_elem(); e++) for (n = 0; n < N; n++) for (d = 0; d < D; d++) //anisotrope diagonal
          nu_loc(e, n, d) = nu_(e,n, d) / ((tau(e,n) > limiter_tau_) ? limiter_tau_/(tau(e,n)*tau(e,n)) : 1/limiter_tau_);
  else for (e = 0; e < zone.nb_elem(); e++) for (n = 0; n < N; n++) for (d = 0; d < D; d++) //anisotrope complet
          nu_loc(e, n, d, d) = nu_(e,n, d, d) / ((tau(e,n) > limiter_tau_) ? limiter_tau_/(tau(e,n)*tau(e,n)) : 1/limiter_tau_);

  double dt = 1e10;
  DoubleTrav flux(N);
  for (e = 0; e < zone.nb_elem(); e++)
    {
      for (flux = 0, i = 0; i < e_f.dimension(1) && (f = e_f(e, i)) >= 0; i++) for (n = 0; n < N; n++)
          flux(n) += zone.nu_dot(&nu_loc, e, n, &nf(f, 0), &nf(f, 0)) / vf(f);
      for (n = 0; n < N; n++) if ((!alp || (*alp)(e, n) > 1e-3) && flux(n)) /* sous 0.5e-6, on suppose que l'evanescence fait le job */
          dt = std::min(dt, pe(e) * ve(e) * (alp ? (*alp)(e, n) : 1) * (lambda(!cL * e, n) / diffu(!cD * e, n)) / flux(n));
      if (dt < 0) abort();
    }
  return Process::mp_min(dt);

}


void Op_Diff_Tau_CoviMAC_Elem::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  Op_Diff_Turbulent_CoviMAC_Elem::dimensionner_blocs(matrices, semi_impl);
}

void Op_Diff_Tau_CoviMAC_Elem::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const DoubleTab& tau = equation().probleme().get_champ("tau").passe();

  update_phif();
  const std::string& nom_inco = equation().inconnue().le_nom().getString();
  int n_ext = 1; // Ici on n'a qu'un seul milieu
  int i, j, e, eb, f, fb, n, M;
  std::vector<Matrice_Morse *> mat(n_ext); //matrices
  std::vector<int> N; //composantes
  std::vector<std::reference_wrapper<const Zone_CoviMAC>> zone; //zones
  std::vector<std::reference_wrapper<const Conds_lim>> cls; //conditions aux limites
  std::vector<std::reference_wrapper<const IntTab>> fcl, e_f, f_e, f_s; //tableaux "fcl", "elem_faces", "faces_voisins"
  std::vector<std::reference_wrapper<const DoubleVect>> fs; //surfaces
  std::vector<std::reference_wrapper<const DoubleTab>> inco, nf, xp, xs, xv; //inconnues, normales aux faces, positions elems / faces / sommets
  for (i = 0, M = 0; i < n_ext; M = std::max(M, N[i]), i++)
    {
      std::string nom_mat = nom_inco;
      mat[i] = !semi_impl.count(nom_inco) && matrices.count(nom_mat) ? matrices.at(nom_mat) : NULL;
      zone.push_back(std::ref(ref_cast(Zone_CoviMAC, equation().zone_dis().valeur())));
      f_e.push_back(std::ref(zone[i].get().face_voisins())), e_f.push_back(std::ref(zone[i].get().elem_faces())), f_s.push_back(std::ref(zone[i].get().face_sommets()));
      fs.push_back(std::ref(zone[i].get().face_surfaces())), nf.push_back(std::ref(zone[i].get().face_normales()));
      xp.push_back(std::ref(zone[i].get().xp())), xv.push_back(std::ref(zone[i].get().xv())), xs.push_back(std::ref(zone[i].get().zone().domaine().coord_sommets()));
      cls.push_back(std::ref(equation().zone_Cl_dis().les_conditions_limites()));
      const Champ_P0_CoviMAC& ch = ref_cast(Champ_P0_CoviMAC, equation().inconnue().valeur());
      inco.push_back(std::ref(semi_impl.count(nom_mat) ? semi_impl.at(nom_mat) : ch.valeurs()));
      N.push_back(inco[i].get().line_size()), fcl.push_back(std::ref(ch.fcl()));
    }
  const Zone_CoviMAC& zone0 = zone[0];

  /* avec phif : flux hors Echange_contact -> mat[0] seulement */
  DoubleTrav flux(N[0]);
  DoubleTrav secmem_loc(secmem);
  for (f = 0; f < zone0.nb_faces(); f++)
    {
      for (flux = 0, i = phif_d(f); i < phif_d(f + 1); i++)
        {
          if ((fb = (eb = phif_e(i)) - zone0.nb_elem_tot()) < 0) //element
            {
              for (n = 0; n < N[0]; n++) flux(n) += phif_c(i, n) * fs[0](f) * tau(eb, n);
              if (mat[0]) for (j = 0; j < 2 && (e = f_e[0](f, j)) >= 0; j++) if (e < zone[0].get().nb_elem()) for (n = 0; n < N[0]; n++) //derivees
                      (*mat[0])(N[0] * e + n, N[0] * eb + n) += (j ? 1 : -1) * phif_c(i, n) * fs[0](f) / ((tau(e,n) > limiter_tau_) ? limiter_tau_/(tau(e,n)*tau(e,n)) : 1/limiter_tau_) ;
            }
          else if (fcl[0](fb, 0) == 1 || fcl[0](fb, 0) == 2) for (n = 0; n < N[0]; n++) //Echange_impose_base
              flux(n) += (phif_c(i, n) ? phif_c(i, n) * fs[0](f) * ref_cast(Echange_impose_base, cls[0].get()[fcl[0](fb, 1)].valeur()).T_ext(fcl[0](fb, 2), n) : 0);
          else if (fcl[0](fb, 0) == 4) for (n = 0; n < N[0]; n++) //Neumann non homogene
              flux(n) += (phif_c(i, n) ? phif_c(i, n) * fs[0](f) * ref_cast(Neumann_paroi, cls[0].get()[fcl[0](fb, 1)].valeur()).flux_impose(fcl[0](fb, 2), n) : 0) * ((tau(f_s[0](f,0),n) > limiter_tau_) ? limiter_tau_/(tau(f_s[0](f,0),n)*tau(f_s[0](f,0),n)) : 1/limiter_tau_); // Pour compenser la correction de flux qui vient avec le tau**2 a la fin
          else if (fcl[0](fb, 0) == 6) for (n = 0; n < N[0]; n++) //Dirichlet
              flux(n) += (phif_c(i, n) ? phif_c(i, n) * fs[0](f) * ref_cast(Dirichlet, cls[0].get()[fcl[0](fb, 1)].valeur()).val_imp(fcl[0](fb, 2), n) : 0);
          else if (fcl[0](fb, 0) != 5) Cerr <<  "condslim" << fcl[0](fb, 0) ;
        }

      for (j = 0; j < 2 && (e = f_e[0](f, j)) >= 0; j++) if (e < zone[0].get().nb_elem()) for (n = 0; n < N[0]; n++) //second membre -> amont/aval
            secmem_loc(e, n) += (j ? -1 : 1) * flux(n) ;
    }
  Matrice_Morse* mat2 = (matrices.count("tau")) ? matrices.at("tau") : nullptr;

  for (e = 0 ; e<zone0.nb_elem() ; e++)  for (n = 0; n < N[0]; n++)
      {
        secmem(e,n) += secmem_loc(e,n) / ((tau(e,n) > limiter_tau_) ? limiter_tau_/(tau(e,n)*tau(e,n)) : 1/limiter_tau_);
        if (!(mat2==nullptr)) (*mat2)(N[0] * e + n, N[0] * e + n) -= 2 * secmem_loc(e,n) * ((tau(e,n) > limiter_tau_) ? tau(e,n)/limiter_tau_ : 0)  ;
      }

}