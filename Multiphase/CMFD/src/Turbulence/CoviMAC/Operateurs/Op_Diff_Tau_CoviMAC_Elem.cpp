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

void Op_Diff_Tau_CoviMAC_Elem::modifier_nu(DoubleTab& mu) const
{
  Op_Diff_Turbulent_CoviMAC_Elem::modifier_nu(mu);
  const DoubleTab& tau = equation().probleme().get_champ("tau").passe();
  int i, nl = mu.dimension_tot(0), n, N = mu.dimension(1), d, D = dimension;
  // Ici on divise mu par 1/tau^2
  if (mu.nb_dim() == 2) for (i = 0; i < nl; i++) for (n = 0; n < N; n++) //isotrope
        mu(i, n) *= (tau(i,n) > limiter_tau_) ? limiter_tau_*limiter_tau_/(tau(i,n)*tau(i,n)) : 1;
  else if (mu.nb_dim() == 3) for (i = 0; i < nl; i++) for (n = 0; n < N; n++) for (d = 0; d < D; d++) //anisotrope diagonal
          mu(i, n, d) *= (tau(i,n) > limiter_tau_) ? limiter_tau_*limiter_tau_/(tau(i,n)*tau(i,n)) : 1;
  else for (i = 0; i < nl; i++) for (n = 0; n < N; n++) for (d = 0; d < D; d++) //anisotrope complet
          mu(i, n, d, d) *= (tau(i,n) > limiter_tau_) ? limiter_tau_*limiter_tau_/(tau(i,n)*tau(i,n)) : 1;
}

void Op_Diff_Tau_CoviMAC_Elem::remodifier_nu(DoubleTab& mu) const
{
  const DoubleTab& tau = equation().probleme().get_champ("tau").passe();
  int i, nl = mu.dimension_tot(0), n, N = mu.dimension(1), d, D = dimension;
  // On annule l'effet du précédent modifier_nu
  if (mu.nb_dim() == 2) for (i = 0; i < nl; i++) for (n = 0; n < N; n++) //isotrope
        mu(i, n) /= (tau(i,n) > limiter_tau_) ? limiter_tau_*limiter_tau_/(tau(i,n)*tau(i,n)) : 1;
  else if (mu.nb_dim() == 3) for (i = 0; i < nl; i++) for (n = 0; n < N; n++) for (d = 0; d < D; d++) //anisotrope diagonal
          mu(i, n, d) /= (tau(i,n) > limiter_tau_) ? limiter_tau_*limiter_tau_/(tau(i,n)*tau(i,n)) : 1;
  else for (i = 0; i < nl; i++) for (n = 0; n < N; n++) for (d = 0; d < D; d++) //anisotrope complet
          mu(i, n, d, d) /= (tau(i,n) > limiter_tau_) ? limiter_tau_*limiter_tau_/(tau(i,n)*tau(i,n)) : 1;
  // On prepare le suivant
  if (mu.nb_dim() == 2) for (i = 0; i < nl; i++) for (n = 0; n < N; n++) //isotrope
        mu(i, n) *= (tau(i,n) > limiter_tau_) ? limiter_tau_*limiter_tau_*limiter_tau_/(tau(i,n)*tau(i,n)*tau(i,n)) : 1;
  else if (mu.nb_dim() == 3) for (i = 0; i < nl; i++) for (n = 0; n < N; n++) for (d = 0; d < D; d++) //anisotrope diagonal
          mu(i, n, d) *= (tau(i,n) > limiter_tau_) ? limiter_tau_*limiter_tau_*limiter_tau_/(tau(i,n)*tau(i,n)*tau(i,n)) : 1;
  else for (i = 0; i < nl; i++) for (n = 0; n < N; n++) for (d = 0; d < D; d++) //anisotrope complet
          mu(i, n, d, d) *= (tau(i,n) > limiter_tau_) ? limiter_tau_*limiter_tau_*limiter_tau_/(tau(i,n)*tau(i,n)*tau(i,n)) : 1;
}

void Op_Diff_Tau_CoviMAC_Elem::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  Op_Diff_Turbulent_CoviMAC_Elem::dimensionner_blocs(matrices, semi_impl);
}

void Op_Diff_Tau_CoviMAC_Elem::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Zone_CoviMAC&                      zone = ref_cast(Zone_CoviMAC, equation().zone_dis().valeur());
  const Champ_P0_CoviMAC&                   tau = ref_cast(Champ_P0_CoviMAC, equation().inconnue().valeur());
//  const DoubleTab&                      tab_tau = tau.passe();
  const DoubleTab&                      tab_tau = tau.valeurs();

  DoubleTrav sec_m(tab_tau); //residus
  matrices_t mat_m; //derivees vides
  Op_Diff_Turbulent_CoviMAC_Elem::ajouter_blocs(mat_m, sec_m, semi_impl); // Faire quelque chose si non semi implicite

  int N = tab_tau.dimension(1), ne = zone.nb_elem() ;
  Matrice_Morse *Mtau = (matrices.count("tau") && !semi_impl.count("tau")) ? matrices.at("tau") : nullptr;

  for(int e = 0 ; e < ne ; e++) for (int n=0 ; n<N ; n++)
      {
        double secmem_en = sec_m(e, n) ;
        secmem_en *= (tab_tau(e,n) > limiter_tau_) ? tab_tau(e,n)*tab_tau(e,n)/(limiter_tau_*limiter_tau_) : 1;
        secmem(e, n) += secmem_en;
        if (!(Mtau==nullptr))
          (*Mtau)(N*e+n, N*e+n) -= 2*sec_m(e, n)*((tab_tau(e,n) > limiter_tau_) ? tab_tau(e,n)/(limiter_tau_*limiter_tau_) : 0);
      }

  DoubleTrav sec_m_2(tab_tau); //residus
  matrices_t mat_m_2; //derivees vides
  Op_Diff_Turbulent_CoviMAC_Elem::ajouter_blocs(mat_m_2, sec_m_2, semi_impl); // Faire quelque chose si non semi implicite

  remodifier_nu(nu_);
  for(int e = 0 ; e < ne ; e++) for (int n=0 ; n<N ; n++)
      if (!(Mtau==nullptr))
        (*Mtau)(N*e+n, N*e+n) -= -2*sec_m_2(e, n)*((tab_tau(e,n) > limiter_tau_) ? tab_tau(e,n)*tab_tau(e,n)/(limiter_tau_*limiter_tau_*limiter_tau_) : 0);

}