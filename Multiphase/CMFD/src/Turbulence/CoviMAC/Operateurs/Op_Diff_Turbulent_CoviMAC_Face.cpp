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
// File:        Op_Diff_CoviMAC_Face.cpp
// Directory:   $TRUST_ROOT/src/CoviMAC/Operateurs
// Version:     1
//
//////////////////////////////////////////////////////////////////////////////

#include <Op_Diff_Turbulent_CoviMAC_Face.h>
#include <Pb_Multiphase.h>
#include <Viscosite_turbulente_base.h>

Implemente_instanciable( Op_Diff_Turbulent_CoviMAC_Face, "Op_Diff_Turbulent_CoviMAC_Face|Op_Diff_Turbulente_CoviMAC_Face", Op_Diff_CoviMAC_Face ) ;

Sortie& Op_Diff_Turbulent_CoviMAC_Face::printOn( Sortie& os ) const
{
  Op_Diff_CoviMAC_base::printOn( os );
  return os;
}

Entree& Op_Diff_Turbulent_CoviMAC_Face::readOn( Entree& is )
{
  //lecture de la correlation de viscosite turbulente
  corr.typer_lire(equation().probleme(), "viscosite_turbulente", is);
  return is;
}

void Op_Diff_Turbulent_CoviMAC_Face::modifier_nu(DoubleTab& mu) const
{
  if (!corr.non_nul()) return; //rien a faire
  const DoubleTab& rho = equation().milieu().masse_volumique().passe(),
                   *alpha = sub_type(Pb_Multiphase, equation().probleme()) ? &ref_cast(Pb_Multiphase, equation().probleme()).eq_masse.inconnue().passe() : NULL;
  int i, nl = mu.dimension(0), n, N = equation().inconnue().valeurs().line_size(), cR = rho.dimension(0) == 1, d, D = dimension;
  DoubleTrav nu_t(nl, N); //viscosite turbulente : toujours scalaire
  ref_cast(Viscosite_turbulente_base, corr.valeur()).eddy_viscosity(nu_t); //remplissage par la correlation
  if (mu.nb_dim() == 2) //nu scalaire
    for (i = 0; i < nl; i++) for (n = 0; n < N; n++) mu(i, n) += (alpha ? (*alpha)(i, n) : 1) * rho(!cR * i, n) * nu_t(i, n);
  else if (mu.nb_dim() == 3) //nu anisotrope diagonal
    for (i = 0; i < nl; i++) for (n = 0; n < N; n++) for (d = 0; d < D; d++) mu(i, n, d) += (alpha ? (*alpha)(i, n) : 1) * rho(!cR * i, n) * nu_t(i, n);
  else if (mu.nb_dim() == 4) //nu anisotrope complet
    for (i = 0; i < nl; i++) for (n = 0; n < N; n++) for (d = 0; d < D; d++) mu(i, n, d, d) += (alpha ? (*alpha)(i, n) : 1) * rho(!cR * i, n) * nu_t(i, n);
  else abort();
}