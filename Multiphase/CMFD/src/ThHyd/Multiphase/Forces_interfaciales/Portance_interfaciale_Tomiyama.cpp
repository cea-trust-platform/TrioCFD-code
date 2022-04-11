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
// File:        Portance_interfaciale_Tomiyama.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Portance_interfaciale_Tomiyama.h>
#include <Pb_Multiphase.h>
#include <math.h>

Implemente_instanciable(Portance_interfaciale_Tomiyama, "Portance_interfaciale_Tomiyama", Portance_interfaciale_base);

Sortie& Portance_interfaciale_Tomiyama::printOn(Sortie& os) const
{
  return os;
}

Entree& Portance_interfaciale_Tomiyama::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("constante_gravitation", &g_);
  param.lire_avec_accolades_depuis(is);
  return is;
}

void Portance_interfaciale_Tomiyama::coefficient(const DoubleTab& alpha, const DoubleTab& p, const DoubleTab& T,
                                                 const DoubleTab& rho, const DoubleTab& mu, const DoubleTab& sigma, double Dh,
                                                 const DoubleTab& ndv, int e, DoubleTab& coeff) const
{
  const DoubleTab& diametres = ref_cast(Pb_Multiphase, pb_.valeur()).get_champ("diametre_bulles").valeurs();
  int k, l, N = ndv.dimension(0);

  coeff = 0;
  k = 0; // Phase porteuse
  for (l = 1; l < N; l++)
    {
      double Re = rho(k) * ndv(k,l) * diametres(e, l)/mu(k);
      double Eo = g_ * std::abs(rho(k)-rho(l)) * diametres(e, l)*diametres(e, l)/sigma(k,l);
      double f_Eo = .00105*Eo*Eo*Eo - .0159*Eo*Eo - .0204*Eo + .474;
      double Cl;
      if (Eo<4) Cl = std::min( .288*std::tanh( .121*Re ), f_Eo) ;
      else      Cl = f_Eo ;

      coeff(l, k) = Cl * rho(k) * alpha(l) ;
      coeff(k, l) =  coeff(k, 0);
    }
}
