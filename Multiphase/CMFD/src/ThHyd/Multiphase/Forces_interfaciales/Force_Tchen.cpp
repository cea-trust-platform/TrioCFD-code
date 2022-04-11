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

#include <Force_Tchen.h>
#include <Pb_Multiphase.h>
#include <Champ_Face_CoviMAC.h>
#include <Op_Grad_CoviMAC_Face.h>
#include <math.h>

Implemente_instanciable(Force_Tchen, "Force_Tchen|Tchen_force", Source_base);

Sortie& Force_Tchen::printOn(Sortie& os) const
{
  return os;
}

Entree& Force_Tchen::readOn(Entree& is)
{
  return is;
}

void Force_Tchen::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
}//Vitesse locale et pas de droit de deriver QDM/P,T,alpha

void Force_Tchen::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Champ_Face_CoviMAC& ch = ref_cast(Champ_Face_CoviMAC, equation().inconnue().valeur());
  Matrice_Morse *mat = matrices.count(ch.le_nom().getString()) ? matrices.at(ch.le_nom().getString()) : NULL;

  const Zone_CoviMAC& zone = ref_cast(Zone_CoviMAC, equation().zone_dis().valeur());
  const IntTab& f_e = zone.face_voisins(), &fcl = ch.fcl();
  const DoubleVect& pe = zone.porosite_elem(), &pf = zone.porosite_face(), &ve = zone.volumes(), &vf = zone.volumes_entrelaces();

  const DoubleTab& inco = ch.valeurs(), &pvit = ch.passe(), &mu_f = ref_cast(Op_Grad_CoviMAC_Face, ref_cast(Navier_Stokes_std, equation()).operateur_gradient().valeur()).mu_f(),
                   &alpha = ref_cast(Pb_Multiphase, equation().probleme()).eq_masse.inconnue().passe(),
                    &rho   = equation().milieu().masse_volumique().passe();

  double pas_tps = equation().probleme().schema_temps().pas_de_temps();
  int N = inco.line_size(), D = dimension, nf_tot = zone.nb_faces_tot(), nf = zone.nb_faces();

  int k = 0 ; //Phase porteuse

  /* elements */
  for (int e = 0; e < zone.nb_elem(); e++)
    for (int d = 0 ; d <D ; d++) for (int l = 1 ; l<N ; l++) //phase secondaire
        {
          double fac = pe(e) * ve(e) * alpha(e, l) * rho(e, k) ;

          secmem(nf_tot + D * e + d, l) += fac * (inco(nf_tot + D * e + d, k)-pvit(nf_tot + D * e + d, k))/pas_tps;
          secmem(nf_tot + D * e + d, k) -= fac * (inco(nf_tot + D * e + d, k)-pvit(nf_tot + D * e + d, k))/pas_tps;
          if (mat)
            {
              (*mat)( N *(nf_tot + D * e + d) + l, N *(nf_tot + D * e + d) + k) -= fac/pas_tps ;
              (*mat)( N *(nf_tot + D * e + d) + k, N *(nf_tot + D * e + d) + k) += fac/pas_tps ;
            }
        }

  /* faces */
  for (int f = 0; f < nf; f++) if (fcl(f, 0) < 2)
      for (int l = 1 ; l<N ; l++) //phase secondaire
        {
          double alpha_loc = alpha(f_e(f, 0), l) * mu_f(f, 0) + alpha(f_e(f, 1), l) * mu_f(f, 1);
          double rho_loc = rho(f_e(f, 0), k) * mu_f(f, 0) + rho(f_e(f, 1), k) * mu_f(f, 1);
          double fac = pf(f) * vf(f) * alpha_loc * rho_loc ;

          secmem(f, l) += fac * (inco(f, k)-pvit(f, k))/pas_tps;
          secmem(f, k) -= fac * (inco(f, k)-pvit(f, k))/pas_tps;
          if (mat)
            {
              (*mat)( f + l, f + k) -= fac/pas_tps ;
              (*mat)( f + k, f + k) += fac/pas_tps ;
            }
        }

}