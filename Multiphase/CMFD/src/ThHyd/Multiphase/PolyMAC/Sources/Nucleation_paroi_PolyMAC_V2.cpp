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
// File:        Nucleation_paroi_PolyMAC_V2.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Nucleation_paroi_PolyMAC_V2.h>
#include <Pb_Multiphase.h>
#include <Champ_P0_PolyMAC_V2.h>
#include <Op_Diff_PolyMAC_V2_Elem.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <math.h>

Implemente_instanciable(Nucleation_paroi_PolyMAC_V2, "Nucleation_paroi_PolyMAC_V2", Source_base);

Sortie& Nucleation_paroi_PolyMAC_V2::printOn(Sortie& os) const
{
  return os;
}

Entree& Nucleation_paroi_PolyMAC_V2::readOn(Entree& is)
{
  return is;
}

void Nucleation_paroi_PolyMAC_V2::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const // 100% explicit for now
{
}

void Nucleation_paroi_PolyMAC_V2::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Zone_PolyMAC_V2& zone = ref_cast(Zone_PolyMAC_V2, equation().zone_dis().valeur());
  const DoubleVect& pe = zone.porosite_elem(), &ve = zone.volumes();
  const IntTab& f_e = zone.face_voisins();

  const DoubleTab& rho = pbm.milieu().masse_volumique().passe(),
                   &press = pbm.eq_qdm.pression().passe(),
                    &qpi = ref_cast(Op_Diff_PolyMAC_V2_Elem, equation().probleme().operateur(0).l_op_base()).q_pi(),
                     &dnuc = ref_cast(Op_Diff_PolyMAC_V2_Elem, equation().probleme().operateur(0).l_op_base()).d_nucleation();

  int N = ref_cast(Pb_Multiphase, equation().probleme()).nb_phases(), Np = equation().probleme().get_champ("pression").valeurs().line_size();

  const Milieu_composite& milc = ref_cast(Milieu_composite, equation().probleme().milieu());

  /* elements */
  int n_l = 0 ; // phase porteuse
  for (int f = 0; f < zone.premiere_face_int(); f++) // On n'injecte que dans les elems de bord
    for (int k = 0, mp = 0 ; k<N ; k++ , mp += (Np > 1)) if (k != n_l) if (milc.has_saturation(n_l, k)) //phase vapeur
          {
            Saturation_base& sat = milc.get_saturation(n_l, k);
            int e = f_e(f, 0);
            double fac = pe(e) * ve(e) ;
            secmem(e , k) += fac * 6. * qpi(e, n_l, k) / ( dnuc(k) * rho(e, k) * sat.Lvap(press(e, mp))) ;
          }

}