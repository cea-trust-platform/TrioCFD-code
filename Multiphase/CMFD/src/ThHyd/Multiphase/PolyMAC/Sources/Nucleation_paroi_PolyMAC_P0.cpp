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
// File:        Nucleation_paroi_PolyMAC_P0.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Nucleation_paroi_PolyMAC_P0.h>
#include <Pb_Multiphase.h>
#include <Champ_Elem_PolyMAC_P0.h>
#include <Op_Diff_PolyMAC_P0_Elem.h>
#include <Milieu_composite.h>
#include <Flux_parietal_base.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <math.h>

Implemente_instanciable(Nucleation_paroi_PolyMAC_P0, "Nucleation_paroi_Elem_PolyMAC_P0", Source_base);

Sortie& Nucleation_paroi_PolyMAC_P0::printOn(Sortie& os) const
{
  return os;
}

Entree& Nucleation_paroi_PolyMAC_P0::readOn(Entree& is)
{
  Pb_Multiphase *pbm = sub_type(Pb_Multiphase, equation().probleme()) ? &ref_cast(Pb_Multiphase, equation().probleme()) : NULL;

  if (!pbm || pbm->nb_phases() == 1) Process::exit(que_suis_je() + " : not needed for single-phase flow!");
  for (int n = 0; n < pbm->nb_phases(); n++) //recherche de n_l, n_g : phase {liquide,gaz}_continu en priorite
    if (pbm->nom_phase(n).debute_par("liquide") && (n_l < 0 || pbm->nom_phase(n).finit_par("continu")))  n_l = n;
  if (n_l < 0) Process::exit(que_suis_je() + " : liquid phase not found!");
  if (n_l != 0) Process::exit(que_suis_je() + " : liquid phase must be the first declared phase !");

  if (!pbm->has_correlation("flux_parietal")) Process::exit("Nucleation_paroi_PolyMAC_P0 : wall heat flux correlation needed !");
  const Flux_parietal_base& correlation_fp = ref_cast(Flux_parietal_base, pbm->get_correlation("flux_parietal").valeur());

  if (!correlation_fp.calculates_bubble_nucleation_diameter()) Process::exit("Nucleation_paroi_PolyMAC_P0 : wall heat flux correlation must calculate the nucleated bubble diameter !");

  return is;
}

void Nucleation_paroi_PolyMAC_P0::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const // 100% explicit for now
{
}

void Nucleation_paroi_PolyMAC_P0::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Pb_Multiphase& pbm = ref_cast(Pb_Multiphase, equation().probleme());
  const Zone_PolyMAC_P0& zone = ref_cast(Zone_PolyMAC_P0, equation().zone_dis().valeur());
  const DoubleVect& pe = zone.porosite_elem(), &ve = zone.volumes();
  const IntTab& f_e = zone.face_voisins();

  const DoubleTab& rho = pbm.milieu().masse_volumique().passe(),
                   &press = pbm.eq_qdm.pression().passe(),
                    &qpi = ref_cast(Op_Diff_PolyMAC_P0_Elem, pbm.equation(2).operateur(0).l_op_base()).q_pi(),
                     &dnuc = ref_cast(Op_Diff_PolyMAC_P0_Elem, pbm.equation(2).operateur(0).l_op_base()).d_nucleation();

  int N = pbm.nb_phases(), Np = pbm.get_champ("pression").valeurs().line_size();

  const Milieu_composite& milc = ref_cast(Milieu_composite, pbm.milieu());

  /* elements */
  for (int f = 0; f < zone.premiere_face_int(); f++) // On n'injecte que dans les elems de bord
    for (int k = 0, mp = 0 ; k<N ; k++ , mp += (Np > 1))
      if (k != n_l)
        if (milc.has_saturation(n_l, k)) //phase vapeur
          {
            Saturation_base& sat = milc.get_saturation(n_l, k);
            int e = f_e(f, 0);
            double fac = pe(e) * ve(e) ;
            secmem(e , k) += fac * 6. * qpi(e, n_l, k) / ( std::max(dnuc(e,k), 1.e-8) * rho(e, k) * sat.Lvap(press(e, mp))) ;
          }

}
