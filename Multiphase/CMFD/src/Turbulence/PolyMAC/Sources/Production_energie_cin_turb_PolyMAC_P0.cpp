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
// File:        Production_energie_cin_turb_PolyMAC_P0.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/PolyMAC_P0
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Production_energie_cin_turb_PolyMAC_P0.h>

#include <Zone_PolyMAC_P0.h>
#include <Champ_Elem_PolyMAC_P0.h>
#include <Matrix_tools.h>
#include <Probleme_base.h>
#include <grad_Champ_Face_PolyMAC_P0.h>
#include <Champ_Uniforme.h>
#include <Flux_interfacial_base.h>
#include <Milieu_composite.h>
#include <Operateur_Diff.h>
#include <Op_Diff_Turbulent_PolyMAC_P0_Face.h>
#include <Navier_Stokes_std.h>
#include <Viscosite_turbulente_base.h>
#include <TRUSTTab_parts.h>


Implemente_instanciable(Production_energie_cin_turb_PolyMAC_P0,"Production_energie_cin_turb_Elem_PolyMAC_P0", Source_base);
// XD Production_energie_cin_turb source_base Production_energie_cin_turb 0 Production source term for the TKE equation


Sortie& Production_energie_cin_turb_PolyMAC_P0::printOn(Sortie& os) const
{
  return os;
}

Entree& Production_energie_cin_turb_PolyMAC_P0::readOn(Entree& is)
{
  return is;
}

void Production_energie_cin_turb_PolyMAC_P0::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
// empty : no derivative for the turbulent kinetic energy production to add in the blocks
}

void Production_energie_cin_turb_PolyMAC_P0::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Zone_PolyMAC_P0&                      zone = ref_cast(Zone_PolyMAC_P0, equation().zone_dis().valeur());
  const Probleme_base&                       pb = ref_cast(Probleme_base, equation().probleme());
  const Navier_Stokes_std&               eq_qdm = ref_cast(Navier_Stokes_std, pb.equation(0));
  const DoubleTab&                     tab_grad = pb.get_champ("gradient_vitesse").passe();
  const Op_Diff_Turbulent_PolyMAC_P0_Face& Op_diff = ref_cast(Op_Diff_Turbulent_PolyMAC_P0_Face, eq_qdm.operateur(0).l_op_base());
  const Viscosite_turbulente_base&    visc_turb = ref_cast(Viscosite_turbulente_base, Op_diff.correlation().valeur());
  const DoubleTab&                      tab_rho = equation().probleme().get_champ("masse_volumique").passe();
  const DoubleTab&                      tab_alp = equation().probleme().get_champ("alpha").passe();
  const DoubleVect& pe = zone.porosite_elem(), &ve = zone.volumes();

  int Nph = pb.get_champ("vitesse").valeurs().dimension(1), nb_elem = zone.nb_elem(), D = dimension, nf_tot = zone.nb_faces_tot() ;
  int N = equation().inconnue()->valeurs().line_size();

  DoubleTrav Rij(0, Nph, D, D);
  MD_Vector_tools::creer_tableau_distribue(eq_qdm.pression()->valeurs().get_md_vector(), Rij); //Necessary to compare size in reynolds_stress()
  visc_turb.reynolds_stress(Rij);

  for(int e = 0 ; e < nb_elem ; e++)
    for(int n = 0; n<N ; n++)
      {
        double secmem_en = 0;
        for (int d_U = 0; d_U < D; d_U++)
          for (int d_X = 0; d_X < D; d_X++)
            secmem_en += Rij(e, n, d_X, d_U) * tab_grad(nf_tot + d_X + e * D , D * n + d_U) ;
        secmem_en *= (-1) * pe(e) * ve(e) * tab_alp(e, n) * tab_rho(e, n) ;
        secmem(e, n) += std::max(secmem_en, 0.);
      }
}
