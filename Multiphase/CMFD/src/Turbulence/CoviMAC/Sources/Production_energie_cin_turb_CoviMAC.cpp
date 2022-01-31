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
// File:        Production_energie_cin_turb_CoviMAC.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/CoviMAC
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Production_energie_cin_turb_CoviMAC.h>

#include <Zone_CoviMAC.h>
#include <Champ_P0_CoviMAC.h>
#include <Matrix_tools.h>
#include <Probleme_base.h>
#include <grad_Champ_Face_CoviMAC.h>
#include <Champ_Uniforme.h>
#include <Flux_interfacial_base.h>
#include <Milieu_composite.h>
#include <Operateur_Diff.h>
#include <Op_Diff_Turbulent_CoviMAC_Face.h>
#include <Navier_Stokes_std.h>
#include <Viscosite_turbulente_base.h>

Implemente_instanciable(Production_energie_cin_turb_CoviMAC,"Production_energie_cin_turb_P0_CoviMAC", Source_base);

Sortie& Production_energie_cin_turb_CoviMAC::printOn(Sortie& os) const
{
  return os;
}

Entree& Production_energie_cin_turb_CoviMAC::readOn(Entree& is)
{
  return is;
}

void Production_energie_cin_turb_CoviMAC::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
// empty : no derivative for the turbulent kinetic energy production to add in the blocks
}

void Production_energie_cin_turb_CoviMAC::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Zone_CoviMAC&                      zone = ref_cast(Zone_CoviMAC, equation().zone_dis().valeur());
  const Probleme_base&                       pb = ref_cast(Probleme_base, equation().probleme());
  const Navier_Stokes_std&               eq_qdm = ref_cast(Navier_Stokes_std, pb.equation(0));
  const grad_Champ_Face_CoviMAC&           grad = ref_cast(grad_Champ_Face_CoviMAC, eq_qdm.get_champ("gradient_vitesse"));
  const DoubleTab&                     tab_grad = grad.valeurs();
  const Op_Diff_Turbulent_CoviMAC_Face& Op_diff = ref_cast(Op_Diff_Turbulent_CoviMAC_Face, eq_qdm.operateur(0).l_op_base());
  const Viscosite_turbulente_base&    visc_turb = ref_cast(Viscosite_turbulente_base, Op_diff.corr.valeur());
  const DoubleTab&                      tab_rho = equation().probleme().get_champ("masse_volumique").passe();
  const DoubleTab&                      tab_alp = equation().probleme().get_champ("alpha").passe();

  int Nph = pb.get_champ("vitesse").valeurs().dimension(1), nf_tot = zone.nb_faces_tot(), ne = zone.nb_elem(), D = dimension ;
  int N = equation().inconnue()->valeurs().line_size();

  DoubleTrav Rij(0, Nph, D, D);
  MD_Vector_tools::creer_tableau_distribue(eq_qdm.pression()->valeurs().get_md_vector(), Rij); //Necessary to compare size in reynolds_stress()
  visc_turb.reynolds_stress(Rij);

  for(int e = 0 ; e < ne ; e++) for(int n = 0; n<N ; n++)
      {
        double secmem_en = 0;
        for (int d_U = 0; d_U < D; d_U++) for (int d_X = 0; d_X < D; d_X++)
            secmem_en += Rij(e, n, d_X, d_U) * tab_grad(nf_tot + d_X + e * D , D * n + d_U) ;
        secmem_en *= (-1) * tab_alp(e, n) * tab_rho(e, n) ;
        secmem(e, n) += secmem_en;
      }
}
