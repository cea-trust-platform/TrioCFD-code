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

#include <Op_Diff_Turbulent_PolyMAC_P0_Face.h>
#include <Echelle_temporelle_turbulente.h>
#include <Taux_dissipation_turbulent.h>
#include <Viscosite_turbulente_base.h>
#include <Domaine_PolyMAC_P0.h>
#include <Navier_Stokes_std.h>
#include <Pb_Multiphase.h>
#include <Synonyme_info.h>

Implemente_instanciable(Production_energie_cin_turb_PolyMAC_P0,"Production_energie_cin_turb_Elem_PolyMAC_P0", Source_Production_energie_cin_turb);
Add_synonym(Production_energie_cin_turb_PolyMAC_P0, "Production_energie_cin_turb_Elem_PolyMAC_P0P1NC");
Add_synonym(Production_energie_cin_turb_PolyMAC_P0, "Production_energie_cin_turb_Elem_PolyVEF_P0");
Add_synonym(Production_energie_cin_turb_PolyMAC_P0, "Production_energie_cin_turb_Elem_PolyVEF_P0P1");
Add_synonym(Production_energie_cin_turb_PolyMAC_P0, "Production_energie_cin_turb_Elem_PolyVEF_P0P1NC");

Sortie& Production_energie_cin_turb_PolyMAC_P0::printOn(Sortie& os) const {return Source_Production_energie_cin_turb::printOn(os);}
Entree& Production_energie_cin_turb_PolyMAC_P0::readOn(Entree& is) { return Source_Production_energie_cin_turb::readOn(is);}

void Production_energie_cin_turb_PolyMAC_P0::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Domaine_PolyMAC_P0&             domaine = ref_cast(Domaine_PolyMAC_P0, equation().domaine_dis().valeur());
  const Probleme_base&                       pb = ref_cast(Probleme_base, equation().probleme());
  const Navier_Stokes_std&               eq_qdm = ref_cast(Navier_Stokes_std, pb.equation(0));
  const Viscosite_turbulente_base&    visc_turb = ref_cast(Viscosite_turbulente_base, ref_cast(Op_Diff_Turbulent_PolyMAC_P0_Face, eq_qdm.operateur(0).l_op_base()).correlation().valeur());
  const DoubleVect& pe = equation().milieu().porosite_elem(), &ve = domaine.volumes();

  std::string Type_diss = ""; // omega or tau dissipation
  for (int i = 0 ; i < equation().probleme().nombre_d_equations() ; i++)
    {
      if      sub_type(Echelle_temporelle_turbulente, equation().probleme().equation(i)) Type_diss = "tau";
      else if sub_type(Taux_dissipation_turbulent, equation().probleme().equation(i)) Type_diss = "omega";
    }

  int Nph = pb.get_champ("vitesse").valeurs().dimension(1), nb_elem = domaine.nb_elem(), D = dimension, nf_tot = domaine.nb_faces_tot() ;
  int N = equation().inconnue()->valeurs().line_size();
  int e, n;

  const DoubleTab& tab_rho = equation().probleme().get_champ("masse_volumique").passe(),
                   &tab_grad = pb.get_champ("gradient_vitesse").passe(),
                    &k = equation().probleme().get_champ("k").valeurs(),
                     *diss = equation().probleme().has_champ(Type_diss) ? &equation().probleme().get_champ(Type_diss).valeurs() : nullptr,
                      *pdiss = semi_impl.count(Type_diss) ? &semi_impl.at(Type_diss) : (equation().probleme().has_champ(Type_diss) ? &equation().probleme().get_champ(Type_diss).passe() : nullptr );

  if (Type_diss == "")
    {
      DoubleTrav nut(0, Nph);
      MD_Vector_tools::creer_tableau_distribue(eq_qdm.pression()->valeurs().get_md_vector(), nut); //Necessary to compare size in eddy_viscosity()
      visc_turb.eddy_viscosity(nut);

      for( e = 0 ; e < nb_elem ; e++)
        for( n = 0; n<N ; n++)
          {
            double secmem_en = 0.;
            for (int d_U = 0; d_U < D; d_U++)
              for (int d_X = 0; d_X < D; d_X++)
                secmem_en += ( tab_grad(nf_tot + d_U + e * D , D * n + d_X) + tab_grad(nf_tot + d_X + e * D , D * n + d_U) ) * tab_grad(nf_tot + d_X + e * D , D * n + d_U) ;
            secmem_en *= pe(e) * ve(e) * tab_rho(e, n) * nut(e, n) ;

            secmem(e, n) += std::max(secmem_en, 0.);
          }
    }

  else
    {
      for( e = 0 ; e < nb_elem ; e++)
        for(n = 0; n<N ; n++)
          {
            double grad_grad = 0.;
            for (int d_U = 0; d_U < D; d_U++)
              for (int d_X = 0; d_X < D; d_X++)
                grad_grad += ( tab_grad(nf_tot + d_U + e * D , D * n + d_X) + tab_grad(nf_tot + d_X + e * D , D * n + d_U) ) * tab_grad(nf_tot + d_X + e * D , D * n + d_U) ;

            double fac = std::max(grad_grad, 0.) * pe(e) * ve(e) ;

            if (Type_diss == "tau")
              {
                secmem(e, n) += fac * k(e, n) * (*diss)(e, n) ;
                for (auto &&i_m : matrices)
                  {
                    Matrice_Morse& mat = *i_m.second;
                    if (i_m.first == "k")         mat(N * e + n,  N * e + n) -= fac * (*diss)(e,n) + (*pdiss)(e,n)-(*pdiss)(e,n);
                    if (i_m.first == "tau")       mat(N * e + n,  N * e + n) -= fac * k(e,n);
                  }
              }
            else if (Type_diss == "omega")
              {
                secmem(e, n) += fac * k(e, n) / (*pdiss)(e, n)*(2-(*diss)(e, n)/(*pdiss)(e, n)) ;
                for (auto &&i_m : matrices)
                  {
                    Matrice_Morse& mat = *i_m.second;
                    if (i_m.first == "omega")     mat(N * e + n,  N * e + n) += fac * k(e,n) / ((*pdiss)(e, n)*(*pdiss)(e, n));
                    if (i_m.first == "k")         mat(N * e + n,  N * e + n) -= fac * 1./ (*pdiss)(e, n)*(2-(*diss)(e, n)/(*pdiss)(e, n)) ;
                  }
              }
          }
    }
}
