/****************************************************************************
* Copyright (c) 2023, CEA
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
// File:        Source_Transport_K_Omega_VDF_Elem.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Sources
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Transport_K_Omega_VDF_Elem.h>
#include <Modele_turbulence_hyd_K_Omega.h>
#include <Milieu_base.h>
#include <TRUSTTrav.h>

Implemente_instanciable_sans_constructeur(Source_Transport_K_Omega_VDF_Elem,
                                          "Source_Transport_K_Omega_VDF_P0_VDF",
                                          Source_Transport_K_Omega_VDF_Elem_base);

Sortie& Source_Transport_K_Omega_VDF_Elem::printOn(Sortie& s) const { return s << que_suis_je() ; }

Entree& Source_Transport_K_Omega_VDF_Elem::readOn(Entree& is)
{
  Source_Transport_K_Omega_VDF_Elem_base::verifier_pb_komega(mon_equation->probleme(), que_suis_je());
  return Source_Transport_K_Omega_VDF_Elem_base::readOn(is);
}

void Source_Transport_K_Omega_VDF_Elem::associer_pb(const Probleme_base& pb)
{
  Source_Transport_K_Omega_VDF_Elem_base::associer_pb(pb);
  eqn_K_Omega = ref_cast(Transport_K_Omega, equation());
}

const DoubleTab& Source_Transport_K_Omega_VDF_Elem::get_visc_turb() const
{
  return eqn_K_Omega->modele_turbulence().viscosite_turbulente()->valeurs();
}

void Source_Transport_K_Omega_VDF_Elem::calculer_terme_production(const Champ_Face_VDF& vitesse, const DoubleTab& visco_turb, const DoubleTab& vit, DoubleVect& P) const
{
  const DoubleTab& K_Omega = eqn_K_Omega->inconnue().valeurs();
  if (axi) calculer_terme_production_K_Axi(le_domaine_VDF.valeur(), vitesse, P, K_Omega, visco_turb);
  else calculer_terme_production_K_for_komega(le_domaine_VDF.valeur(), le_domaine_Cl_VDF.valeur(), P, K_Omega, vit, vitesse, visco_turb);
}

void Source_Transport_K_Omega_VDF_Elem::fill_resu(const DoubleVect& P, DoubleTab& resu) const
{
  const DoubleVect& volumes = le_domaine_VDF->volumes();
  const DoubleVect& porosite_vol = le_domaine_Cl_VDF->equation().milieu().porosite_elem();
  const DoubleTab& K_Omega = eqn_K_Omega->inconnue().valeurs();
  // const double K_MIN = eqn_K_Omega->modele_turbulence().get_K_MIN();
  for (int elem = 0; elem < le_domaine_VDF->nb_elem(); elem++)
    {
      // cAlan : enfin. A adapter.
      double volporo = volumes(elem)*porosite_vol(elem);
      resu(elem, 0) += (P(elem) - BETA_K*K_Omega(elem, 0)*K_Omega(elem, 1))*volporo;
      resu(elem, 1) += (ALPHA_OMEGA*P(elem)*K_Omega(elem, 1)/K_Omega(elem, 0)
                        - BETA_OMEGA*K_Omega(elem, 1)*K_Omega(elem, 1));
      // if (K_eps(elem,0) >= LeK_MIN)
      //   resu(elem,1) += (C1*P(elem)- C2*K_eps(elem,1))*volumes(elem)*porosite_vol(elem)*K_eps(elem,1)/K_eps(elem,0);
    }
}

void Source_Transport_K_Omega_VDF_Elem::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  Source_Transport_K_Omega_VDF_Elem_base::ajouter_komega(secmem);

  const std::string& nom_inco = equation().inconnue().le_nom().getString();
  Matrice_Morse* mat = matrices.count(nom_inco) ? matrices.at(nom_inco) : nullptr;
  if(!mat) return;

  const DoubleTab& val = equation().inconnue().valeurs();
  const DoubleVect& porosite = le_domaine_Cl_VDF->equation().milieu().porosite_elem(), &volumes = le_domaine_VDF->volumes();
  const int size = val.dimension(0);
  // on implicite le -eps et le -eps^2/k
  // cAlan : impliciter omega ?

  for (int c = 0; c < size; c++)
    {
      // -eps*vol  donne +vol dans la bonne case
      if (val(c, 0) > DMINFLOAT)
        {
          // cAlan : a adapter
          double coef_k = porosite(c)*volumes(c)*val(c, 1)/val(c, 0);
          (*mat)(c*2, c*2) += coef_k;
          double coef_omega = ALPHA_OMEGA*coef_k;
          // if (is_modele_fonc) coef_eps*=F2(c);
          (*mat)(c*2+1,c*2+1) += coef_omega;
        }
    }


}
