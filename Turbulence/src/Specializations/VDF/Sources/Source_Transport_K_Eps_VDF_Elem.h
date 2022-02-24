/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
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
// File:        Source_Transport_K_Eps_VDF_Elem.h
// Directory:   $TRUST_ROOT/src/VDF/Turbulence
// Version:     /main/16
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Source_Transport_K_Eps_VDF_Elem_included
#define Source_Transport_K_Eps_VDF_Elem_included

#include <Source_Transport_VDF_Elem_base.h>
#include <Ref_Transport_K_Eps.h>

//.DESCRIPTION class Source_Transport_K_Eps_VDF_Elem
// Cette classe represente le terme source qui figure dans l'equation de transport du couple (k,eps) dans le cas ou les equations
// de Navier-Stokes ne sont pas couplees a la thermique ou a l'equation de convection-diffusion d'une concentration.
class Source_Transport_K_Eps_VDF_Elem : public Source_Transport_VDF_Elem_base
{
  Declare_instanciable_sans_constructeur(Source_Transport_K_Eps_VDF_Elem);
public:
  Source_Transport_K_Eps_VDF_Elem(double cte1 = C1__, double cte2 = C2__ ) : Source_Transport_VDF_Elem_base(cte1,cte2) // TODO : FIXME : comprends pas ...
  {
    C1 = cte1;
    C2 = cte2;
  }

  void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const ;
  inline DoubleTab& ajouter(DoubleTab& resu) const { return Source_Transport_VDF_Elem_base::ajouter_keps(resu); }

protected:
  REF(Transport_K_Eps)  mon_eq_transport_K_Eps;
  virtual void associer_pb(const Probleme_base& pb);

private:
  const DoubleTab& get_visc_turb() const;
  const Modele_Fonc_Bas_Reynolds& get_modele_fonc_bas_reyn() const ;
  void calculer_terme_production(const Champ_Face&, const DoubleTab& , const DoubleTab& , DoubleVect&) const;
  void calcul_D_E(const DoubleTab& , const DoubleTab& , const Champ_Don& , DoubleTab& , DoubleTab& ) const;
  void calcul_F1_F2(const Champ_base& , DoubleTab& , DoubleTab& , DoubleTab& , DoubleTab& ) const;
  void fill_resu_bas_rey(const DoubleVect& , const DoubleTab& , const DoubleTab& , const DoubleTab& , const DoubleTab& , DoubleTab& ) const;
  void fill_resu(const DoubleVect& , DoubleTab& ) const;
};

#endif /* Source_Transport_K_Eps_VDF_Elem_included */
