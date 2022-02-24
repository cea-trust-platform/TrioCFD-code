/****************************************************************************
* Copyright (c) 2019, CEA
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
// File:        Source_Transport_K_Eps_Realisable_VDF_Elem.h
// Directory:   $TRUST_ROOT/src/VDF/Turbulence
// Version:     /main/14
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Source_Transport_K_Eps_Realisable_VDF_Elem_included
#define Source_Transport_K_Eps_Realisable_VDF_Elem_included

#include <Source_Transport_Realisable_VDF_Elem_base.h>
#include <Ref_Transport_K_Eps_Realisable.h>

class Source_Transport_K_Eps_Realisable_VDF_Elem : public Source_Transport_Realisable_VDF_Elem_base
{
  Declare_instanciable_sans_constructeur(Source_Transport_K_Eps_Realisable_VDF_Elem);
public :
  Source_Transport_K_Eps_Realisable_VDF_Elem( double cte2 = C21_R__ ) : Source_Transport_Realisable_VDF_Elem_base(cte2) { C2 = cte2; }
  void mettre_a_jour(double temps);
  virtual void associer_pb(const Probleme_base& );

  inline DoubleTab& ajouter(DoubleTab& resu) const { return Source_Transport_Realisable_VDF_Elem_base::ajouter_keps_real(resu); }

protected :
  REF(Transport_K_Eps_Realisable) eqn_keps_Rea;

private:
  const DoubleTab& get_visc_turb() const;
  const Modele_Fonc_Realisable_base& get_modele_fonc() const;
  void calculer_terme_production_real(const Champ_Face&, const DoubleTab& , const DoubleTab& , DoubleTrav&) const;
  void fill_resu_real(const int , const DoubleTab& , const DoubleTrav& , const DoubleTrav& , const DoubleTrav& , double& , DoubleTab& ) const;
  void fill_coeff_matrice(const int , const DoubleTab& , const DoubleVect& , const DoubleVect& , double& , Matrice_Morse& ) const;
};

#endif /* Source_Transport_K_Eps_Realisable_VDF_Elem_included */
