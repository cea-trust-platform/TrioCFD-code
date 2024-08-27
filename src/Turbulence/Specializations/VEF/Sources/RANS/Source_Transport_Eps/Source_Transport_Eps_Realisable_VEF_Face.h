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
// File:        Source_Transport_Eps_Realisable_VEF_Face.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Sources
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Source_Transport_Eps_Realisable_VEF_Face_included
#define Source_Transport_Eps_Realisable_VEF_Face_included

#include <Source_Transport_Realisable_VEF_Face_base.h>
#include <Transport_K_ou_Eps_Realisable.h>
#include <TRUST_Ref.h>

class Source_Transport_Eps_Realisable_VEF_Face : public Source_Transport_Realisable_VEF_Face_base
{
  Declare_instanciable_sans_constructeur(Source_Transport_Eps_Realisable_VEF_Face);
public :
  Source_Transport_Eps_Realisable_VEF_Face(double cte2 = C21_R__) : Source_Transport_Realisable_VEF_Face_base(cte2) { }

  DoubleTab& ajouter(DoubleTab& ) const override;
  void associer_pb(const Probleme_base& ) override;
  void mettre_a_jour(double temps) override;

protected :
  REF(Transport_K_ou_Eps_Realisable) eqn_k_Rea, eqn_eps_Rea;

private:
  const DoubleTab& get_visc_turb() const override;
  const Modele_Fonc_Realisable_base& get_modele_fonc() const override;
  void calculer_terme_production_real(const DoubleTab&, const DoubleTab&, DoubleTrav&) const override;
  void fill_resu_real(const int, const DoubleVect&, const DoubleTrav&, const DoubleTab&, const DoubleTab&, const double, DoubleTab&) const override;
  inline int get_size_k_eps() const override { return eqn_k_Rea->inconnue()->valeurs().dimension(0); }
  void fill_coeff_matrice(const int , const DoubleVect& , const DoubleVect& , const double , Matrice_Morse& ) const override;
};

#endif /* Source_Transport_Eps_Realisable_VEF_Face_included */
