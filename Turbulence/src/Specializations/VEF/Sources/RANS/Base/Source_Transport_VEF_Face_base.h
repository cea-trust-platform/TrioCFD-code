/****************************************************************************
* Copyright (c) 2022, CEA
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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Source_Transport_VEF_Face_base.h
// Directory : $TURBULENCE_ROOT/src/Specializations/VEF/Sources/RANS/Base
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Source_Transport_VEF_Face_base_included
#define Source_Transport_VEF_Face_base_included

#include <Calcul_Production_K_VEF.h>
#include <Source_Transport_proto.h>
#include <Ref_Transport_K_Eps.h>
#include <Ref_Zone_Cl_VEF.h>
#include <Ref_Zone_VEF.h>

class Source_Transport_VEF_Face_base : public Source_base, public Calcul_Production_K_VEF, public Source_Transport_proto
{
  Declare_base_sans_constructeur(Source_Transport_VEF_Face_base);
public :
  Source_Transport_VEF_Face_base() { }
  Source_Transport_VEF_Face_base(double cs1, double cs2) : Source_Transport_proto(cs1,cs2) { }

  void associer_pb(const Probleme_base& pb) override;
  void associer_domaines(const Zone_dis& ,const Zone_Cl_dis& ) override;
  DoubleTab& calculer(DoubleTab& ) const override;
  DoubleTab& ajouter(DoubleTab& ) const override = 0; // XXX XXX XXX Elie Saikali : like that !!;

  inline void mettre_a_jour(double temps) override { Calcul_Production_K_VEF::mettre_a_jour(temps); }
  inline void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const override { /* Do Nothing */ }

protected :
  DoubleTab& ajouter_keps(DoubleTab& ) const;
  DoubleTab& ajouter_anisotherme(DoubleTab& ) const;
  DoubleTab& ajouter_concen(DoubleTab& ) const;
  DoubleTab& ajouter_anisotherme_concen(DoubleTab& ) const;

  REF(Zone_VEF) le_dom_VEF;
  REF(Zone_Cl_VEF) le_dom_Cl_VEF;

private:
  // methodes a surcharger sinon throw !!
  virtual const DoubleTab& get_visc_turb() const { return not_implemented<DoubleTab&>(__func__); }
  virtual const DoubleTab& get_cisaillement_paroi() const { return not_implemented<DoubleTab&>(__func__); }
  virtual const DoubleTab& get_K_pour_production() const { return not_implemented<DoubleTab&>(__func__); }
  virtual const Modele_Fonc_Bas_Reynolds& get_modele_fonc_bas_reyn() const { return not_implemented<Modele_Fonc_Bas_Reynolds&>(__func__); }
  virtual void calcul_tabs_bas_reyn(const DoubleTrav& , const DoubleTab& , const DoubleTab& , const Champ_Don& , const Champ_base& , DoubleTab& , DoubleTab& , DoubleTab& , DoubleTab& ) const { return not_implemented<void>(__func__); }
  virtual const Nom get_type_paroi() const { return not_implemented<Nom>(__func__); }
  virtual void calcul_tenseur_reyn(const DoubleTab& , const DoubleTab& , DoubleTab& ) const { return not_implemented<void>(__func__); }
  virtual void fill_resu_bas_rey(const DoubleVect& , const DoubleTrav& , const DoubleTab& , const DoubleTab& , const DoubleTab& , const DoubleTab& , DoubleTab& ) const { return not_implemented<void>(__func__); }
  virtual void fill_resu(const DoubleVect& , const DoubleTrav& , DoubleTab& ) const { return not_implemented<void>(__func__); }
  virtual void fill_resu_anisotherme(const DoubleVect& , const DoubleVect& , DoubleTab& ) const { return not_implemented<void>(__func__); }
  virtual void fill_resu_concen(const DoubleTrav& , const DoubleVect& , DoubleTab& ) const { return not_implemented<void>(__func__); } // XXX on peut faire une methode unique avec fill_resu_anisotherme ...
  virtual void fill_resu_anisotherme_concen(const DoubleTrav& , const DoubleTrav& , const DoubleVect& , DoubleTab& ) const { return not_implemented<void>(__func__); }
};

#endif /* Source_Transport_VEF_Face_base_included */
