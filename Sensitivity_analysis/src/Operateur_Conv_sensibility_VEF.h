/****************************************************************************
* Copyright (c) 2020, CEA
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
// File      : Operateur_Conv_sensibility_VEF.h
// Directory : $$Sensitivity_analysis/src
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Operateur_Conv_sensibility_VEF_included
#define Operateur_Conv_sensibility_VEF_included

#include <Operateur_Conv_sensibility.h>
#include <Ref_Zone_VEF.h>
#include <Ref_Zone_Cl_VEF.h>

/*! @brief : class Operateur_Conv_sensibility_VEF
 *
 *  <Description of class Operateur_Conv_sensibility_VEF>
 *
 *
 *
 */
class Operateur_Conv_sensibility_VEF : public Operateur_Conv_sensibility
{

  Declare_instanciable( Operateur_Conv_sensibility_VEF ) ;

public :
  void associer (const Zone_dis& , const Zone_Cl_dis& ,const Champ_Inc& ) override;
  DoubleTab& ajouter(const DoubleTab&, DoubleTab& ) const override;
  void ajouter_Lstate_sensibility_Amont(const DoubleTab&, const DoubleTab&, DoubleTab& ) const; //L(U0)U1
  void ajouter_Lsensibility_state_Amont(const DoubleTab&, const DoubleTab&, DoubleTab& ) const;//L(U1)U0
  void calcul_vc(const ArrOfInt&, ArrOfDouble& , const ArrOfDouble& s, const DoubleTab& , const DoubleTab& ,int  ) const;
  virtual void remplir_fluent(DoubleVect& ) const;
  double calculer_dt_stab() const override;
  void  add_diffusion_term(const DoubleTab&, DoubleTab&) const;
  void add_diffusion_scalar_term(const DoubleTab&, DoubleTab&, double diffu=1.) const;
  inline double viscA(int face_i, int face_j, int num_elem, double diffu=1.) const;
  void ajouter_conv_term(const Champ_Inc_base&, const DoubleTab&, DoubleTab&, DoubleTab& ) const;
  double application_LIMITEUR(double, double, Motcle&) const;

protected :
  REF(Zone_VEF) la_zone_vef;
  REF(Zone_Cl_VEF) la_zcl_vef;
  mutable DoubleVect fluent;           // tableau qui sert pour le calcul du pas de temps de stabilite
  mutable ArrOfInt traitement_pres_bord_;
  mutable ArrOfInt est_une_face_de_dirichlet_;

};


#endif /* Operateur_Conv_sensibility_VEF_included */
