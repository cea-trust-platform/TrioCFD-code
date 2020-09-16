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
// File:        Op_Conv_ALE_VEF.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/ALE/src
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Op_Conv_ALE_VEF_included
#define Op_Conv_ALE_VEF_included

#include <Op_VEF_Face.h>
#include <Op_Conv_ALE.h>
#include <Ref_Zone_VEF.h>
#include <Ref_Zone_Cl_VEF.h>


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
class Op_Conv_ALE_VEF : public Op_Conv_ALE, public Op_VEF_Face
{
  Declare_instanciable(Op_Conv_ALE_VEF);

public :
  virtual void associer (const Zone_dis& , const Zone_Cl_dis& ,const Champ_Inc& );
  double calculer_dt_stab() const ;
  virtual void remplir_fluent_ALEincluded(DoubleVect& ) const;
  virtual void calculateALEMeshVelocityGradientOnFaces( DoubleTab& ) const;
  virtual void calculateALEjacobian(DoubleTab& jacobianALE) const;
  double application_LIMITEUR(double, double, Motcle&) const;
  virtual DoubleTab& ajouterALE(const DoubleTab&, DoubleTab& ) const;
  void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const;
  void modifier_pour_Cl(Matrice_Morse& matrice, DoubleTab& secmem) const;

  inline void dimensionner(Matrice_Morse& ) const;

protected :
  REF(Zone_VEF) la_zone_vef;
  REF(Zone_Cl_VEF) la_zcl_vef;
  mutable ArrOfInt traitement_pres_bord_;
  mutable ArrOfInt est_une_face_de_dirichlet_;
};

inline  void Op_Conv_ALE_VEF::dimensionner(Matrice_Morse& matrice) const
{
  Op_VEF_Face::dimensionner(la_zone_vef.valeur(),la_zcl_vef.valeur(), matrice);
}


#endif
