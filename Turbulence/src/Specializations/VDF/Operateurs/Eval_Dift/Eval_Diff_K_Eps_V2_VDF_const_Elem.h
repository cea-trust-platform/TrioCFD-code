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
// File:        Eval_Diff_K_Eps_V2_VDF_const_Elem.h
// Directory:   $TRUST_ROOT/src/VDF/Turbulence
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Eval_Diff_K_Eps_V2_VDF_const_Elem_included
#define Eval_Diff_K_Eps_V2_VDF_const_Elem_included

#include <Eval_Diff_K_Eps_VDF_const.h>
#include <CL_Types_include.h>
#include <Evaluateur_VDF.h>
#include <Eval_VDF_Elem.h>
#include <Ref_Champ_Inc.h>

class Eval_Diff_K_Eps_V2_VDF_const_Elem : public Eval_Diff_K_Eps_Bas_Re_VDF_const, public Eval_VDF_Elem, public Evaluateur_VDF
{

public:
  static constexpr bool CALC_FLUX_FACES_ECH_EXT_IMP = false, CALC_FLUX_FACES_ECH_GLOB_IMP = false, CALC_FLUX_FACES_PAR = false,
                        CALC_FLUX_FACES_SORTIE_LIB = true, CALC_FLUX_FACES_SYMM = true, CALC_FLUX_FACES_PERIO = false;

  template <typename Type_Double, typename BC> inline void flux_face(const DoubleTab&, const int , const BC&, int, Type_Double& ) const { /* do nothing */}
  template <typename Type_Double> inline void flux_face(const DoubleTab&, const int , const Dirichlet_entree_fluide&, int, Type_Double& ) const;
  template <typename Type_Double> inline void flux_face(const DoubleTab&, const int , const Dirichlet_paroi_fixe&, int, Type_Double& ) const;
  template <typename Type_Double> inline void flux_face(const DoubleTab&, const int , const int, const int, const Echange_externe_impose&, const int, Type_Double& ) const { /* do nothing */}
  template <typename Type_Double> inline void flux_faces_interne(const DoubleTab&, const int ,  Type_Double& ) const;

  template <typename Type_Double, typename BC> inline void coeffs_face(const int, const int, const BC&, Type_Double& , Type_Double&  ) const { /* do nothing */}
  template <typename Type_Double> inline void coeffs_face(const int, const int, const Dirichlet_entree_fluide&, Type_Double& , Type_Double&  ) const;
  template <typename Type_Double> inline void coeffs_face(const int, const int, const Dirichlet_paroi_fixe&, Type_Double& , Type_Double&  ) const;
  template <typename Type_Double> inline void coeffs_face(const int, const int, const int, const int, const Echange_externe_impose&, Type_Double& , Type_Double&  ) const { /* do nothing */}
  template <typename Type_Double> inline void coeffs_faces_interne(int, Type_Double& , Type_Double&  ) const;

  template <typename Type_Double, typename BC> inline void secmem_face(const int, const BC&, const int, Type_Double& ) const { /* do nothing */}
  template <typename Type_Double> inline void secmem_face(const int, const Dirichlet_entree_fluide&, const int, Type_Double& ) const;
  template <typename Type_Double> inline void secmem_face(const int, const Dirichlet_paroi_fixe&, const int, Type_Double& ) const;
  template <typename Type_Double> inline void secmem_face(const int, const int, const int, const Echange_externe_impose&, const int, Type_Double& ) const { /* do nothing */}
  template <typename Type_Double> inline void secmem_faces_interne(const int, Type_Double& ) const;

private:
  REF(Champ_Inc) KEps;
};

template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::flux_face(const DoubleTab& inco, int face, const Dirichlet_entree_fluide& la_cl, int num1,Type_Double& flux) const
{
  const int n0 = elem_(face,0), n1 = elem_(face,1);
  const double dist = dist_norm_bord(face);
  if (n0 != -1)
    {
      flux(0) = (la_cl.val_imp(face-num1,0)-inco(n0,0))
                *surface(face)*porosite(face)*(db_diffusivite+dv_diffusivite_turbulente(n0)/Prdt_K)/dist;
      flux(1) = (la_cl.val_imp(face-num1,1)-inco(n0,1))
                *surface(face)*porosite(face)*(db_diffusivite+dv_diffusivite_turbulente(n0)/Prdt_Eps)/dist;
    }
  else  // n1 != -1
    {
      flux(0) = (inco(n1,0)-la_cl.val_imp(face-num1,0))
                *surface(face)*porosite(face)*(db_diffusivite+dv_diffusivite_turbulente(n1)/Prdt_K)/dist;
      flux(1) = (inco(n1,1)-la_cl.val_imp(face-num1,1))
                *surface(face)*porosite(face)*(db_diffusivite+dv_diffusivite_turbulente(n1)/Prdt_Eps)/dist;
    }
}

template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::coeffs_face(int face, int num1,const Dirichlet_entree_fluide& la_cl, Type_Double& aii, Type_Double& ajj) const
{
  const int i = elem_(face,0), j = elem_(face,1);
  if (i != -1)
    {
      aii(0) = surface(face)*porosite(face)*(db_diffusivite+dv_diffusivite_turbulente(i)/Prdt_K);
      aii(1) = surface(face)*porosite(face)*(db_diffusivite+dv_diffusivite_turbulente(i)/Prdt_Eps);
      ajj(0) =  ajj(1) = 0;
    }
  else  // j != -1
    {

      ajj(0) = surface(face)*porosite(face)*(db_diffusivite+dv_diffusivite_turbulente(j)/Prdt_K);
      ajj(1) = surface(face)*porosite(face)*(db_diffusivite+dv_diffusivite_turbulente(j)/Prdt_Eps);
      aii(0) =  aii(1) = 0;
    }
}

template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::secmem_face(int face, const Dirichlet_entree_fluide& la_cl, int num1,Type_Double& flux) const
{
  const int i = elem_(face,0), j = elem_(face,1);
  const double dist = dist_norm_bord(face);

  if (i != -1)
    {
      flux(0) = la_cl.val_imp(face-num1,0)
                *surface(face)*porosite(face)*(db_diffusivite+dv_diffusivite_turbulente(i)/Prdt_K)/dist;
      flux(1) = la_cl.val_imp(face-num1,1)
                *surface(face)*porosite(face)*(db_diffusivite+dv_diffusivite_turbulente(i)/Prdt_Eps)/dist;
    }
  else // j != -1
    {
      flux(0) = -la_cl.val_imp(face-num1,0)
                *surface(face)*porosite(face)*(db_diffusivite+dv_diffusivite_turbulente(j)/Prdt_K)/dist;
      flux(1) = -la_cl.val_imp(face-num1,1)
                *surface(face)*porosite(face)*(db_diffusivite+dv_diffusivite_turbulente(j)/Prdt_Eps)/dist;
    }
}

template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::flux_face(const DoubleTab& inco, int face , const Dirichlet_paroi_fixe& la_cl, int num1, Type_Double& flux) const
{
  assert(dv_diffusivite_turbulente.ref_count() >=2);
  assert(diffusivite_turbulente_->valeurs().addr() == dv_diffusivite_turbulente.addr() );
  const int n0 = elem_(face,0), n1 = elem_(face,1);
  const double dist = dist_norm_bord(face), coef = surface(face)*porosite(face)/dist;
  double temp;

  if (n0 != -1)
    {
      flux(0) = -inco(n0,0)*coef*(db_diffusivite+dv_diffusivite_turbulente(n0)/Prdt_K);
      temp = 2*db_diffusivite*inco(n0,0)/(dist*dist);
      flux(1) = (temp-inco(n0,1))*coef*(db_diffusivite+dv_diffusivite_turbulente(n0)/Prdt_Eps);
    }
  else  // n1 != -1
    {
      flux(0) = inco(n1,0)*coef*(db_diffusivite+dv_diffusivite_turbulente(n1)/Prdt_K);
      temp = 2*db_diffusivite*inco(n1,0)/(dist*dist);
      flux(1) = (inco(n1,1)-temp)*coef*(db_diffusivite+dv_diffusivite_turbulente(n1)/Prdt_Eps);
    }
}

template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::coeffs_face(int face,int num1, const Dirichlet_paroi_fixe& la_cl, Type_Double& aii, Type_Double& ajj) const
{
  assert(dv_diffusivite_turbulente.ref_count() >=2);
  assert(diffusivite_turbulente_->valeurs().addr() == dv_diffusivite_turbulente.addr() );
  const int i = elem_(face,0), j = elem_(face,1);
  const double dist = dist_norm_bord(face), coef = surface(face)*porosite(face)/dist;

  if (i != -1)
    {
      aii(0) = coef*(db_diffusivite+dv_diffusivite_turbulente(i)/Prdt_K);
      aii(1) = coef*(db_diffusivite+dv_diffusivite_turbulente(i)/Prdt_Eps);
      ajj(0) =  ajj(1) = 0;
    }
  else // j != -1
    {
      ajj(0) = coef*(db_diffusivite+dv_diffusivite_turbulente(j)/Prdt_K);
      ajj(1) = coef*(db_diffusivite+dv_diffusivite_turbulente(j)/Prdt_Eps);
      aii(0) =  aii(1) = 0;
    }
}

template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::secmem_face(int face, const Dirichlet_paroi_fixe& la_cl, int num1, Type_Double& flux) const
{
  assert(dv_diffusivite_turbulente.ref_count() >=2);
  assert(diffusivite_turbulente_->valeurs().addr() == dv_diffusivite_turbulente.addr() );
  const int i = elem_(face,0), j = elem_(face,1);
  const double dist = dist_norm_bord(face), coef = surface(face)*porosite(face)/dist;
  double temp;

  if (i != -1)
    {
      flux(0) = 0.;
      temp = 2*db_diffusivite*KEps->valeurs()(i,0)/(dist*dist);
      flux(1) = temp*coef*(db_diffusivite+dv_diffusivite_turbulente(i)/Prdt_Eps);
    }
  else   // j != -1
    {
      flux(0) = 0.;
      temp = 2*db_diffusivite*KEps->valeurs()(j,0)/(dist*dist);
      flux(1) = -temp*coef*(db_diffusivite+dv_diffusivite_turbulente(j)/Prdt_Eps);
    }
}

template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::flux_faces_interne(const DoubleTab& inco, int face,Type_Double& flux) const
{
  assert(dv_diffusivite_turbulente.ref_count() >=2);
  assert(diffusivite_turbulente_->valeurs().addr() == dv_diffusivite_turbulente.addr() );
  const int n0 = elem_(face,0), n1 = elem_(face,1);
  const double dist = la_zone->dist_norm(face), coef = surface(face)*porosite(face)/dist;
  const double diffu = 0.5*(dv_diffusivite_turbulente(n0)+dv_diffusivite_turbulente(n1));
  flux(0) = (db_diffusivite+diffu/Prdt_K)*coef*(inco(n1,0) - inco(n0,0));
  flux(1) = (db_diffusivite+diffu/Prdt_Eps)*coef*(inco(n1,1) - inco(n0,1));
}

template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::coeffs_faces_interne(int face, Type_Double& aii, Type_Double& ajj ) const
{
  assert(dv_diffusivite_turbulente.ref_count() >=2);
  assert(diffusivite_turbulente_->valeurs().addr() == dv_diffusivite_turbulente.addr() );
  const int i = elem_(face,0), j = elem_(face,1);
  const double dist = la_zone->dist_norm(face), coef = surface(face)*porosite(face)/dist;
  const double diffu = 0.5*(dv_diffusivite_turbulente(i)+dv_diffusivite_turbulente(j));
  aii(0) = ajj(0) = (db_diffusivite+diffu/Prdt_K)*coef;
  aii(1) = ajj(1) = (db_diffusivite+diffu/Prdt_Eps)*coef;
}

template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::secmem_faces_interne( int face, Type_Double& flux ) const
{
  assert(dv_diffusivite_turbulente.ref_count() >=2);
  assert(diffusivite_turbulente_->valeurs().addr() == dv_diffusivite_turbulente.addr() );
  flux(0) = 0;
  flux(1) = 0;
}

#endif /* Eval_Diff_K_Eps_V2_VDF_const_Elem_included */
