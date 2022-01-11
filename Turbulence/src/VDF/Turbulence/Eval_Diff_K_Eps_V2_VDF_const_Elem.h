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

#include <Eval_Diff_K_Eps_Bas_Re_VDF_const.h>
#include <Champ_Fonc.h>
#include <Eval_VDF_Elem.h>

class Eval_Diff_K_Eps_V2_VDF_const_Elem : public Eval_Diff_K_Eps_Bas_Re_VDF_const, public Eval_VDF_Elem
{

public:
  static constexpr bool CALC_FLUX_FACES_ECH_EXT_IMP = false, CALC_FLUX_FACES_ECH_GLOB_IMP = false, CALC_FLUX_FACES_PAR = false,
                        CALC_FLUX_FACES_SORTIE_LIB = true, CALC_FLUX_FACES_SYMM = true, CALC_FLUX_FACES_PERIO = false;

  inline Eval_Diff_K_Eps_V2_VDF_const_Elem();

  template <typename Type_Double> inline void flux_face(const DoubleTab&, const int , const Symetrie&, int, Type_Double& ) const;
  template <typename Type_Double> inline void flux_face(const DoubleTab&, const int , const Periodique&, int, Type_Double& ) const;
  template <typename Type_Double> inline void flux_face(const DoubleTab&, const int , const Neumann_sortie_libre&, int, Type_Double& ) const;
  template <typename Type_Double> inline void flux_face(const DoubleTab&, const int , const Dirichlet_entree_fluide&, int, Type_Double& ) const;
  template <typename Type_Double> inline void flux_face(const DoubleTab&, const int , const Dirichlet_paroi_fixe&, int, Type_Double& ) const;
  template <typename Type_Double> inline void flux_face(const DoubleTab&, const int , const Dirichlet_paroi_defilante&, int, Type_Double& ) const;
  template <typename Type_Double> inline void flux_face(const DoubleTab&, const int , const Neumann_paroi_adiabatique&, int, Type_Double& ) const;
  template <typename Type_Double> inline void flux_face(const DoubleTab&, const int , const Neumann_paroi&, int, Type_Double& ) const;
  template <typename Type_Double> inline void flux_face(const DoubleTab&, const int , const int, const int, const Echange_externe_impose&, const int, Type_Double& ) const;
  template <typename Type_Double> inline void flux_face(const DoubleTab&, const int , const Echange_global_impose&, int, Type_Double& ) const;
  template <typename Type_Double> inline void flux_faces_interne(const DoubleTab&, const int ,  Type_Double& ) const;

  template <typename Type_Double> inline void coeffs_face(const int, const int, const Symetrie&, Type_Double& , Type_Double&  ) const;
  template <typename Type_Double> inline void coeffs_face(const int, const int, const Neumann_sortie_libre&, Type_Double& , Type_Double&  ) const;
  template <typename Type_Double> inline void coeffs_face(const int, const int, const Dirichlet_entree_fluide&, Type_Double& , Type_Double&  ) const;
  template <typename Type_Double> inline void coeffs_face(const int, const int, const Dirichlet_paroi_fixe&, Type_Double& , Type_Double&  ) const;
  template <typename Type_Double> inline void coeffs_face(const int, const int, const Dirichlet_paroi_defilante&, Type_Double& , Type_Double&  ) const;
  template <typename Type_Double> inline void coeffs_face(const int, const int, const Neumann_paroi_adiabatique&, Type_Double& , Type_Double&  ) const;
  template <typename Type_Double> inline void coeffs_face(const int, const int, const Neumann_paroi&, Type_Double& , Type_Double&  ) const;
  template <typename Type_Double> inline void coeffs_face(const int, const int, const int, const int, const Echange_externe_impose&, Type_Double& , Type_Double&  ) const;
  template <typename Type_Double> inline void coeffs_face(const int, const int, const Echange_global_impose&, Type_Double& , Type_Double&  ) const;
  template <typename Type_Double> inline void coeffs_face(const int, const int, const Periodique&, Type_Double& , Type_Double&  ) const;
  template <typename Type_Double> inline void coeffs_faces_interne(int, Type_Double& , Type_Double&  ) const;

  template <typename Type_Double> inline void secmem_face(const int, const Symetrie&, const int, Type_Double& ) const;
  template <typename Type_Double> inline void secmem_face(const int, const Neumann_sortie_libre&, const int, Type_Double& ) const;
  template <typename Type_Double> inline void secmem_face(const int, const Dirichlet_entree_fluide&, const int, Type_Double& ) const;
  template <typename Type_Double> inline void secmem_face(const int, const Dirichlet_paroi_fixe&, const int, Type_Double& ) const;
  template <typename Type_Double> inline void secmem_face(const int, const Dirichlet_paroi_defilante&, const int, Type_Double& ) const;
  template <typename Type_Double> inline void secmem_face(const int, const Neumann_paroi_adiabatique&, const int, Type_Double& ) const;
  template <typename Type_Double> inline void secmem_face(const int, const Neumann_paroi&, const int, Type_Double& ) const;
  template <typename Type_Double> inline void secmem_face(const int, const int, const int, const Echange_externe_impose&, const int, Type_Double& ) const;
  template <typename Type_Double> inline void secmem_face(const int, const Echange_global_impose&, const int, Type_Double& ) const;
  template <typename Type_Double> inline void secmem_face(const int, const Periodique&, const int, Type_Double& ) const;
  template <typename Type_Double> inline void secmem_faces_interne(const int, Type_Double& ) const;

private:
  REF(Champ_Inc) KEps;
};

inline Eval_Diff_K_Eps_V2_VDF_const_Elem::Eval_Diff_K_Eps_V2_VDF_const_Elem() : Eval_Diff_K_Eps_Bas_Re_VDF_const() {}


////////////////////////////////////////////////////////////////
// Fonctions de calcul des flux pour une grandeur vectorielle
///////////////////////////////////////////////////////////////

//// flux_face avec Dirichlet_entree_fluide
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::flux_face(const DoubleTab& inco, int face,
                                                         const Dirichlet_entree_fluide& la_cl,
                                                         int num1,Type_Double& flux) const
{
  // Cerr << " coucou dans Dirichlet_entree_fluide Eval_Diff_K_Eps_V2_VDF_const_Elem::flux_face " << finl;
  int n0 = elem_(face,0);
  int n1 = elem_(face,1);
  //int k;
  double dist = dist_norm_bord(face);
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


//// coeffs_face avec Dirichlet_entree_fluide
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::coeffs_face(int face, int num1,const Dirichlet_entree_fluide& la_cl,
                                                           Type_Double& aii, Type_Double& ajj) const
{
  //int k;
  int i = elem_(face,0);
  int j = elem_(face,1);
  //double dist = dist_norm_bord(face);

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

//// secmem_face avec Dirichlet_entree_fluide
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::secmem_face(int face, const Dirichlet_entree_fluide& la_cl,
                                                           int num1,Type_Double& flux) const
{
  int i = elem_(face,0);
  int j = elem_(face,1);
  //int k;
  double dist = dist_norm_bord(face);

  if (i != -1)
    {

      flux(0) = la_cl.val_imp(face-num1,0)
                *surface(face)*porosite(face)*(db_diffusivite+dv_diffusivite_turbulente(i)/Prdt_K)/dist;
      flux(1) = la_cl.val_imp(face-num1,1)
                *surface(face)*porosite(face)*(db_diffusivite+dv_diffusivite_turbulente(i)/Prdt_Eps)/dist;
      /*   for (k=0; k<flux.size(); k++) */
      /*       flux(k) = la_cl.val_imp(face-num1,k) */
      /*         *surface(face)*porosite(face)*(db_diffusivite+dv_diffusivite_turbulente(i)/Prdt_K)/dist; */
    }
  else // j != -1
    {

      flux(0) = -la_cl.val_imp(face-num1,0)
                *surface(face)*porosite(face)*(db_diffusivite+dv_diffusivite_turbulente(j)/Prdt_K)/dist;
      flux(1) = -la_cl.val_imp(face-num1,1)
                *surface(face)*porosite(face)*(db_diffusivite+dv_diffusivite_turbulente(j)/Prdt_Eps)/dist;
      /* for (k=0; k<flux.size(); k++)  */
      /*       flux(k) = -la_cl.val_imp(face-num1,k) */
      /*         *surface(face)*(porosite(face)*(db_diffusivite+dv_diffusivite_turbulente(j)/Prdt_K))/dist; */
    }
}

//// flux_face avec Dirichlet_paroi_defilante
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::flux_face(const DoubleTab&, int ,
                                                         const Dirichlet_paroi_defilante&,
                                                         int, Type_Double& ) const
{
  ;
}


//// coeffs_face avec Dirichlet_paroi_defilante
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::coeffs_face(int , int,
                                                           const Dirichlet_paroi_defilante&,
                                                           Type_Double&, Type_Double& ) const
{
  ;
}

//// secmem_face avec Dirichlet_paroi_defilante
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::secmem_face(int, const Dirichlet_paroi_defilante&,
                                                           int, Type_Double& ) const
{
  ;
}

//// flux_face avec Dirichlet_paroi_fixe
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::flux_face(const DoubleTab& inco, int face ,
                                                         const Dirichlet_paroi_fixe& la_cl,
                                                         int num1, Type_Double& flux) const
{

  // Cerr << " coucou dans Dirichlet_paroi_fixe Eval_Diff_K_Eps_V2_VDF_const_Elem::flux_face " << finl;

  assert(dv_diffusivite_turbulente.ref_count() >=2);
  assert(diffusivite_turbulente_->valeurs().addr() == dv_diffusivite_turbulente.addr() );
  // Cerr << " DANS EVAL_DIFF POUR V2 " << finl;
  int n0 = elem_(face,0);
  int n1 = elem_(face,1);
  /*   Cerr << "n0 " << n0 << finl; */
  /*   Cerr << "n1 " << n1 << finl; */
  //const Type_Double& porosite_surf = la_zone->porosite_face();
  //const Type_Double& volume_entrelaces = la_zone->volumes_entrelaces();
  double dist = dist_norm_bord(face);
  /*   Cerr << "dist " << dist << finl; */
  /*   Cerr << "db_diffusivite = " << db_diffusivite << finl; */

  double temp;
  double coef = surface(face)*porosite(face)/dist;

  if (n0 != -1)
    {
      // k=0 a la paroi
      //  Cerr << "dv_diffusivite_turbulente(n0)  = " << dv_diffusivite_turbulente(n0)<< finl;
      flux(0) = -inco(n0,0)*coef*(db_diffusivite+dv_diffusivite_turbulente(n0)/Prdt_K);
      //Cerr << " pour elem " << n0 << " le flux de k = " << flux(0) << finl;
      // epsilon = 2*db_diffusivite*k/dist^2
      temp = 2*db_diffusivite*inco(n0,0)/(dist*dist);
      flux(1) = (temp-inco(n0,1))*coef*(db_diffusivite+dv_diffusivite_turbulente(n0)/Prdt_Eps);
      // Cerr << " pour elem " << n0 << " le flux de eps = " << flux(1) << finl;
    }
  else  // n1 != -1
    {
      // k=0 a la paroi
      // Cerr << "dv_diffusivite_turbulente(n1)  = " << dv_diffusivite_turbulente(n1)<< finl;
      flux(0) = inco(n1,0)*coef*(db_diffusivite+dv_diffusivite_turbulente(n1)/Prdt_K);
      // Cerr << " pour elem " << n1 << " le flux de k = " << flux(0) << finl;
      // epsilon = 2*db_diffusivite*k/dist^2
      temp = 2*db_diffusivite*inco(n1,0)/(dist*dist);
      //Cerr << " valeur de epsilon a la paroi " << temp << finl;
      flux(1) = (inco(n1,1)-temp)*coef*(db_diffusivite+dv_diffusivite_turbulente(n1)/Prdt_Eps);
      //Cerr << " pour elem " << n1 << " le flux de eps = " << flux(1) << finl;
    }
}


//// coeffs_face avec Dirichlet_paroi_fixe
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::coeffs_face(int face,int num1, const Dirichlet_paroi_fixe& la_cl,
                                                           Type_Double& aii, Type_Double& ajj) const
{
  assert(dv_diffusivite_turbulente.ref_count() >=2);
  assert(diffusivite_turbulente_->valeurs().addr() == dv_diffusivite_turbulente.addr() );
  int i = elem_(face,0);
  int j = elem_(face,1);
  //const Type_Double& porosite_surf = la_zone->porosite_face();
  //const Type_Double& volume_entrelaces = la_zone->volumes_entrelaces();
  double dist = dist_norm_bord(face);
  double coef = surface(face)*porosite(face)/dist;

  if (i != -1)
    {
      aii(0) = coef*(db_diffusivite+dv_diffusivite_turbulente(i)/Prdt_K);
      aii(1) = coef*(db_diffusivite+dv_diffusivite_turbulente(i)/Prdt_Eps);
      ajj(0) =  ajj(1) = 0;
    }
  else
    {
      // j != -1

      ajj(0) = coef*(db_diffusivite+dv_diffusivite_turbulente(j)/Prdt_K);
      ajj(1) = coef*(db_diffusivite+dv_diffusivite_turbulente(j)/Prdt_Eps);
      aii(0) =  aii(1) = 0;
    }
}

//// secmem_face avec Dirichlet_paroi_fixe
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::secmem_face(int face, const Dirichlet_paroi_fixe& la_cl,
                                                           int num1, Type_Double& flux) const
{
  assert(dv_diffusivite_turbulente.ref_count() >=2);
  assert(diffusivite_turbulente_->valeurs().addr() == dv_diffusivite_turbulente.addr() );

  //Cerr << "dans Eval_Diff_K_Eps_Bas_Re_VDF_const_Elem::secmem_face " << finl;

  int i = elem_(face,0);
  int j = elem_(face,1);
  //const Type_Double& porosite_surf = la_zone->porosite_face();
  //const Type_Double& volume_entrelaces = la_zone->volumes_entrelaces();
  double dist = dist_norm_bord(face);
  double coef = surface(face)*porosite(face)/dist;
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

//// flux_face avec Echange_externe_impose
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::flux_face(const DoubleTab& inco, int boundary_index, int face, int local_face,
                                                         const Echange_externe_impose& la_cl,
                                                         int num1,Type_Double& flux) const
{
}

//// coeffs_face avec Echange_externe_impose
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::coeffs_face(int boundary_index, int face, int local_face, int num1,
                                                           const Echange_externe_impose& la_cl,
                                                           Type_Double& aii, Type_Double& ajj) const
{
}

//// secmem
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::secmem_face(int boundary_index, int face, int local_face, const Echange_externe_impose& la_cl,
                                                           int num1,Type_Double& flux) const
{
}

//// flux_face avec Echange_global_impose
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::flux_face(const DoubleTab& inco, int face,
                                                         const Echange_global_impose& la_cl,
                                                         int num1,Type_Double& flux) const
{
}

//// coeffs_face avec Echange_global_impose
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::coeffs_face(int face, int num1,
                                                           const Echange_global_impose& la_cl,
                                                           Type_Double& aii, Type_Double& ajj ) const
{
}

//// secmem_face avec Echange_global_impose
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::secmem_face(int face, const Echange_global_impose& la_cl,
                                                           int num1,Type_Double& flux) const
{
}

//// flux_face avec Neumann_paroi
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::flux_face(const DoubleTab& , int face,
                                                         const Neumann_paroi& la_cl,
                                                         int num1,Type_Double& flux) const
{
  ;
}

//// coeffs_face avec Neumann_paroi
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::coeffs_face(int , int,
                                                           const Neumann_paroi& ,
                                                           Type_Double& , Type_Double& ) const
{
  ;
}

//// secmem_face avec Neumann_paroi
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::secmem_face(int face, const Neumann_paroi& la_cl,
                                                           int num1, Type_Double& flux) const
{
  ;
}


//// flux_face avec Neumann_paroi_adiabatique
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::flux_face(const DoubleTab&, int ,
                                                         const Neumann_paroi_adiabatique&,
                                                         int, Type_Double& ) const
{
  ;
}


//// coeffs_face avec Neumann_paroi_adiabatique
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::coeffs_face(int , int,
                                                           const Neumann_paroi_adiabatique&,
                                                           Type_Double&, Type_Double& ) const
{
  ;
}

//// secmem_face avec Neumann_paroi_adiabatique
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::secmem_face(int, const Neumann_paroi_adiabatique&,
                                                           int, Type_Double& ) const
{
  ;
}

//// flux_face avec Neumann_sortie_libre
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::flux_face(const DoubleTab& inco, int face,
                                                         const Neumann_sortie_libre& la_cl,
                                                         int num1, Type_Double& flux ) const
{
  // Cerr << "coucou dans Neumann_sortie_libre Eval_Diff_K_Eps_V2_VDF_const_Elem::flux_face" << finl;
  flux = 0 ;
}

//// coeffs_face avec Neumann_sortie_libre
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::coeffs_face(int face, int num1,const Neumann_sortie_libre& la_cl,
                                                           Type_Double& aii, Type_Double& ajj) const
{
  aii=ajj=0;
}

//// secmem_face avec Neumann_sortie_libre
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::secmem_face(int face, const Neumann_sortie_libre& la_cl,
                                                           int num1 , Type_Double& flux) const
{
  flux = 0 ;
}

//// flux_face avec Symetrie
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::flux_face(const DoubleTab& inco, int face,
                                                         const Symetrie& la_cl,
                                                         int num1, Type_Double& flux) const
{
  //  Cerr << "coucou dans symetrie Eval_Diff_K_Eps_V2_VDF_const_Elem::flux_face" << finl;
  flux = 0;
}

//// coeffs_face avec Symetrie
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::coeffs_face(int face,int num1, const Symetrie& la_cl,
                                                           Type_Double& aii, Type_Double& ajj) const
{
  aii=ajj=0;
}

//// secmem_face avec Symetrie
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::secmem_face(int face, const Symetrie& la_cl,
                                                           int num1, Type_Double& flux) const
{
  flux = 0;
}

//// flux_face avec Periodique
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::flux_face(const DoubleTab& inco, int face,
                                                         const Periodique& la_cl,
                                                         int, Type_Double& flux) const
{

}

//// coeffs_face avec Periodique
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::coeffs_face(int face,int, const Periodique& la_cl,
                                                           Type_Double& aii, Type_Double& ajj ) const
{
}


//// secmem_face avec Periodique
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::secmem_face(int face, const Periodique& la_cl,
                                                           int, Type_Double& flux) const
{
  ;
}

//// flux_faces_interne
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::flux_faces_interne(const DoubleTab& inco,
                                                                  int face,Type_Double& flux) const
{

  // Cerr << "Eval_Diff_K_Eps_V2_VDF_const_Elem::flux_faces_interne" << finl;
  assert(dv_diffusivite_turbulente.ref_count() >=2);
  assert(diffusivite_turbulente_->valeurs().addr() == dv_diffusivite_turbulente.addr() );
  int n0 = elem_(face,0);
  int n1 = elem_(face,1);
  //const Type_Double& porosite_surf = la_zone->porosite_face();
  //const Type_Double& volume_entrelaces = la_zone->volumes_entrelaces();
  //const Type_Double& volumes = la_zone->volumes();
  //const Type_Double& porosite_vol = la_zone->porosite_elem();
  double dist = la_zone->dist_norm(face);
  double coef = surface(face)*porosite(face)/dist;
  double diffu = 0.5*(dv_diffusivite_turbulente(n0)+dv_diffusivite_turbulente(n1));
  flux(0) = (db_diffusivite+diffu/Prdt_K)*coef*(inco(n1,0) - inco(n0,0));
  //  flux(0) = (inco(n1,0) - inco(n0,0))*coef;
  //        Cerr << "1.visco = " << flux(0) << " " << volume_entrelaces(face)*porosite_surf(face)/dist/dist << " " << surface(face)*porosite(face)/dist << " " << volumes(n0)*porosite_vol(n0) << finl;
  //        Cerr << dv_diffusivite_turbulente(n1) << " " << coef << " " << dist << " " << inco(n1,0) << " " <<  inco(n0,0) << " " << n0 << " " << n1 << finl;
  flux(1) = (db_diffusivite+diffu/Prdt_Eps)*coef*(inco(n1,1) - inco(n0,1));
  //Cerr << "1.visco = " << db_diffusivite << finl;
  //Cerr << "1.visco turb = " <<  dv_diffusivite_turbulente(n0) << finl;

}


//// coeffs_faces_interne
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::coeffs_faces_interne(int face, Type_Double& aii, Type_Double& ajj ) const
{
  assert(dv_diffusivite_turbulente.ref_count() >=2);
  assert(diffusivite_turbulente_->valeurs().addr() == dv_diffusivite_turbulente.addr() );
  // Cerr << "Eval_Diff_K_Eps_V2_VDF_const_Elem::coeffs_faces_interne" << finl;
  int i = elem_(face,0);
  int j = elem_(face,1);
  //const Type_Double& porosite_surf = la_zone->porosite_face();
  //const Type_Double& volume_entrelaces = la_zone->volumes_entrelaces();
  double dist = la_zone->dist_norm(face);
  double coef = surface(face)*porosite(face)/dist;
  //  double coef = volume_entrelaces(face)*porosite_surf(face)/dist/dist;
  double diffu = 0.5*(dv_diffusivite_turbulente(i)+dv_diffusivite_turbulente(j));
  aii(0) = ajj(0) = (db_diffusivite+diffu/Prdt_K)*coef;
  aii(1) = ajj(1) = (db_diffusivite+diffu/Prdt_Eps)*coef;
  //Cerr << " aii(0) = ajj(0) = " << aii(0) << " aii(1) = ajj(1) = " << aii(1) << finl;
  //Cerr << "2.visco = " << db_diffusivite << finl;
  //        Cerr << "2.visco turb = " <<  dv_diffusivite_turbulente(i) << " " <<   dv_diffusivite_turbulente(j) << finl;

}


//// secmem_faces_interne
//
template <typename Type_Double>
inline void Eval_Diff_K_Eps_V2_VDF_const_Elem::secmem_faces_interne( int face, Type_Double& flux ) const
{
  assert(dv_diffusivite_turbulente.ref_count() >=2);
  assert(diffusivite_turbulente_->valeurs().addr() == dv_diffusivite_turbulente.addr() );
  flux(0) = 0;
  flux(1) = 0;
}

#endif
