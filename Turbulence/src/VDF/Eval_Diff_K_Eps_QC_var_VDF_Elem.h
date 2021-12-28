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
// File:        Eval_Diff_K_Eps_QC_var_VDF_Elem.h
// Directory:   $TRUST_ROOT/src/VDF/Turbulence
// Version:     1
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Eval_Diff_K_Eps_QC_var_VDF_Elem_included
#define Eval_Diff_K_Eps_QC_var_VDF_Elem_included

#include <Eval_Diff_K_Eps_VDF.h>
#include <CL_Types_include.h>
#include <Champ_Fonc.h>
#include <Eval_VDF_Elem.h>
#include <Fluide_Quasi_Compressible.h>
#include <Zone_VDF.h>

class Eval_Diff_K_Eps_QC_var_VDF_Elem : public Eval_Diff_K_Eps_VDF, public Eval_VDF_Elem
{
public:
  inline Eval_Diff_K_Eps_QC_var_VDF_Elem();

  inline void flux_face(const DoubleTab&, int , const Symetrie&, int, ArrOfDouble& flux) const;
  inline void flux_face(const DoubleTab&, int , const Periodique&, int, ArrOfDouble& flux) const;
  inline void flux_face(const DoubleTab&, int , const Neumann_sortie_libre&, int, ArrOfDouble& flux) const;
  inline void flux_face(const DoubleTab&, int , const Dirichlet_entree_fluide&, int, ArrOfDouble& flux) const;
  inline void flux_face(const DoubleTab&, int , const Dirichlet_paroi_fixe&, int, ArrOfDouble& flux) const;
  inline void flux_face(const DoubleTab&, int , const Dirichlet_paroi_defilante&, int, ArrOfDouble& flux) const;
  inline void flux_face(const DoubleTab&, int , const Neumann_paroi_adiabatique&, int, ArrOfDouble& flux) const;
  inline void flux_face(const DoubleTab&, int , const Neumann_paroi&, int, ArrOfDouble& flux) const;
  inline void flux_face(const DoubleTab&, int , int, int, const Echange_externe_impose&, int, ArrOfDouble& flux) const;
  inline void flux_face(const DoubleTab&, int , const Echange_global_impose&, int, ArrOfDouble& flux) const;
  inline void flux_faces_interne(const DoubleTab&, int ,  ArrOfDouble& flux) const;

  inline void coeffs_face(int,int, const Symetrie&, ArrOfDouble& aii, ArrOfDouble& ajj ) const;
  inline void coeffs_face(int, int,const Neumann_sortie_libre&, ArrOfDouble& aii, ArrOfDouble& ajj ) const;
  inline void coeffs_face(int,int, const Dirichlet_entree_fluide&, ArrOfDouble& aii, ArrOfDouble& ajj ) const;
  inline void coeffs_face(int,int, const Dirichlet_paroi_fixe&, ArrOfDouble& aii, ArrOfDouble& ajj ) const;
  inline void coeffs_face(int,int, const Dirichlet_paroi_defilante&, ArrOfDouble& aii, ArrOfDouble& ajj ) const;
  inline void coeffs_face(int,int, const Neumann_paroi_adiabatique&, ArrOfDouble& aii, ArrOfDouble& ajj ) const;
  inline void coeffs_face(int,int, const Neumann_paroi&, ArrOfDouble& aii, ArrOfDouble& ajj ) const;
  inline void coeffs_face(int,int,int,int, const Echange_externe_impose&, ArrOfDouble& aii, ArrOfDouble& ajj ) const;
  inline void coeffs_face(int,int, const Echange_global_impose&, ArrOfDouble& aii, ArrOfDouble& ajj ) const;
  inline void coeffs_face(int,int, const Periodique&, ArrOfDouble& aii, ArrOfDouble& ajj ) const;
  inline void coeffs_faces_interne(int, ArrOfDouble& aii, ArrOfDouble& ajj ) const;

  inline void secmem_face(int, const Symetrie&, int, ArrOfDouble& ) const;
  inline void secmem_face(int, const Neumann_sortie_libre&, int, ArrOfDouble& ) const;
  inline void secmem_face(int, const Dirichlet_entree_fluide&, int, ArrOfDouble& ) const;
  inline void secmem_face(int, const Dirichlet_paroi_fixe&, int, ArrOfDouble& ) const;
  inline void secmem_face(int, const Dirichlet_paroi_defilante&, int, ArrOfDouble& ) const;
  inline void secmem_face(int, const Neumann_paroi_adiabatique&, int, ArrOfDouble& ) const;
  inline void secmem_face(int, const Neumann_paroi&, int, ArrOfDouble& ) const;
  inline void secmem_face(int, int, int, const Echange_externe_impose&, int, ArrOfDouble& ) const;
  inline void secmem_face(int, const Echange_global_impose&, int, ArrOfDouble& ) const;
  inline void secmem_face(int, const Periodique&, int, ArrOfDouble& ) const;
  inline void secmem_faces_interne(int, ArrOfDouble& ) const;

  inline void associer(const Champ_Don& );
  inline void mettre_a_jour( );

protected:
  DoubleVect dv_diffusivite;

};

//
// Fonctions inline de la classe Eval_Diff_K_Eps_QC_var_VDF_Elem
//
// Description:
// associe le champ de diffusivite
inline void Eval_Diff_K_Eps_QC_var_VDF_Elem::associer(const Champ_Don& diffu)
{
  diffusivite_ = diffu.valeur();
  dv_diffusivite.ref(diffu.valeurs());

}

// Description:
// mise a jour de DoubleVect diffusivite
inline void  Eval_Diff_K_Eps_QC_var_VDF_Elem::mettre_a_jour( )
{
  (diffusivite_->valeurs().echange_espace_virtuel());
  dv_diffusivite.ref(diffusivite_->valeurs());
  dv_diffusivite_turbulente.ref(diffusivite_turbulente_->valeurs());
  //Cerr<<"dv_diffusivite "<<dv_diffusivite<<finl;
  //Process::exit();
}

inline Eval_Diff_K_Eps_QC_var_VDF_Elem::Eval_Diff_K_Eps_QC_var_VDF_Elem()
  : Eval_Diff_K_Eps_VDF() {}
// DEBUT DES DEFINES
#define CLASSNAME Eval_Diff_K_Eps_QC_var_VDF_Elem
#define nu_1(i,k) (dv_diffusivite(i)+dv_diffusivite_turbulente(i)/Prdt[k])
#define nu_2(i,k) dv_diffusivite(i);Process::exit();
#define f_heq(d0,i,d1,j,k) heq=0.5*(nu_1(i,k) + nu_1(j,k))/(d1+d0);
#undef D_AXI
#undef DEQUIV
#undef MULTD
#define ISQUASI
#include <Vect_corps_base.h>
#undef CLASSNAME
#undef f_heq
#undef nu_1
#undef nu_2
#undef DEQUIV
#undef MULTD
#undef D_AXI
#undef Dist_face_elem0
#undef Dist_face_elem1
#undef Dist_norm_bord_externe
#undef Dist_norm_bord
#undef MODIF_DEQ
#undef ISQUASI
#endif
