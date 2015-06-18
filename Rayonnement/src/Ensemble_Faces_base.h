/****************************************************************************
* Copyright (c) 2015, CEA
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
// File:        Ensemble_Faces_base.h
// Directory:   $TRUST_ROOT/src/Rayonnement
// Version:     /main/12
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Ensemble_Faces_base_included
#define Ensemble_Faces_base_included

#include <Cond_Lim_Rayo.h>
#include <Ref_Cond_lim_base.h>
#include <Motcle.h>

class Ensemble_Faces_base : public Objet_U
{

  Declare_instanciable(Ensemble_Faces_base);

public:

  inline const Cond_Lim_Rayo& cond_lim_rayo() const;
  inline Cond_Lim_Rayo& cond_lim_rayo();
  inline const Cond_lim_base& la_cl_base() const;
  inline Cond_lim_base& la_cl_base();
  void associer_les_cl(Cond_lim_base& ) ;
  double surface(int ) const;
  double teta_i(int ) ;
  int contient(int) const;
  inline  const IntVect& Table_faces () const ;
  void lire(const Nom&,const Nom&, const Domaine& );
  inline int nb_faces_bord() const;
  int is_ok() const;
protected:

  IntVect num_face_Ensemble;//contient_;
  REF(Cond_lim_base) les_cl_base;
  Cond_Lim_Rayo* la_cond_lim_rayo_;
  int nb_faces_bord_;
  DoubleTab positions_;
};

inline int Ensemble_Faces_base::nb_faces_bord() const
{
  return nb_faces_bord_;
}
// les verifs de type cl seront faites dans associer_cl (?)
inline const Cond_lim_base& Ensemble_Faces_base::la_cl_base() const
{
  return les_cl_base;
}
inline  Cond_lim_base& Ensemble_Faces_base::la_cl_base()
{
  return les_cl_base;
}
inline const Cond_Lim_Rayo& Ensemble_Faces_base::cond_lim_rayo() const
{
  //return (const Cond_Lim_Rayo&) les_cl_base.valeur();
  assert(la_cond_lim_rayo_!=0);
  return *la_cond_lim_rayo_;
}
inline Cond_Lim_Rayo& Ensemble_Faces_base::cond_lim_rayo()
{
  //return (Cond_Lim_Rayo&) les_cl_base.valeur();
  assert(la_cond_lim_rayo_!=0);
  return *la_cond_lim_rayo_;

}
inline  const IntVect& Ensemble_Faces_base::Table_faces () const
{
  return num_face_Ensemble;
}
#endif
