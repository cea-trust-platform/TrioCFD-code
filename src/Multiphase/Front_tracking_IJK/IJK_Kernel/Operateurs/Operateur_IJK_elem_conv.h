/****************************************************************************
* Copyright (c) 2023, CEA
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
// File      : Operateur_IJK_elem_conv.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/IJK_Kernel/Operateurs
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Operateur_IJK_elem_conv_included
#define Operateur_IJK_elem_conv_included

#include <TRUST_Deriv.h>
#include <Operateur_IJK_elem_conv_base.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Operateur_IJK_elem_conv
//
// <Description of class Operateur_IJK_elem_conv>
//
/////////////////////////////////////////////////////////////////////////////

class Operateur_IJK_elem_conv : public DERIV(Operateur_IJK_elem_conv_base_double)
{
  Declare_instanciable( Operateur_IJK_elem_conv ) ;
public :
  inline void initialize(const IJK_Splitting& splitting);
  inline void compute_set(IJK_Field_double& dx);
  inline void compute_add(IJK_Field_double& dx);
  inline void calculer(const IJK_Field_double& field,
                       const IJK_Field_double& vx,
                       const IJK_Field_double& vy,
                       const IJK_Field_double& vz,
                       IJK_Field_double& result);
  inline void ajouter(const IJK_Field_double& field,
                      const IJK_Field_double& vx,
                      const IJK_Field_double& vy,
                      const IJK_Field_double& vz,
                      IJK_Field_double& result);
  Nom get_convection_op_type(Motcle word);
  /*
   * ReadOn
   */
  Entree& typer_convection_op(Entree& is);
  void typer_convection_op(const char * convection_op);
  int lire_motcle_non_standard(const Motcle& mot, Entree& is) override;
  void set_param(Param& param);
  /*
   * Getters
   */

  /*
   * Setters
   */
  inline void set_corrige_flux(Corrige_flux_FT& corrige_flux);

protected:
  Motcles convection_op_words_;
  Nom prefix_;
  Nom suffix_;
  int convection_rank_;
  Nom convection_op_;
  Nom convection_op_option_;
};

inline void Operateur_IJK_elem_conv::initialize(const IJK_Splitting& splitting)
{
  valeur().initialize(splitting);
}

inline void Operateur_IJK_elem_conv::compute_set(IJK_Field_double& dx)
{
  valeur().compute_set(dx);
}

inline void Operateur_IJK_elem_conv::compute_add(IJK_Field_double& dx)
{
  valeur().compute_add(dx);
}

inline void Operateur_IJK_elem_conv::set_corrige_flux(Corrige_flux_FT& corrige_flux)
{
  valeur().set_corrige_flux(corrige_flux);
}

inline void Operateur_IJK_elem_conv::calculer(const IJK_Field_double& field,
                                              const IJK_Field_double& vx,
                                              const IJK_Field_double& vy,
                                              const IJK_Field_double& vz,
                                              IJK_Field_double& result)
{
  valeur().calculer(field, vx, vy, vz, result);
}

inline void Operateur_IJK_elem_conv::ajouter(const IJK_Field_double& field,
                                             const IJK_Field_double& vx,
                                             const IJK_Field_double& vy,
                                             const IJK_Field_double& vz,
                                             IJK_Field_double& result)
{
  valeur().ajouter(field, vx, vy, vz, result);
}

#endif /* Operateur_IJK_elem_conv_included */
