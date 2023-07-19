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
// File      : Operateur_IJK_elem_diff.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/IJK_Kernel/Operateurs
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Operateur_IJK_elem_diff_included
#define Operateur_IJK_elem_diff_included

#include <OpDiffTurbIJKScalar.h>
#include <TRUST_Deriv.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Operateur_IJK_elem_diff
//
// <Description of class Operateur_IJK_elem_diff>
//
/////////////////////////////////////////////////////////////////////////////

class Operateur_IJK_elem_diff : public DERIV(OpDiffIJKScalarGeneric_double)
{

  Declare_instanciable( Operateur_IJK_elem_diff ) ;

public :
  inline void typer_diffusion_op(const char * diffusion_op);
  inline void initialize(const IJK_Splitting& splitting);
  inline void compute_set(IJK_Field_double& dx);
  inline void compute_add(IJK_Field_double& dx);
  inline void calculer(const IJK_Field_double& field,
                       IJK_Field_double& result,
                       const IJK_Field_local_double& boundary_flux_kmin,
                       const IJK_Field_local_double& boundary_flux_kmax);
  inline void ajouter(const IJK_Field_double& field,
                      IJK_Field_double& result,
                      const IJK_Field_local_double& boundary_flux_kmin,
                      const IJK_Field_local_double& boundary_flux_kmax);
  inline void set_conductivity_coefficient(const double& uniform_lambda,
                                           const IJK_Field_local_double& lambda,
                                           IJK_Field_local_double& coeff_field_x,
                                           IJK_Field_local_double& coeff_field_y,
                                           IJK_Field_local_double& coeff_field_z);
  inline void set_uniform_lambda(const double& uniform_lambda);
  inline void set_lambda(const IJK_Field_local_double& lambda);
  inline void set_coeff_x_y_z(IJK_Field_local_double& coeff_field_x,
                              IJK_Field_local_double& coeff_field_y,
                              IJK_Field_local_double& coeff_field_z);
  //  inline void typer(const char * type);
protected:
  Motcles diffusion_op_words_;
  Nom prefix_;
  Nom suffix_;
  int diffusion_rank_;
};

inline void Operateur_IJK_elem_diff::set_conductivity_coefficient(const double& uniform_lambda,
                                                                  const IJK_Field_local_double& lambda,
                                                                  IJK_Field_local_double& coeff_field_x,
                                                                  IJK_Field_local_double& coeff_field_y,
                                                                  IJK_Field_local_double& coeff_field_z)
{
  switch(diffusion_rank_)
    {
    case 0 :
      {
        // Standard
        set_lambda(lambda);
        break;
      }
    case 1 :
      {
        // Uniform
        set_uniform_lambda(uniform_lambda);
        break;
      }
    case 2 :
      {
        // Anisotropic
        set_coeff_x_y_z(coeff_field_x, coeff_field_y, coeff_field_z);
        break;
      }
    case 3 :
      {
        // Vectorial
        set_coeff_x_y_z(coeff_field_x, coeff_field_y, coeff_field_z);
        break;
      }
    case 4 :
      {
        // VectorialAnisotropic
        set_coeff_x_y_z(coeff_field_x, coeff_field_y, coeff_field_z);
        break;
      }
    case 5 :
      {
        // StructuralOnly
        set_coeff_x_y_z(coeff_field_x, coeff_field_y, coeff_field_z);
        break;
      }
    default :
      {
        Cerr << "ERROR : The conductivity can not be set properly" << finl;
        abort();
      }
    }
}

inline void Operateur_IJK_elem_diff::typer_diffusion_op(const char * diffusion_op)
{
  Cerr << "Read and Cast Diffusion operators :" << finl;
  Motcle diffusion_key(diffusion_op);
  diffusion_rank_ = diffusion_op_words_.search(diffusion_key);
  Nom type = "";
  type += prefix_;
  switch(diffusion_rank_)
    {
    case 0 :
      {
        type += "";
        break;
      }
    case 1 :
      {
        type += "Uniform";
        break;
      }
    case 2 :
      {
        type += "Anisotropic";
        break;
      }
    case 3 :
      {
        type += "Vectorial";
        break;
      }
    case 4 :
      {
        type += "VectorialAnisotropic";
        break;
      }
    case 5 :
      {
        type += "StructuralOnly";
        break;
      }
    default :
      {
        Cerr << "ERROR : Diffusion operators that are already implemented are:" << finl;
        Cerr << diffusion_op_words_ << finl;
        abort();
      }
    }
  type += suffix_;
  typer(type);
}

inline void Operateur_IJK_elem_diff::initialize(const IJK_Splitting& splitting)
{
  valeur().initialize(splitting);
}

inline void Operateur_IJK_elem_diff::compute_set(IJK_Field_double& dx)
{
  valeur().compute_set(dx);
}

inline void Operateur_IJK_elem_diff::compute_add(IJK_Field_double& dx)
{
  valeur().compute_add(dx);
}

inline void Operateur_IJK_elem_diff::calculer(const IJK_Field_double& field,
                                              IJK_Field_double& result,
                                              const IJK_Field_local_double& boundary_flux_kmin,
                                              const IJK_Field_local_double& boundary_flux_kmax)
{
  return valeur().calculer(field,
                           result,
                           boundary_flux_kmin,
                           boundary_flux_kmax);
}


inline void Operateur_IJK_elem_diff::ajouter(const IJK_Field_double& field,
                                             IJK_Field_double& result,
                                             const IJK_Field_local_double& boundary_flux_kmin,
                                             const IJK_Field_local_double& boundary_flux_kmax)
{
  return valeur().ajouter(field,
                          result,
                          boundary_flux_kmin,
                          boundary_flux_kmax);
}


inline void Operateur_IJK_elem_diff::set_uniform_lambda(const double& uniform_lambda)
{
  return valeur().set_uniform_lambda(uniform_lambda);
}

inline void Operateur_IJK_elem_diff::set_lambda(const IJK_Field_local_double& lambda)
{
  return valeur().set_lambda(lambda);
}

inline void Operateur_IJK_elem_diff::set_coeff_x_y_z(IJK_Field_local_double& coeff_field_x,
                                                     IJK_Field_local_double& coeff_field_y,
                                                     IJK_Field_local_double& coeff_field_z)
{
  return valeur().set_coeff_x_y_z(coeff_field_x, coeff_field_y, coeff_field_z);
}

#endif /* Operateur_IJK_elem_diff_included */
