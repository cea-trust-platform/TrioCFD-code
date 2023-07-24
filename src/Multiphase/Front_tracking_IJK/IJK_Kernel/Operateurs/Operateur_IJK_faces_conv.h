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
// File      : Deriv_Operateur_IJK_faces_conv.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/IJK_Kernel/Operateurs
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Operateur_IJK_faces_conv_included
#define Operateur_IJK_faces_conv_included

#include <TRUST_Deriv.h>
#include <Operateur_IJK_faces_conv_base.h>
#include <OpConvAmontIJK.h>

class Operateur_IJK_faces_conv : public DERIV( Operateur_IJK_faces_conv_base_double )
{
  Declare_instanciable( Operateur_IJK_faces_conv ) ;

public:
  inline double compute_dtstab_convection_local(IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz);
  inline void compute_set(IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz);
  inline void compute_add(IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz);
  inline void initialize(const IJK_Splitting& splitting);
  inline void calculer(const IJK_Field_double& inputx, const IJK_Field_double& inputy, const IJK_Field_double& inputz,
                       const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                       IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz);
  inline void ajouter(const IJK_Field_double& inputx, const IJK_Field_double& inputy, const IJK_Field_double& inputz,
                      const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                      IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz);
  inline void set_velocity_components(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz);
  inline void set_bc(const Boundary_Conditions& bc);
  inline void set_bc_thermique(const Boundary_Conditions_Thermique& bc_th);
  inline void ajouter_avec_u_div_rhou(const IJK_Field_double& rhovx, const IJK_Field_double& rhovy, const IJK_Field_double& rhovz,
                                      const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                      IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz,
                                      IJK_Field_double& div_rho_u);
  inline void calculer_avec_u_div_rhou(const IJK_Field_double& rhovx, const IJK_Field_double& rhovy, const IJK_Field_double& rhovz,
                                       const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                       IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz,
                                       IJK_Field_double& div_rho_u);

  /*
   * ReadOn
   */
  int lire_motcle_non_standard(const Motcle& mot, Entree& is) override;
  void set_param(Param& param);
  void typer_convection_op(const char * convection_op);
  Entree& typer_convection_op(Entree& is);
  Nom get_convection_op_type( Motcle word );
  /*
   * Getters
   */
  inline Nom get_convection_op_option() { return convection_option_; };
  inline Nom get_convection_op() { return convection_op_; };
  /*
   * Setters
   */
protected:
  Motcles convection_op_words_;
  Motcles convection_op_options_;
  Nom prefix_;
  Nom suffix_;
  Nom convection_op_;
  Nom convection_option_;
  int convection_rank_;
  bool is_cast_;
};

inline double Operateur_IJK_faces_conv::compute_dtstab_convection_local(IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  return valeur().compute_dtstab_convection_local(dvx, dvy, dvz);
}

inline void Operateur_IJK_faces_conv::compute_set(IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  valeur().compute_set(dvx, dvy, dvz);
}

inline void Operateur_IJK_faces_conv::compute_add(IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  valeur().compute_add(dvx, dvy, dvz);
}

inline void Operateur_IJK_faces_conv::initialize(const IJK_Splitting& splitting)
{
  if (!is_cast_)
    typer_convection_op("quick");
  valeur().initialize(splitting);
}

inline void Operateur_IJK_faces_conv::calculer(const IJK_Field_double& inputx, const IJK_Field_double& inputy, const IJK_Field_double& inputz,
                                               const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                               IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  valeur().calculer(inputx, inputy, inputz, vx, vy, vz, dvx, dvy, dvz);
}

inline void Operateur_IJK_faces_conv::ajouter(const IJK_Field_double& inputx, const IJK_Field_double& inputy, const IJK_Field_double& inputz,
                                              const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                              IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  valeur().ajouter(inputx, inputy, inputz, vx, vy, vz, dvx, dvy, dvz);
}

inline void Operateur_IJK_faces_conv::set_velocity_components(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz)
{
  valeur().set_velocity_components(vx, vy, vz);
}

inline void Operateur_IJK_faces_conv::set_bc(const Boundary_Conditions& bc)
{
  valeur().set_bc(bc);
}

inline void Operateur_IJK_faces_conv::set_bc_thermique(const Boundary_Conditions_Thermique& bc_th)
{
  valeur().set_bc_thermique(bc_th);
}

inline void Operateur_IJK_faces_conv::ajouter_avec_u_div_rhou(const IJK_Field_double& rhovx, const IJK_Field_double& rhovy, const IJK_Field_double& rhovz,
                                                              const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                              IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz,
                                                              IJK_Field_double& div_rho_u)
{
  valeur().ajouter_avec_u_div_rhou(rhovx, rhovy, rhovz, vx, vy, vz, dvx, dvy, dvz, div_rho_u);
}

inline void Operateur_IJK_faces_conv::calculer_avec_u_div_rhou(const IJK_Field_double& rhovx, const IJK_Field_double& rhovy, const IJK_Field_double& rhovz,
                                                               const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                               IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz,
                                                               IJK_Field_double& div_rho_u)
{
  valeur().calculer_avec_u_div_rhou(rhovx, rhovy, rhovz, vx, vy, vz, dvx, dvy, dvz, div_rho_u);
}

#endif /* Operateur_IJK_faces_conv_included */
