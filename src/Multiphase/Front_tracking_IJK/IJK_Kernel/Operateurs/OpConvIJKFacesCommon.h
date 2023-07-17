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

#ifndef OpConvIJKFacesCommon
#define OpConvIJKFacesCommon
#include <IJK_Field.h>
#include <IJK_Field_simd_tools.h>
#include <TRUSTTab.h>
#include <IJK_Splitting.h>
#include <Operateur_IJK_base.h>
#include <Operateur_IJK_data_channel.h>
#include <Boundary_Conditions.h>
#include <Boundary_Conditions_Thermique.h>

class OpConvIJKFacesCommon_double : public Operateur_IJK_faces_base_double
{
  Declare_base(OpConvIJKFacesCommon_double);
public:
  virtual void initialize(const IJK_Splitting& splitting);
  virtual void set_bc(const Boundary_Conditions& bc) { ; };
  virtual void set_bc_thermique(const Boundary_Conditions_Thermique& bc_th) { ; };
//  void calculer(const IJK_Field_double& inputx, const IJK_Field_double& inputy, const IJK_Field_double& inputz,
//                const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
//                IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz);
//  void ajouter(const IJK_Field_double& inputx, const IJK_Field_double& inputy, const IJK_Field_double& inputz,
//               const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
//               IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz);
  void calculer(const IJK_Field_double& inputx, const IJK_Field_double& inputy, const IJK_Field_double& inputz,
                const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz);
  void ajouter(const IJK_Field_double& inputx, const IJK_Field_double& inputy, const IJK_Field_double& inputz,
               const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
               IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz);
  void set_velocity_components(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz) {vx_ = &vx; vy_ = &vy; vz_ = &vz; };
  virtual void calculer_avec_u_div_rhou(const IJK_Field_double& rhovx, const IJK_Field_double& rhovy, const IJK_Field_double& rhovz,
                                        const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                        IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz,
                                        IJK_Field_double& div_rho_u) { ; };
  virtual void ajouter_avec_u_div_rhou(const IJK_Field_double& rhovx, const IJK_Field_double& rhovy, const IJK_Field_double& rhovz,
                                       const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                       IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz,
                                       IJK_Field_double& div_rho_u) { ; };
protected:

  Operateur_IJK_data_channel channel_data_;

  // Pointers to input velocity field
  const IJK_Field_double *vx_;
  const IJK_Field_double *vy_;
  const IJK_Field_double *vz_;
  bool perio_k_ ;
  // Pointers to input convected field
  const IJK_Field_double *inputx_;
  const IJK_Field_double *inputy_;
  const IJK_Field_double *inputz_;

  inline const IJK_Field_double& get_input(DIRECTION _DIR_)
  {
    switch(_DIR_)
      {
      case DIRECTION::X:
        return *inputx_;
      case DIRECTION::Y:
        return *inputy_;
      case DIRECTION::Z:
        return *inputz_;
      default:
        Cerr << "Error in OpConvIJKFacesCommon::get_input: wrong direction..." << finl;
        Process::exit();
      }
    // for compilation only...
    return *inputx_;
  }

  inline const IJK_Field_double& get_v(DIRECTION _DIR_)
  {
    switch(_DIR_)
      {
      case DIRECTION::X:
        return *vx_;
      case DIRECTION::Y:
        return *vy_;
      case DIRECTION::Z:
        return *vz_;
      default:
        Cerr << "Error in OpConvIJKFacesCommon::get_v: wrong direction..." << finl;
        Process::exit();
      }
    // for compilation only...
    return *vx_;
  }
};

#endif
