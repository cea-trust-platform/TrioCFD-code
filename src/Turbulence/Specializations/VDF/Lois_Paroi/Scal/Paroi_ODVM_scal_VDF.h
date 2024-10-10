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
// File:        Paroi_ODVM_scal_VDF.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Lois_Paroi/Scal
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Paroi_ODVM_scal_VDF_included
#define Paroi_ODVM_scal_VDF_included

#include <Paroi_scal_hyd_base_VDF.h>

#include <TRUST_Vector.h>
#include <Eq_ODVM.h>
#include <TRUSTTabs_forward.h>
#include <Milieu_base.h>
#include <TRUST_Ref.h>



class Convection_Diffusion_std;
class Champ_Fonc_base;

/*! @brief CLASS: Paroi_ODVM_scal_VDF
 *
 * .SECTION  voir aussi
 *  Turbulence_paroi_base
 *
 */
class Paroi_ODVM_scal_VDF : public Paroi_scal_hyd_base_VDF
{
  Declare_instanciable_sans_constructeur(Paroi_ODVM_scal_VDF);

public:
  int calculer_scal(Champ_Fonc_base& ) override;
  int init_lois_paroi() override;
  double get_Tf0(int num_face) const
  {
    return Tf0(num_face);
  };
  const DoubleVect& get_Tf0() const
  {
    return Tf0;
  };

protected:
  int N, CHT, Timp, diff_unif, stats, compt, check;
  double gamma, temps_deb_stats, dt_impr, Tau,
         temps_dern_post_inst, facteur;

  VECT(Eq_ODVM) eq_odvm ;


private:
  DoubleVect Tmean;                // Mean filtered temperature.
  DoubleVect Trms;                // Mean RMS temperature fluctuations.
  DoubleVect Tpm;                // Mean of temperature fluctuations to check.
  DoubleVect Tf0;
  DoubleVect qb_mean;

  DoubleVect Tm_rms;
  DoubleVect Tm_mean;
  DoubleVect Tp_rms;

  DoubleVect tab_u_star;

};


#endif
