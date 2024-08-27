/****************************************************************************
* Copyright (c) 2021, CEA
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
// File:        Paroi_std_hyd_EF.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/EF/Lois_Paroi/Hydr
//
//////////////////////////////////////////////////////////////////////////////

// .NOM Paroi_std_hyd_EF :
// .ENTETE TRUST EF Turbulence
// .LIBRAIRIE libhydEFturb
// .FILE Paroi_std_hyd_EF.h
// .FILE Paroi_std_hyd_EF.cpp

#ifndef Paroi_std_hyd_EF_included
#define Paroi_std_hyd_EF_included

#include <Paroi_hyd_base_EF.h>
#include <distances_EF.h>
#include <Paroi_log_QDM.h>
#include <Domaine_dis.h>
#include <Domaine_Cl_dis.h>
class Champ_Fonc_base;

/*! @brief CLASS: Paroi_std_hyd_EF
 *
 * .SECTION  voir aussi
 *  Turbulence_paroi_base
 *
 */
class Paroi_std_hyd_EF : public Paroi_hyd_base_EF, public Paroi_log_QDM
{

  Declare_instanciable_sans_constructeur(Paroi_std_hyd_EF);

public:

  Paroi_std_hyd_EF();
  void set_param(Param& param) override;
  int init_lois_paroi() override;
  int calculer_hyd(DoubleTab& ) override;
  int calculer_hyd_BiK(DoubleTab& , DoubleTab& ) override;
  int calculer_hyd(DoubleTab& , DoubleTab& ) override;

  void imprimer_ustar(Sortie& ) const override;
  double calculer_u_plus(const int ,const double ,const  double erugu );

  virtual int calculer_k_eps(double& , double& , double , double , double , double);


protected:

  DoubleVect uplus_;

  DoubleVect seuil_LP;
  IntVect iterations_LP;

  double u_star_impose_;
  int is_u_star_impose_;
  virtual int init_lois_paroi_hydraulique();
};


/*! @brief cette classe permet de specifier des options a la loi de paroi standard.
 *
 * Elle est reservee aux experts.
 *
 */
class Loi_expert_hydr_EF : public Paroi_std_hyd_EF
{
  Declare_instanciable(Loi_expert_hydr_EF);
};

#endif
