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
// File:        Source_Qdm_VDF_Phase_field.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Mutliphase/Phase_field/src/VDF
// Version:     /main/12
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Source_Qdm_VDF_Phase_field_included
#define Source_Qdm_VDF_Phase_field_included

#include <Navier_Stokes_std.h>
#include <Terme_Source_Qdm.h>
#include <Ref_Domaine_VDF.h>
#include <Ref_Domaine_Cl_VDF.h>

/*! @brief class Source_Qdm_VDF_Phase_field
 *
 * @sa Source_base
 */
class Source_Qdm_VDF_Phase_field : public Source_base, public Terme_Source_Qdm
{

  Declare_instanciable(Source_Qdm_VDF_Phase_field);

public:
  void associer_pb(const Probleme_base& ) override;
  DoubleTab& methode_1(DoubleTab& ) const;
  DoubleTab& methode_2(DoubleTab& ) const;
  DoubleTab& methode_3(DoubleTab& ) const;
  DoubleTab& methode_4(DoubleTab& ) const;
  DoubleTab& ajouter(DoubleTab& ) const override;
  DoubleTab& calculer(DoubleTab& ) const override;
  void mettre_a_jour(double ) override;

protected:
  int terme_source;
  int compressible;
  double rho0;
  DoubleTab grad_mutilde_;
  DoubleTab grad_div_alpha_rho_gradC_;
  DoubleTab grad_alpha_gradC_carre_;
  DoubleTab gradC_;
  int boussi_;

  REF(Domaine_VDF) le_dom_VDF;
  REF(Domaine_Cl_VDF) le_dom_Cl_VDF;
  void associer_domaines(const Domaine_dis& ,const Domaine_Cl_dis& ) override;
  REF(Probleme_base) le_probleme2;
  DoubleTab& (Source_Qdm_VDF_Phase_field::*methode)(DoubleTab& ) const;
};

#endif
