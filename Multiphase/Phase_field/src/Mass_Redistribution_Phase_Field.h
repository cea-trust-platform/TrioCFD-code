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

/*
 * Mass_Redistribution_Phase_Field.h
 *
 *  Created on: May 6, 2022
 *      Author: Shambhavi Nandan
 */

#ifndef Mass_Redistribution_Phase_Field_included
#define Mass_Redistribution_Phase_Field_included

#include <Domaine_VDF.h>
#include <Domaine_Cl_VDF.h>
#include <Equation.h>
#include <Process.h>
#include <vector>
#include <MD_Vector_tools.h>
#include <TRUST_Ref.h>

class Mass_Redistribution_Phase_Field
{

public:
  Mass_Redistribution_Phase_Field() = delete;
  ~Mass_Redistribution_Phase_Field() = delete;

  static void impose_mass_redistribution(const Domaine_VDF&, DoubleTab&, DoubleVect, DoubleVect);
  static DoubleTab c_ini;


private:
  REF(Domaine_VDF) le_dom_vdf;
  REF(Domaine_Cl_VDF) la_zcl_vdf;

public:
  static double epsilon_mass_redistribute;
};



#endif /* MASS_REDISTRIBUTION_PHASE_FIELD_H_ */
