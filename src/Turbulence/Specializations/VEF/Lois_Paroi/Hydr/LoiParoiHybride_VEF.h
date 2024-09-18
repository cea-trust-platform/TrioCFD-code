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
// File:        LoiParoiHybride_VEF.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Lois_Paroi/Hydr
//
//////////////////////////////////////////////////////////////////////////////



#ifndef LoiParoiHybride_VEF_included
#define LoiParoiHybride_VEF_included

#include <Paroi_hyd_base_VEF.h>
#include <LoiParoiHybride.h>
#include <TRUST_Ref.h>

class Domaine_Cl_VEF;
class Domaine_VEF;
class Domaine_dis_base;
class Domaine_Cl_dis_base;



/*! @brief CLASS: LoiParoiHybride_VEF
 *
 * .SECTION  voir aussi
 *  Turbulence_paroi_base
 *
 */
class LoiParoiHybride_VEF : public Paroi_hyd_base_VEF, LoiParoiHybride
{

  Declare_instanciable(LoiParoiHybride_VEF);

public:

  void associer(const Domaine_dis_base&, const Domaine_Cl_dis_base& ) override;
  int init_lois_paroi() override;
  int calculer_hyd(DoubleTab& ) override;
  int calculer_hyd(DoubleTab& , DoubleTab& ) override;
  int calculer_hyd_BiK(DoubleTab& , DoubleTab& ) override;

protected:


private:

  REF(Domaine_VEF) le_dom_VEF;
  REF(Domaine_Cl_VEF) le_dom_Cl_VEF;

};

///////////////////////////////////////////////////////////
//
//  Fonctions inline de la classe Paroi_std_hyd_VEF
//
///////////////////////////////////////////////////////////



#endif
