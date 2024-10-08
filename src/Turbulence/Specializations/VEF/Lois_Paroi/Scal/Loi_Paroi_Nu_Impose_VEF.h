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
// File:        Loi_Paroi_Nu_Impose_VEF.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Lois_Paroi/Scal
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Loi_Paroi_Nu_Impose_VEF_included
#define Loi_Paroi_Nu_Impose_VEF_included

#include <Domaine_VEF.h>
#include <Paroi_std_scal_hyd_VEF.h>

#include <Parser_U.h>

/*! @brief classe Loi_Paroi_Nu_Impose_VEF
 *
 *
 *
 * .SECTION  voir aussi
 *   Paroi_std_hyd_VEF
 *
 */
class Loi_Paroi_Nu_Impose_VEF : public Paroi_scal_hyd_base_VEF
{

  Declare_instanciable(Loi_Paroi_Nu_Impose_VEF);

public:
  int calculer_scal(Champ_Fonc_base& ) override ;
  void imprimer_nusselt(Sortie&) const override;

protected:

  OWN_PTR(Champ_Don_base) diam_hydr;
  Parser_U nusselt;

};


///////////////////////////////////////////////////////////
//
//  Fonctions inline de la classe Loi_Paroi_Nu_Impose_VEF
//
///////////////////////////////////////////////////////////


#endif

