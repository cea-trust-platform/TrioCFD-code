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
// File:        Source_rayo_semi_transp_VEF_P1NC.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src/VEF
// Version:     /main/10
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Source_rayo_semi_transp_VEF_P1NC_included
#define Source_rayo_semi_transp_VEF_P1NC_included

#include <Source_rayo_semi_transp_base.h>

/*! @brief classe Source_rayo_semi_transp_VEF_P1NC
 *
 *  .SECTION
 *
 *
 */
class Source_rayo_semi_transp_VEF_P1NC: public Source_rayo_semi_transp_base
{
  Declare_instanciable(Source_rayo_semi_transp_VEF_P1NC);

public :

  DoubleTab& ajouter(DoubleTab& resu) const override;
  DoubleTab& calculer(DoubleTab& resu) const override;

  //  void recherche_et_associe_le_modele();

protected :
  void associer_domaines(const Domaine_dis_base& ,const Domaine_Cl_dis_base& ) override;
  void associer_pb(const Probleme_base& ) override;

};

#endif
