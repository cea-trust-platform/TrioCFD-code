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
// File:        Source_DC_VDF_NS.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/VDF
// Version:     /main/12
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Source_DC_VDF_NS_included
#define Source_DC_VDF_NS_included

#include <Source_Correction_Deficitaire.h>
#include <Probleme_base.h>


class Probleme_base;

//////////////////////////////////////////////////////////////////////////////
//
// CLASS: Source_DC_VDF_NS
// Cette classe modelise les termes sources de corrections deficitaires
// pour les problemes en erreur definis sur les grilles fines.
//
//////////////////////////////////////////////////////////////////////////////

class Source_DC_VDF_NS : public Source_Correction_Deficitaire
{

  Declare_instanciable(Source_DC_VDF_NS);

public:

  void associer_domaines(const Domaine_dis_base&, const Domaine_Cl_dis_base& ) override;



  DoubleTab& calculer_residu(Connectivites_IndGros& connect, Restriction_base& rest, Equation_base& eq_fine) override ;

  inline DoubleTab& calculer_residu(Connectivites_base& connect, LIST(OWN_PTR(Prolongement_base))& P, Equation_base& eqG, Nom&) override ;




};

inline DoubleTab& Source_DC_VDF_NS::calculer_residu(Connectivites_IndGros& connect, Restriction_base& rest, Equation_base& eq_fine)
{
  Cerr<<"N'est pas codee avec ces arguments dans la classe Source_DC_VDF_NS !!"<<finl;
  exit();
  throw;
}



#endif
