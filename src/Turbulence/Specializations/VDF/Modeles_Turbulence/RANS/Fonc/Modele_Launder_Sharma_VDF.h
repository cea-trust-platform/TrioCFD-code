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
// File:        Modele_Launder_Sharma_VDF.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Modeles_Turbulence/RANS/Fonc
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modele_Launder_Sharma_VDF_included
#define Modele_Launder_Sharma_VDF_included

#include <Modele_Jones_Launder_VDF.h>


#include <TRUSTTabs_forward.h>
class Domaine_Cl_VDF;
class Champ_Face_VDF;

class Modele_Launder_Sharma_VDF : public Modele_Jones_Launder_VDF
{

  Declare_instanciable(Modele_Launder_Sharma_VDF);

public :

  DoubleTab& Calcul_Fmu (DoubleTab&,const Domaine_dis_base&,const Domaine_Cl_dis_base&,const DoubleTab&,const Champ_Don&) const override;

  DoubleTab& Calcul_Fmu_BiK ( DoubleTab&,const Domaine_dis_base&,const Domaine_Cl_dis_base&,const DoubleTab&,const DoubleTab&,const Champ_Don& )const override ;

  void associer(const Domaine_dis_base& , const Domaine_Cl_dis_base& ) override;
  Entree& lire(const Motcle&, Entree&);

protected:

  REF(Domaine_VDF) le_dom_VDF;
  REF(Domaine_Cl_VDF) le_dom_Cl_VDF;
};

#endif

