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
// File:        Source_Correction_Deficitaire.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Algos
// Version:     /main/11
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Source_Correction_Deficitaire_included
#define Source_Correction_Deficitaire_included

#include <Source_base.h>
#include <Connectivites_IndGros.h>
#include <Restriction_base.h>
#include <Prolongement_base.h>

#include <Domaine_Cl_dis.h>
class Equation_base;
class Probleme_base;

//////////////////////////////////////////////////////////////////////////////
//
// CLASS: Source_Correction_Deficitaire
// Cette classe modelise les termes sources de corrections deficitaires
// pour les problemes definis sur les grilles grossieres.
//
//////////////////////////////////////////////////////////////////////////////

class Source_Correction_Deficitaire : public Source_base
{

  Declare_base(Source_Correction_Deficitaire);

public:

  inline DoubleTab& ajouter(DoubleTab& ) const override ;
  inline DoubleTab& calculer(DoubleTab& ) const override ;
  void completer() override;
  void associer_pb(const Probleme_base& ) override;
  void associer_domaines(const Domaine_dis_base&, const Domaine_Cl_dis& ) override = 0;
  void mettre_a_jour(double temps) override
  {
    ;
  }


  virtual DoubleTab& calculer_residu(Connectivites_IndGros& connect, Restriction_base& rest, Equation_base& eq_fine) = 0;


  virtual DoubleTab& calculer_residu(Connectivites_base& connect, LIST(OWN_PTR(Prolongement_base))& P, Equation_base& eqG, Nom&) = 0;



  inline void setSource(DoubleTab& tab)
  {
    la_correction = tab;
  }
  inline  DoubleTab& getSource()
  {
    return la_correction ;
  }
  inline  DoubleTab& get_residu()
  {
    return le_residu;
  }

protected:

  DoubleTab la_correction;
  DoubleTab le_residu;

  // Iterateur_Source_VDF iter;
};

/*! @brief
 *
 */
inline DoubleTab& Source_Correction_Deficitaire::ajouter(DoubleTab& resu) const
{
  //Cout << " Dans Source_Correction_Deficitaire le_residu = " << le_residu << finl;
  resu += le_residu;
  return resu;
}

inline DoubleTab& Source_Correction_Deficitaire::calculer(DoubleTab& resu) const
{
  //Cout << " Dans Source_Correction_Deficitaire calculer le_residu = " << le_residu << finl;
  resu = 0;
  resu = le_residu;
  return resu;
}




#endif
