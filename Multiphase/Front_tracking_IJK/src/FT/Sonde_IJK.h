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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Sonde_IJK.h
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Sonde_IJK_included
#define Sonde_IJK_included

#include <Sonde.h>
// #include <IJK_FT.h>
#include <Ref_IJK_FT_double.h>
#include <Ref_IJK_Field_double.h>
#include <Postraitement.h>

/////////////////////////////////////////////////////////////////////////////
//
// .NAME        : Sonde_IJK
// .DESCRIPTION : class Sonde_IJK
//
// <Description of class Sonde_IJK>
//
/////////////////////////////////////////////////////////////////////////////

class Sonde_IJK : public Sonde
{

  Declare_instanciable( Sonde_IJK ) ;

public :
  inline Sonde_IJK(const Nom& );
  void completer(const IJK_FT_double& ijk_ft);
  void initialiser();
  void nommer(const Nom& nom) override
  {
    nom_=nom;
  }
  void mettre_a_jour(const double temps, const double dt);
  void postraiter();
protected :

  REF(IJK_FT_double) ref_ijk_ft_;
  REF(IJK_Field_double) ref_ijk_field_;
  Postraitement post_bidon_;
  Nom nom_;                                // le nom de la sonde
  int ncomp;                           // Numero de la composante a sonder
  SFichier* le_fichier;
  int nbre_points_tot;
};


inline Sonde_IJK::Sonde_IJK(const Nom& nom)
  : nom_(nom), ncomp(-1),le_fichier(0)
{}

#endif /* Sonde_IJK_included */
