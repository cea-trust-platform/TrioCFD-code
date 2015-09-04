/****************************************************************************
* Copyright (c) 2015, CEA
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
// File:        Source_rayo_semi_transp_base.h
// Directory:   $TRUST_ROOT/src/Rayonnement_semi_transp
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Source_rayo_semi_transp_base_included
#define Source_rayo_semi_transp_base_included

#include <Source_base.h>
#include <Ref_Modele_rayo_semi_transp.h>
class Modele_rayo_semi_transp;
//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//    classe Source_rayo_semi_transp_base
//    Classe de base de la hierarchie des termes sources de l'eqution de
//    temperature pour les problemes de rayonnement semi transparent.
// .SECTION voir aussi
//////////////////////////////////////////////////////////////////////////////
class Source_rayo_semi_transp_base : public Source_base
{
  Declare_base(Source_rayo_semi_transp_base);

public :

  inline Modele_rayo_semi_transp& Modele();
  inline const Modele_rayo_semi_transp& Modele() const;
  virtual void associer_modele_rayo(Modele_rayo_semi_transp& modele);
  void mettre_a_jour(double temps)
  {
    ;
  }

protected :

  REF(Modele_rayo_semi_transp) le_modele_;
};

inline Modele_rayo_semi_transp& Source_rayo_semi_transp_base::Modele()
{
  return le_modele_.valeur();
}

inline const Modele_rayo_semi_transp& Source_rayo_semi_transp_base::Modele() const
{
  return le_modele_.valeur();
}

#endif
