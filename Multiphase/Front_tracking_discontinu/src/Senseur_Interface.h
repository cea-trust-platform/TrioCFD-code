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
// File:        Senseur_Interface.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/2
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Senseur_Interface_included
#define Senseur_Interface_included
#include <Objet_U.h>
#include <Ref_Probleme_base.h>
#include <Ref_Equation_base.h>
#include <TRUSTArray.h>
class Probleme_base;

// Senseur permettant de trouver l'intersection entre un segment
//  fourni par l'utilisateur et une interface FT.
class Senseur_Interface : public Objet_U
{
  Declare_instanciable(Senseur_Interface);
public:
  virtual void associer_pb(const Probleme_base& pb)
  {
    probleme_ = pb;
  };
  virtual int calculer_position(ArrOfDouble& pos) const;
protected:
  Nom equation_interface_;
  REF(Probleme_base) probleme_;
  REF(Equation_base) equation_;
  ArrOfDouble segment_senseur_1_;
  ArrOfDouble segment_senseur_2_;
  // nombre de points tests du senseur (il doit y en avoir plus qu'un par maille)
  int nb_points_tests_;
};
#endif
