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
// File:        Terme_Source_Constituant_Vortex_VEF_Face.h
// Directory:   $TRUST_ROOT/src/Front_tracking_discontinu
// Version:     /main/6
//
//////////////////////////////////////////////////////////////////////////////
#ifndef Terme_Source_Constituant_Vortex_VEF_Face_included
#define Terme_Source_Constituant_Vortex_VEF_Face_included
#include <Terme_Source_Constituant_VEF_Face.h>
#include <Senseur_Interface.h>

class Terme_Source_Constituant_Vortex_VEF_Face : public Terme_Source_Constituant_VEF_Face
{
  Declare_instanciable(Terme_Source_Constituant_Vortex_VEF_Face);
public:
  DoubleTab& calculer(DoubleTab& tab) const;
  DoubleTab& ajouter(DoubleTab& tab) const;
  void associer_pb (const Probleme_base& pb);

  void mettre_a_jour(double temps);
  void ajouter_terme_div_u(DoubleVect& secmem_pression, double dt) const;
protected:
  double rayon_spot_;
  ArrOfDouble delta_spot_;
  // si >=0. on renormalise le terme source discret pour que son integrale
  // de volume soit egal a integrale_.
  double integrale_;
  double debit_;
  Senseur_Interface senseur_;

  ArrOfDouble position_spot_;
};

#endif
