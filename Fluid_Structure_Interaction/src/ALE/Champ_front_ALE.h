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
// File:        Champ_front_ALE.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/ALE/src
// Version:     /main/10
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Champ_front_ALE_included
#define Champ_front_ALE_included

#include <Ch_front_var_instationnaire_dep.h>
#include <TRUST_Vector.h>
#include <Parser_U.h>
//////////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//     class Champ_front_ALE
//
//     Cette classe represente un champ sur une frontiere in case of ALE calculation.
//
// .SECTION voir aussi
//   Champ_Input_Proto
/////////////////////////////////////////////////////////////////////////////////
class Champ_front_ALE : public Ch_front_var_instationnaire_dep
{
  Declare_instanciable(Champ_front_ALE);

public:
  inline DoubleTab& get_vit_som_bord_ALE();
  Champ_front_base& affecter_(const Champ_front_base& ch) override;
  int initialiser(double temps, const Champ_Inc_base& inco) override;
  void mettre_a_jour(double temps) override;
  virtual void remplir_vit_som_bord_ALE(double);

protected :
  DoubleTab vit_som_bord_ALE;
  VECT(Parser_U) fxyzt;

private:


};

inline DoubleTab& Champ_front_ALE::get_vit_som_bord_ALE()
{
  return vit_som_bord_ALE;
}
#endif

