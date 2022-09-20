/****************************************************************************
* Copyright (c) 2019, CEA
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
// File      : MonofluidVar.h
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_MonofluidVar_included
#define IJK_MonofluidVar_included

#include <IJK_Field.h>


/*! @brief : class IJK_MonofluidVar
 *
 *  Cette classe vise à simplifier l'appel aux propriétés constantes par phase
 *  dans la formulation monofluide.
 *
 *
 *
 */
class IJK_MonofluidVar
{
public :
  IJK_MonofluidVar() {}
  ~IJK_MonofluidVar() {}
  void initialize(const double liqu_val, const double vap_val);

  ///////////////
  // Accessors //
  ///////////////

  double operator()(const double indic_l) const;
  double operator()(const IJK_Field_double& indic_l, const int i, const int j, const int k) const ;
  void operator()(const IJK_Field_double& indic_l, IJK_Field_double& res) const;
  double harmo(const double indic_l) const ;
  double harmo(const IJK_Field_double indic_l, const int i, const int j, const int k) const ;
  void harmo(const IJK_Field_double& indic_l, IJK_Field_double& res) const;

  ///////////////
  // Modifiers //
  ///////////////

  IJK_MonofluidVar& operator*= (const IJK_MonofluidVar& mono2);
  IJK_MonofluidVar& operator/= (const IJK_MonofluidVar& mono2);
  IJK_MonofluidVar& operator+= (const IJK_MonofluidVar& mono2);

  /////////////
  // getters //
  /////////////

  double liqu() const;
  double vap() const;

private:
  double liqu_val_;
  double vap_val_;
};


// Non member functions


IJK_MonofluidVar operator* (const IJK_MonofluidVar& mono1, const IJK_MonofluidVar& mono2);
IJK_MonofluidVar operator/ (const IJK_MonofluidVar& mono1, const IJK_MonofluidVar& mono2);
IJK_MonofluidVar operator+ (const IJK_MonofluidVar& mono1, const IJK_MonofluidVar& mono2);

#endif /* IJK_Energie_included */
