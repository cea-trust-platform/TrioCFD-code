/****************************************************************************
 * Copyright (c) 2015 - 2016, CEA
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 *modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *this list of conditions and the following disclaimer in the documentation
 *and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *may be used to endorse or promote products derived from this software without
 *specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 *FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBdoubleITUTE GOODS OR
 *SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, doubleRICT
 *LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 *OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 *DAMAGE.
 *
 *****************************************************************************/
/////////////////////////////////////////////////////////////////////////////
//
// File      : OpConvIJKQuickScalar_interface.h
// Directory : $IJK_ROOT/src/IJK/OpVDF
//
/////////////////////////////////////////////////////////////////////////////
#include <Corrige_flux_FT.h>
#include <IJK_Splitting.h>
#include <OpConvIJKQuickScalar.h>
#include <Operateur_IJK_base.h>
#include <Operateur_IJK_data_channel.h>

class Corrige_flux_FT;

class OpConvQuickInterfaceIJKScalar_double  : public OpConvQuickIJKScalar_double
{
  Declare_instanciable_sans_constructeur(OpConvQuickInterfaceIJKScalar_double);

public:
  OpConvQuickInterfaceIJKScalar_double() : OpConvQuickIJKScalar_double() {};

protected:
  // Ces methodes devraient cacher les methodes de OpConv..._base
  // Mais pas sur. Je pense que ca ne marche pas si on appelle
  // la methode de OpConv..base elle les trouvera qd meme.
  void compute_set(IJK_Field_double& dx) override;
  //   void compute_add(IJK_Field_double& dx) final;

public:
//  virtual void initialize(const IJK_Splitting& splitting,
//                          Corrige_flux_FT& corrige_flux);
//  virtual void initialize(const IJK_Splitting& splitting,
//  												const IJK_Field& indicatrice,
//                  				Corrige_flux_FT& corrige_flux);
  // using OpConvIJKQuickScalar_double::calculer;
  // using OpConvIJKQuickScalar_double::ajouter;
  // OpConvIJKQuickScalarInterface_double() :
  // correction_{}
  // {
  // }
  // void calculer(const IJK_Field_double& convect, const IJK_Field_double&
  // disc_prop, const IJK_Field_double& vx, const IJK_Field_double& vy, const
  // IJK_Field_double& vz, IJK_Field_double& res);
};
