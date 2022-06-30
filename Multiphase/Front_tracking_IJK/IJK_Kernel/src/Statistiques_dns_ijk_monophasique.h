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
// File      : Statistiques_dns_ijk_monophasique.h
// Directory : $IJK_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////
#ifndef Statistiques_dns_ijk_monophasique_H
#define Statistiques_dns_ijk_monophasique_H
#include <FixedVector.h>
#include <IJK_Field.h>



#include <Statistiques_dns_ijk.h>
#include <TRUSTArrays.h>
class IJK_Grid_Geometry;
class Statistiques_dns_ijk_monophasique : public Statistiques_dns_ijk
{
  Declare_instanciable(Statistiques_dns_ijk_monophasique);
public:
  using Statistiques_dns_ijk::initialize;
  void initialize(const IJK_Grid_Geometry&) override;
  void update_stat(const FixedVector<IJK_Field_double, 3>& vitesse,
                   const IJK_Field_double& pression,
                   const IJK_Field_double& temperature,
                   double dt);
};
#endif
