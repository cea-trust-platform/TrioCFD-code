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
// File:        Algorithmes_Transport_FT_Disc.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/7
//
//////////////////////////////////////////////////////////////////////////////
#ifndef Algorithmes_Transport_FT_Disc_included
#define Algorithmes_Transport_FT_Disc_included
/////////////////////////////////////////////////////////////////////////////
//
// File                : Algorithmes_Transport_FT_Disc.h
// Directory        : $TRUST_ROOT/Front_Tracking_Discontinu
//
//////////////////////////////////////////////////////////////////////////////

#include <Objet_U.h>
#include <Deriv.h>

class Zone_VF;
class DoubleTab;
class ArrOfDouble;
class Algorithmes_Transport_FT_Disc;

Declare_deriv(Algorithmes_Transport_FT_Disc);

// .DESCRIPTION        : classe Algorithmes_Transport_FT_Disc
//
//  Classe regroupant des algorithmes utilises pour manipuler le
//  maillage dans l'equation de transport des interfaces.
//  Les methodes sont eventuellement specialisees VDF / VEF
class Algorithmes_Transport_FT_Disc : public Objet_U
{
  Declare_instanciable(Algorithmes_Transport_FT_Disc);
public:
  virtual void calculer_moments_indicatrice(const Zone_VF& zone_vf,
                                            const DoubleTab& indicatrice,
                                            double& volume_total,
                                            double& volume_phase_1,
                                            ArrOfDouble& barycentre,
                                            ArrOfDouble& barycentre_phase_1);

protected:

private:

};

#endif
// ifndef Algorithmes_Transport_FT_Disc_inclu
