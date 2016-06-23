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
// File:        Sauvegarde_Reprise_Maillage_FT.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/6
//
//////////////////////////////////////////////////////////////////////////////
#ifndef Sauvegarde_Reprise_Maillage_FT_included
#define Sauvegarde_Reprise_Maillage_FT_included
#include <arch.h>

class Maillage_FT_Disc;
class Zone_VF;
class Sortie;
class Entree;
class Domaine;
class Sauvegarde_Reprise_Maillage_FT
{
public:
  static void ecrire_xyz(const Maillage_FT_Disc& mesh, const Zone_VF& zone_vf, Sortie& fichier);
  static void lire_xyz(Maillage_FT_Disc& mesh,
                       const Zone_VF * zone_vf,
                       Entree * fichier,
                       const Domaine * domaine_src);
protected:

  int init_temps_physique(Maillage_FT_Disc& mesh);
  void   init_sommets(Maillage_FT_Disc& mesh, const Zone_VF * zone_vf);
  void   init_facettes(Maillage_FT_Disc& mesh, const Zone_VF * zone_vf);

};
#endif
