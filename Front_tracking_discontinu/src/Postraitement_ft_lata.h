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
// File:        Postraitement_ft_lata.h
// Directory:   $TRUST_ROOT/src/Front_tracking_discontinu
// Version:     /main/11
//
//////////////////////////////////////////////////////////////////////////////
#ifndef Postraitement_ft_lata_included
#define Postraitement_ft_lata_included

#include <Postraitement_lata.h>
#include <Ref_Transport_Interfaces_FT_Disc.h>

class Motcle;
class Maillage_FT_Disc;
class Fichier_Lata;
class Comm_Group;

class Postraitement_ft_lata : public Postraitement_lata
{
  Declare_instanciable_sans_constructeur(Postraitement_ft_lata);
public:
  Postraitement_ft_lata();
  void set_param(Param& param);
  int lire_motcle_non_standard(const Motcle&, Entree&);
  void postraiter(int forcer);
protected:
  void lire_champ_interface(Entree&);
  void ecrire_maillage(Fichier_Lata&, const Maillage_FT_Disc&,
                       int fichiers_multiples) const;

  // L'equation de transport contenant les interfaces a postraiter
  REF(Transport_Interfaces_FT_Disc) refequation_interfaces;
  // Quels champs d'interface faut-il postraiter ?
  Motcles liste_champs_i_aux_sommets;
  Motcles liste_champs_i_aux_elements;
};

#endif
