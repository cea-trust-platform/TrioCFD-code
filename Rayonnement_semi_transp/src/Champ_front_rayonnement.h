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
// File:        Champ_front_rayonnement.h
// Directory:   $TRUST_ROOT/src/Rayonnement_semi_transp
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Champ_front_rayonnement_included
#define Champ_front_rayonnement_included

#include <Champ_front_var_instationnaire.h>

//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//     classe Champ_front_rayonnement
//     Cette classe represente un champ sur une frontiere utilise pour
//     le rayonnement.
//     Elle ne fait rien, le tableau de valeurs est mis a jour par
//     Flux_radiatif_V?F::evaluer_cl_rayonnement
//     De ce fait, elle est inclassable dans la hierarchie... Eviter
//     si possible d'utiliser les champs de cette facon : le tableau
//     de valeurs ne devrait etre modifie que par les methodes
//     initialiser et mettre_a_jour (cf Champ_front_base).
//
// .SECTION voir aussi
//////////////////////////////////////////////////////////////////////////////

class Champ_front_rayonnement : public Champ_front_var_instationnaire
{
  Declare_instanciable(Champ_front_rayonnement);

public :

  Champ_front_base& affecter_(const Champ_front_base& ch);
  void mettre_a_jour(double);

};

#endif
