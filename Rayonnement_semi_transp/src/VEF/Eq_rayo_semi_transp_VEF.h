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
// File:        Eq_rayo_semi_transp_VEF.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src/VEF
// Version:     /main/15
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Eq_rayo_semi_transp_VEF_included
#define Eq_rayo_semi_transp_VEF_included

#include <Equation_rayonnement_base.h>


/////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//    classe Eq_rayo_semi_transp
//    Cette classe represente l'equation de rayonnement pour l'irradiance
//    dans un milieu semi transparent.
//    elle est associee au modele de rayonnement semi transparent
//    Elle definit la methode resoudre et calcule l'irradiance.
//
// .SECTION voir aussi
//    Equation_rayonnement
//////////////////////////////////////////////////////////////////////


class Eq_rayo_semi_transp_VEF: public Equation_rayonnement_base
{
  Declare_instanciable(Eq_rayo_semi_transp_VEF);

public:

  void modifier_matrice();
  void resoudre(double temps);
  void evaluer_cl_rayonnement(double temps);
  void completer();
  void typer_op_grad();
  void assembler_matrice();

  int nb_colonnes_tot();
  int nb_colonnes();
};


#endif
