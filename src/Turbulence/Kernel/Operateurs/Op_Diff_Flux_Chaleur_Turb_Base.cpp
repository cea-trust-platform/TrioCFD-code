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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Op_Diff_Flux_Chaleur_Turb_Base.cpp
// Directory:   $TURBULENCE_ROOT/src/Kernel/Operateurs
//
//////////////////////////////////////////////////////////////////////////////

#include <Op_Diff_Flux_Chaleur_Turb_Base.h>
#include <Discretisation_base.h>

// GF je ne sais pas si cette classe sert vraiment, mais elle herite pas de Op_Diff_base
Implemente_base(Op_Diff_Flux_Chaleur_Turb_Base,"Op_Diff_Flux_Chaleur_Turb_Base",Operateur_base);
Implemente_instanciable(Op_Diff_Flux_Chaleur_Turb_negligeable,"Op_Diff_Flux_Chaleur_Turb_negligeable",Op_Diff_Flux_Chaleur_Turb_Base);
Implemente_instanciable(Op_Diff_Flux_Chaleur_Turb,"Op_Diff_Flux_Chaleur_Turb",OWN_PTR(Op_Diff_Flux_Chaleur_Turb_Base));


////  printOn
//

Sortie& Op_Diff_Flux_Chaleur_Turb_Base::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}

Sortie& Op_Diff_Flux_Chaleur_Turb_negligeable::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}

Sortie& Op_Diff_Flux_Chaleur_Turb::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}

//// readOn
//

Entree& Op_Diff_Flux_Chaleur_Turb_Base::readOn(Entree& s )
{
  return s ;
}

Entree& Op_Diff_Flux_Chaleur_Turb_negligeable::readOn(Entree& s )
{
  return s ;
}

Entree& Op_Diff_Flux_Chaleur_Turb::readOn(Entree& s )
{
  Operateur::lire(s);
  return s;
}


void Op_Diff_Flux_Chaleur_Turb::typer()
{
  if (typ=="negligeable")
    {
      OWN_PTR(Op_Diff_Flux_Chaleur_Turb_Base)::typer("Op_Diff_Flux_Chaleur_Turb_negligeable");
      valeur().associer_diffusivite_turbulente();
    }
  else
    {
      Nom nom_type="Op_Diff_Flux_Chaleur_Turb_";
      nom_type +=equation().discretisation().que_suis_je();
      nom_type += "_";
      Nom type_inco=equation().inconnue().que_suis_je();
      nom_type+=(type_inco.suffix("Champ_"));
      if (axi)
        nom_type += "_Axi";
      OWN_PTR(Op_Diff_Flux_Chaleur_Turb_Base)::typer(nom_type);
      valeur().associer_eqn(equation());
      valeur().associer_diffusivite_turbulente();
      Cerr << valeur().que_suis_je() << finl;

    }
}




