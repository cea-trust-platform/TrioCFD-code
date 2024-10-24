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
// File:        Paroi_rugueuse.cpp
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Incompressible/Cond_Lim
//
//////////////////////////////////////////////////////////////////////////////


#include <Paroi_rugueuse.h>
#include <Paroi_log_QDM.h>
#include <Param.h>

Implemente_instanciable(Paroi_rugueuse,"Paroi_rugueuse",Dirichlet_paroi_fixe);


/*! @brief Ecrit le type de l'objet sur un flot de sortie.
 *
 * @param (Sortie& s) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Paroi_rugueuse::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}


/*! @brief Simple appel a: Dirichlet_homogene::readOn(Entree& )
 *
 * @param (Entree& s) un flot d'entree
 * @return (Entree& s) le flot d'entree modifie
 */
Entree& Paroi_rugueuse::readOn(Entree& s )
{

  Dirichlet_paroi_fixe::readOn(s) ;
  Cout <<"Valeur par defaut de Erugu dans la loi standard :"<<E_DEF<<finl;
  Param param(que_suis_je());
  param.ajouter("erugu",&erugu_,Param::REQUIRED);

  param.lire_avec_accolades(s);
  Cout <<"Valeur utilisateur :"<< erugu_<<finl;
  return s;
}




