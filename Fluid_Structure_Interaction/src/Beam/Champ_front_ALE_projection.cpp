/****************************************************************************
* Copyright (c) 2021, CEA
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
// File      : Champ_front_ALE_projection.cpp
// Directory : $BEAM_MODEL_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Champ_front_ALE_projection.h>
#include <Champ_front_ALE.h>
#include <Zone.h>
#include <Frontiere_dis_base.h>
#include <MD_Vector_tools.h>

Implemente_instanciable( Champ_front_ALE_projection, "Champ_front_ALE_projection", Ch_front_var_instationnaire_dep) ;

Sortie& Champ_front_ALE_projection::printOn( Sortie& os ) const
{

  return os;
}

Entree& Champ_front_ALE_projection::readOn( Entree& is )
{
  int dim;
  is >> dim;
  fixer_nb_comp(dim);

  fxyzt.dimensionner(dim);


  for (int i = 0; i<dim; i++)
    {
      Nom tmp;
      is >> tmp;
      Cerr << "Lecture et interpretation de la fonction " << tmp << finl;
      fxyzt[i].setNbVar(4);
      fxyzt[i].setString(tmp);
      fxyzt[i].addVar("x");
      fxyzt[i].addVar("y");
      fxyzt[i].addVar("z");
      fxyzt[i].addVar("t");
      fxyzt[i].parseString();
      Cerr << "Interpretation de la fonction " << tmp << " Ok" << finl;
    }

  return is;
}


double Champ_front_ALE_projection::evaluate(double temps, double x, double y, int dir)
{

  fxyzt[dir].setVar("x",x);
  fxyzt[dir].setVar("y",y);
  fxyzt[dir].setVar("t",temps);
  return fxyzt[dir].eval();

}
double Champ_front_ALE_projection::evaluate(double temps, double x, double y, double z, int dir)
{
  fxyzt[dir].setVar("x",x);
  fxyzt[dir].setVar("y",y);
  fxyzt[dir].setVar("z",z);
  fxyzt[dir].setVar("t",temps);
  return fxyzt[dir].eval();
}

