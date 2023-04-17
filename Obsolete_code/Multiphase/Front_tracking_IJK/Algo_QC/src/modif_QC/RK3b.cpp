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
// File      : RK3b.cpp
// Directory : $NEW_ALGO_QC_ROOT/src/modif_QC
//
/////////////////////////////////////////////////////////////////////////////

#include <RK3b.h>
#include <Equation.h>

Implemente_instanciable(RK3b,"Runge_Kutta_ordre_3_QC",Schema_Temps_base);


/*! @brief Simple appel a: Schema_Temps_base::printOn(Sortie& ) Ecrit le schema en temps sur un flot de sortie.
 *
 * @param (Sortie& s) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& RK3b::printOn(Sortie& s) const
{
  return  Schema_Temps_base::printOn(s);
}


/*! @brief Lit le schema en temps a partir d'un flot d'entree.
 *
 * Simple appel a: Schema_Temps_base::readOn(Entree& )
 *
 * @param (Entree& s) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& RK3b::readOn(Entree& s)
{
  return Schema_Temps_base::readOn(s) ;
}


/*! @brief Renvoie le nombre de valeurs temporelles a conserver.
 *
 * Ici : n et n+1, donc 2.
 *
 */
int RK3b::nb_valeurs_temporelles() const
{
  return 2 ;
}

/*! @brief Renvoie le nombre de valeurs temporelles futures.
 *
 * Ici : n+1, donc 1.
 *
 */
int RK3b::nb_valeurs_futures() const
{
  return 1 ;
}

/*! @brief Renvoie le le temps a la i-eme valeur future.
 *
 * Ici : t(n+1)
 *
 */
double RK3b::temps_futur(int i) const
{
  assert(i==1);
  return temps_courant()+pas_de_temps();
}

/*! @brief Renvoie le temps que doivent rendre les champs a l'appel de valeurs()
 *
 *     Ici : t(n+1)
 *
 */
double RK3b::temps_defaut() const
{
  return temps_courant()+pas_de_temps();
}

/*! @brief Effectue un pas de temps de Runge Kutta d'ordre 3, sur l'equation passee en parametre.
 *
 *     Le schema de Runge Kutta  d'ordre 3
 *      (cas 7 de Williamson) s'ecrit :
 *      q1=h f(x0)
 *      x1=x0+b1 q1
 *      q2=h f(x1) +a2 q1
 *      x2=x1+b2 q2
 *      q3=h f(x2)+a3 q2
 *      x3=x2+b3 q3
 *       avec a1=0, a2=-5/9, a3=-153/128
 *                               b1=1/3, b2=15/16, b3=8/15
 *
 * @param (Equation_base& eqn) l'equation que l'on veut faire avancer d'un pas de temps
 * @return (int) renvoie toujours 1
 */
int RK3b::faire_un_pas_de_temps_eqn_base(Equation_base& eqn)
{
  Cerr<<" ce schema ne doit servir que pour le QC.... qui a recode le RK3" <<finl;
  exit();

  return 1;
}

