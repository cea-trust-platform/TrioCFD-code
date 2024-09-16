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
// File:        Source_Reaction_Particules.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/6
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Source_Reaction_Particules_included
#define Source_Reaction_Particules_included

#include <Source_base.h>
#include <TRUST_Ref.h>

class Transport_Marqueur_FT;

/*! @brief Classe Source_Reaction_Particules Calcul du terme source a ajouter dans Navier_Stokes pour prendre en compte
 *
 *    l action en retour des particules sur le fluide
 *    La somme des forces (action fluide sur les particules) a ete calculee aux positions
 *    des particules et stockee dans source_stockage_ porte par Transport_Marqueur_FT
 *    Dans chacune des classes filles on procede aux deux etapes suivantes :
 *       - une interpolation aux sommets de source_stockage_
 *        (chaque particule p situee dans un element elem contribue a chacun des sommets de elem)
 *         - l evaluation aux faces du terme source a ajouter dans Navier_Stokes
 *          (pour une face donnee, on ajoute la contribution de la source evaluee en chacun
 *           des sommets de la face)
 *
 *
 */
class Source_Reaction_Particules : public Source_base
{
  Declare_base(Source_Reaction_Particules);
public :

  void mettre_a_jour(double temps) override;
  DoubleTab& calculer(DoubleTab& ) const override;
  void associer_eq_marqueur(const Equation_base& equation);

protected:

  void associer_domaines(const Domaine_dis_base& ,const Domaine_Cl_dis_base& ) override;
  void associer_pb(const Probleme_base& ) override;


  REF(Transport_Marqueur_FT) eq_marqueur_; //Reference a l equation de type Transport_Marqueur_FT
  //permet d acceder source_stockage_
};

#endif
