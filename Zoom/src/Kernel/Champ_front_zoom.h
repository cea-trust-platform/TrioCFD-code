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
// File:        Champ_front_zoom.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Kernel
// Version:     /main/11
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Champ_front_zoom_included
#define Champ_front_zoom_included

#include <Ch_front_var_instationnaire_dep.h>
#include <Ref_Probleme_base.h>
#include <Ref_Pb_MG.h>
#include <Ref_Champ_Inc_base.h>
#include <Frontiere_dis_base.h>
class Motcle;
class Champ_front_base;
class Equation_base;
class Milieu_base;
class Domaine_dis_base;
class Domaine_Cl_dis_base;
class Champ_Inc_base;

/*! @brief classe Champ_front_zoom Classe derivee de Champ_front_var.
 *
 * Elle represente les
 *      champs a l interface d'un probleme deux grilles
 *      On distingue le probleme concerne par ce champ sur la frontiere et
 *      le probleme exterieure dont on veut prendre la trace de l'inconnue
 *
 * @sa Champ_front_var_instationnaire Champ_Inc
 */
class Champ_front_zoom : public Ch_front_var_instationnaire_dep
{

  Declare_instanciable(Champ_front_zoom);

public:
  void creer(Nom&, Nom&, Nom&, Nom&,  Motcle&);
  Champ_front_base& affecter_(const Champ_front_base&) override ;
  int initialiser(double temps, const Champ_Inc_base& inco) override;
  void mettre_a_jour(double temps) override;
  //  Nom& nature_cond_lim();
  Pb_MG& le_pb_MG();
  Probleme_base& le_pb_fin();
  Probleme_base& le_pb_courant();
  Probleme_base& le_pb_exterieur();
  const Champ_Inc_base& inconnue() const ;
  const Equation_base& equation() const;
  //  inline Nom& dirichlet_ou_neumann_();
  const Domaine_dis_base& domaine_dis() const;
  const Domaine_Cl_dis_base& domaine_Cl_dis() const;
  const Milieu_base& milieu() const;
  const Frontiere_dis_base& front_dis_exterieure() const;


protected :
  REF(Pb_MG) le_pb_MG_;
  REF(Probleme_base) le_pb_ext_;
  REF(Probleme_base) le_pb_courant_;
  REF(Champ_Inc_base) l_inconnueF;
  REF(Champ_Inc_base) l_inconnue;
  //  Nom dirichlet_ou_neumann;
  Nom nom_bord;
  Nom nom_pbMG_, nom_pbG_, nom_pbF_, bord_;
  Motcle nom_inco_;
};

#endif


/*Nom& Champ_front_zoom::dirichlet_ou_neumann_()
  {
  return dirichlet_ou_neumann;
  }*/

