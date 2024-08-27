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
// File:        Topologie_Maillage_FT.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Topologie_Maillage_FT_included
#define Topologie_Maillage_FT_included

#include <FTd_tools.h>
#include <TRUSTTabFT.h>
#include <Remailleur_Collision_FT_base.h>
#include <TRUSTTabs_forward.h>

class Maillage_FT_Disc;
class Domaine_VF;
class Motcle;
class Champ_base;
class Remaillage_FT;

/*! @brief : class Topologie_Maillage_FT Cette classe implemente les procedures de remaillage des interfaces pour le Front-Tracking :
 *
 *
 *
 * @sa Transport_Interfaces_FT_Disc Maillage_FT_Disc
 */

class Topologie_Maillage_FT : public Objet_U
{
  Declare_instanciable_sans_constructeur(Topologie_Maillage_FT);
public:
  Topologie_Maillage_FT();

  virtual void remailler_interface(const double temps,
                                   Maillage_FT_Disc& maillage,
                                   Champ_base& indicatrice,
                                   Remaillage_FT& algo_remaillage_local);

  virtual int calculer_composantes_connexes_pour_suppression(const Domaine_VF& domaine_vf,
                                                             const DoubleTab& indicatrice,
                                                             IntVect& num_compo) const;

  virtual double suppression_interfaces(const IntVect& num_compo,
                                        const ArrOfInt& flags_compo_a_supprimer,
                                        Maillage_FT_Disc& maillage,
                                        DoubleTab& indicatrice);

  int get_phase_continue() const;

protected:
  int test_RuptureCoalescenceInterfaces(const Maillage_FT_Disc& maillage);

  int test_collision_facettes(const Maillage_FT_Disc& maillage,
                              ArrOfInt& liste_elements_collision) const;

  int test_intersection_facettes_3D(int fa70, int fa71,
                                    const Maillage_FT_Disc& maillage) const;
  int test_intersection_facettes_2D(int fa70, int fa71,
                                    const Maillage_FT_Disc& maillage) const;

  int active_;

  OWN_PTR(Remailleur_Collision_FT_base) remailleur_Collision_;

  double Erreur_max_coordonnees_; // donnee par Parcours_interface
  int juric_local_;
  int phase_continue_;
};

#endif
