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
// File:        Convection_Diffusion_Concentration_Turbulent_FT_Disc.h
// Directory:   $TrioCFD_ROOT/Front_tracking_discontinu/src
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Convection_Diffusion_Concentration_Turbulent_FT_Disc_included
#define Convection_Diffusion_Concentration_Turbulent_FT_Disc_included

#include <Convection_Diffusion_Concentration_Turbulent.h>
#include <Ref_Equation_base.h>

class ArrOfBit;
/*! @brief Cette equation corrige le champ de concentration pour tenir compte de la presence d'une interface.
 *
 * @sa Convection_Diffusion_Concentration
 */
class Convection_Diffusion_Concentration_Turbulent_FT_Disc : public Convection_Diffusion_Concentration_Turbulent
{
  Declare_instanciable_sans_constructeur(Convection_Diffusion_Concentration_Turbulent_FT_Disc);
public:

  Convection_Diffusion_Concentration_Turbulent_FT_Disc();
  void set_param(Param& titi) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  void     mettre_a_jour(double temps) override;

protected:
  void ramasse_miette_simple(double temps);
  void mettre_a_jour_chimie();
  void multiplier_valeurs_faces(const ArrOfBit marqueurs_faces, double facteur, double& integrale_avant, DoubleTab& tab);
  void marquer_faces_sous_domaine(const Nom& nom_sous_domaine,
                                  ArrOfBit& marqueur,
                                  int sans_faces_non_std_volume_nul) const;

  // Algorithme a utiliser pour tenir compte des interfaces:
  enum Options { RIEN = 0, RAMASSE_MIETTES_SIMPLE = 1 };
  Options option_;
  // Dans quelle phase la masse du constituant se trouve-t-elle (RAMASSE_MIETTES_SIMPLE)
  int phase_a_conserver_;

  // Quelle equation du probleme fournit l'indicatrice de phase ?
  REF(Equation_base) ref_eq_transport_;

  // Liste non vide => cette equation modelise le produit de la reaction avec les equations
  //  dont la liste contient les noms:
  Noms equations_source_chimie_;
  // Constantes cinetiques
  double constante_cinetique_;
  double constante_cinetique_nu_t_;
  Nom nom_equation_nu_t_;
  // Modele pour taux de reactions
  int modele_cinetique_;
  // Domaine de sortie ou on annule la concentration
  Nom nom_domaine_sortie_;
};
#endif
