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
// File:        Equation_rayonnement_base.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Equation_rayonnement_base_included
#define Equation_rayonnement_base_included

#include <Operateur_Diff.h>
#include <Matrice_Morse.h>
#include <Operateur_Grad.h>
#include <TRUST_Ref.h>

class Modele_rayo_semi_transp;
class Motcle;
class Milieu_base;
class Fluide_base;

class Equation_rayonnement_base: public Equation_base
{
  Declare_base_sans_constructeur(Equation_rayonnement_base);
public:

  Equation_rayonnement_base();
  void set_param(Param& titi) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  bool initTimeStep(double dt) override;
  virtual bool solve();
  void associer_milieu_base(const Milieu_base&) override;
  inline void associer_fluide(const Fluide_base&);
  void associer_modele_rayonnement(const Modele_rayo_semi_transp&);
  Milieu_base& milieu() override;
  const Milieu_base& milieu() const override;
  inline const Modele_rayo_semi_transp&  Modele() const;
  inline Modele_rayo_semi_transp& Modele();
  int nombre_d_operateurs() const override;
  const Operateur& operateur(int) const override;
  Operateur& operateur(int) override;
  const Champ_Inc_base& inconnue() const override;
  Champ_Inc_base& inconnue() override;
  void discretiser() override;
  inline Fluide_base& fluide();
  inline const Fluide_base& fluide() const;

  // pas de flux calcule correctement par les operateurs...
  inline int impr(Sortie& os) const override
  {
    return 1;
  };

  virtual int nb_colonnes_tot()=0;
  virtual int nb_colonnes()=0;
  void Mat_Morse_to_Mat_Bloc(Matrice& matrice_tmp);
  void dimensionner_Mat_Bloc_Morse_Sym(Matrice& matrice_tmp);


  const Discretisation_base& discretisation() const;

  Operateur_Grad& operateur_gradient();
  const Operateur_Grad& operateur_gradient() const;
  virtual void typer_op_grad()=0;

  void associer_pb_base(const Probleme_base& pb) override;

  virtual void resoudre(double temps)=0;
  virtual void assembler_matrice()=0;


  /////////////////////////////////////////////////////

  virtual void modifier_matrice()=0;
  virtual void evaluer_cl_rayonnement(double temps)=0;

  void completer() override;

protected:

  REF(Fluide_base) le_fluide;
  REF(Modele_rayo_semi_transp) le_modele;

  Operateur_Diff terme_diffusif;

  OWN_PTR(Champ_Inc_base) irradiance_;
  SolveurSys solveur;
  Matrice_Morse la_matrice;

  Operateur_Grad gradient;


};


/*! @brief Renvoie le modele de rayonnement semi transparent Precondition :
 *
 *  Parametre :
 *     Signification :
 *     Contraintes :
 *     Acces :
 *  Retour
 *     Signification :
 *     Contraintes :
 *  Exception :
 *  Effets de bord :
 *  Postcondition :
 *  Postcondition :
 *
 */
inline Modele_rayo_semi_transp& Equation_rayonnement_base::Modele()
{
  return le_modele.valeur();
}


/*! @brief Renvoie le modele de rayonnement semi transparent Precondition :
 *
 *  Parametre :
 *     Signification :
 *     Contraintes :
 *     Acces :
 *  Retour
 *     Signification :
 *     Contraintes : reference constante
 *  Exception :
 *  Effets de bord :
 *  Postcondition :
 *  Postcondition :
 *
 */
inline const Modele_rayo_semi_transp& Equation_rayonnement_base::Modele() const
{
  return le_modele.valeur();
}


/*! @brief Associe un fluide incompressible semi transparent a l'equation.
 *
 * @param (Fluide_base& un_fluide) le fluide incompressible semi transparent a associer
 */
inline void Equation_rayonnement_base::associer_fluide(const Fluide_base& un_fluide)
{
  le_fluide = un_fluide;
}


/*! @brief renvoie le fluide semi transparent associe a l'equation de rayonnement
 *
 * @return (le fluide associe a l'equation de rayonnement)
 */
inline const Fluide_base& Equation_rayonnement_base::fluide() const
{
  return le_fluide.valeur();
}


/*! @brief renvoie le fluide semi transparent associe a l'equation de rayonnement
 *
 * @return (le fluide associe a l'equation de rayonnement)
 */
inline Fluide_base& Equation_rayonnement_base::fluide()
{
  return le_fluide.valeur();
}


#endif
