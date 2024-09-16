/****************************************************************************
* Copyright (c) 2023, CEA
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
// File:        Op_Diff_K_Omega_base.h
// Directory:   $TURBULENCE_ROOT/src/Kernel/Operateurs
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Op_Diff_K_Omega_base_included
#define Op_Diff_K_Omega_base_included

#include <Operateur_negligeable.h>
#include <Operateur_Diff_base.h>
#include <TRUSTTabs_forward.h>
#include <Matrice_Morse.h>
#include <TRUST_Deriv.h>
#include <Champ_Fonc.h>
#include <TRUST_Ref.h>
#include <Operateur.h>
#include <Champ_Don.h>

class Champ_base;

// cAlan: 2023-01-30 : pardon pour ce vilain copier-coller, mais
// j'aimerais que ça marche avant de modifier toute la structure des Operateurs.


/*! @brief classe Op_Diff_K_Omega_base Classe de base de la hierarchie des operateurs de diffusion
 *
 *     pour l'equation de transport K-omega.
 *     Sert a modeliser le terme diffusif dans l'equation de transport
 *     de K_Omega. On traite les deux equations de transport en une seule
 *     equation vectorielle.
 *     La viscosite_turbulente est a diviser par la constante Prdt_K
 *     pour la diffusion de K et par la constante Prdt_Omega pour la
 *     diffusion de Omega
 *
 * @sa Operateur_base Op_Diff_K_Omega Op_Diff_K_Omega_negligeable, Classe abstraite, Methode abstraite:, void associer_diffusivite_turbulente()
 */
class Op_Diff_K_Omega_base : public Operateur_base
{

  Declare_base(Op_Diff_K_Omega_base);

public:

  virtual void associer_diffusivite_turbulente() =0;
  inline virtual void associer_diffusivite(const Champ_base& )
  {
    Cerr<<get_info()->name()<<" doit coder associer_diffusivite"<<finl;
    Process::exit();
  }

  inline  virtual const Champ_base& diffusivite() const
  {
    Cerr<<get_info()->name()<<" doit coder diffusivite()"<<finl;
    throw;

  }

protected:

  bool diffuse_k_seul;
  bool diffuse_eps_seul;

};


/*! @brief Classe Op_Diff_K_Omega_negligeable Cette classe represente un operateur de diffusion pour
 *
 *     l'equation de transport K-omega negligeable.
 *     Lorsqu'un operateur de ce type est utilise dans une equation
 *     transport K-omega cela revient a negliger le terme de diffusion.
 *
 * @sa Op_Diff_K_Omega_base Operateur_negligeable
 */
class Op_Diff_K_Omega_negligeable : public Operateur_negligeable,
  public Op_Diff_K_Omega_base

{
  Declare_instanciable(Op_Diff_K_Omega_negligeable);

public:

  inline void associer(const Domaine_dis_base&, const Domaine_Cl_dis_base&, const Champ_Inc& ) override;
  inline DoubleTab& ajouter(const DoubleTab& ,  DoubleTab& ) const override;
  inline DoubleTab& calculer(const DoubleTab& , DoubleTab& ) const override;
  inline void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const override;
  inline void contribuer_au_second_membre(DoubleTab& ) const override;
  inline void modifier_pour_Cl(Matrice_Morse&, DoubleTab&) const override;
  inline void dimensionner(Matrice_Morse& ) const override;
  inline void associer_diffusivite_turbulente() override;

protected:
  REF(Champ_base) la_diffusivite;
};



/*! @brief Classe Op_Diff_K_Omega Classe generique de la hierarchie des operateurs de diffusion dans
 *
 *     l'equation de transport K-omega. Un objet Op_Diff_K_Omega peut
 *     referencer n'importe quel objet derivant de Op_Diff_K_Omega_base.
 *
 * @sa Op_Diff_K_Omega_base
 */
class Op_Diff_K_Omega : public Operateur,
  public OWN_PTR(Op_Diff_K_Omega_base)
{

  Declare_instanciable(Op_Diff_K_Omega);

public:

  inline Operateur_base& l_op_base() override;
  inline const Operateur_base& l_op_base() const override;
  inline void associer_diffusivite_turbulente();
  DoubleTab& ajouter(const DoubleTab& , DoubleTab& ) const override;
  DoubleTab& calculer(const DoubleTab& , DoubleTab& ) const override;
  inline Champ_base& associer_diffusivite(const Champ_base&);
  inline const Champ_base& diffusivite() const;
  void typer() override;
  inline int op_non_nul() const override;

protected :

  REF(Champ_base) la_diffusivite;
};


/*! @brief Associe divers objets a un operateurs negligeable: NE FAIT RIEN Simple appel a Operateur_negligeable::associer(const Domaine_dis_base&,
 *
 *                                                      const Domaine_Cl_dis_base&,
 *                                                      const Champ_Inc&)
 *
 * @param (Domaine_dis_base& domaine_dis)
 * @param (Domaine_Cl_dis_base& domaine_cl_dis)
 * @param (Champ_Inc& inco)
 */
inline void Op_Diff_K_Omega_negligeable::associer(const Domaine_dis_base& domaine_dis,
                                                  const Domaine_Cl_dis_base& domaine_cl_dis,
                                                  const Champ_Inc& inco)
{
  Operateur_negligeable::associer(domaine_dis,domaine_cl_dis,inco);
}



/*! @brief Ajoute la contribution de l'operateur a un tableau passe en parametre.
 *
 * Simple appel a Operateur_negligeable::ajouter(const DoubleTab&,DoubleTab&)
 *
 * @param (DoubleTab& x) le tableau sur lequel on applique l'operateur
 * @param (DoubleTab& y) tableau auquel on ajoute la contribution de l'operateur
 * @return (DoubleTab&) le parametre d'entree y non modifie
 */
inline DoubleTab& Op_Diff_K_Omega_negligeable::ajouter(const DoubleTab& x,  DoubleTab& y) const
{
  return Operateur_negligeable::ajouter(x,y);
}


/*! @brief Initialise le parametre tableau avec la contribution de l'operateur negligeable: initialise le tableau a ZERO.
 *
 *     Simple appel a Operateur_negligeable::(calculer(const DoubleTab&, DoubleTab&)
 *
 * @param (DoubleTab& x) le tableau sur lequel on applique l'operateur
 * @param (DoubleTab& y) tableau dans lequel stocke la contribution de l'operateur
 * @return (DoubleTab&) le tableau d'entree y mis a zero
 */
inline DoubleTab& Op_Diff_K_Omega_negligeable::calculer(const DoubleTab& x, DoubleTab& y) const
{
  return Operateur_negligeable::calculer(x,y);
}


/*! @brief NE FAIT RIEN
 *
 */
inline void Op_Diff_K_Omega_negligeable::contribuer_a_avec(const DoubleTab& inco,
                                                           Matrice_Morse& matrice) const
{
  ;
}

/*! @brief on ajoute la contribution du second membre.
 *
 */
inline void Op_Diff_K_Omega_negligeable::contribuer_au_second_membre(DoubleTab& resu) const
{
  ;
}

// Modification des Cl
inline void  Op_Diff_K_Omega_negligeable::modifier_pour_Cl(Matrice_Morse& matrice, DoubleTab& resu) const
{
  ;
}

inline void  Op_Diff_K_Omega_negligeable::dimensionner(Matrice_Morse& matrice) const
{
  ;
}

inline void Op_Diff_K_Omega_negligeable::associer_diffusivite_turbulente()
{
  ;
}

/*! @brief Renvoie l'objet sous-jacent upcaste en Operateur_base
 *
 * @return (Operateur_base&) l'objet sous-jacent upcaste en Operateur_base
 * @throws l'operateur n'a pas ete type
 */
inline Operateur_base& Op_Diff_K_Omega::l_op_base()
{
  if(!non_nul())
    Cerr << "Op_Diff_K_Omega n'a pas ete typer" << finl;
  return valeur();
}

/*! @brief Renvoie l'objet sous-jacent upcaste en Operateur_base (version const)
 *
 * @return (Operateur_base&) l'objet sous-jacent upcaste en Operateur_base
 * @throws l'operateur n'a pas ete type
 */
inline const Operateur_base& Op_Diff_K_Omega::l_op_base() const
{
  if(!non_nul())
    Cerr << "Op_Diff_K_Omega n'a pas ete typer" << finl;
  return valeur();
}


/*! @brief Appel a l'objet sous-jacent.
 *
 */
inline void Op_Diff_K_Omega::associer_diffusivite_turbulente()
{
  valeur().associer_diffusivite_turbulente();
}



inline int Op_Diff_K_Omega::op_non_nul() const
{
  if (non_nul())
    return 1;
  else
    return 0;
}


/*! @brief Appel a l'objet sous-jacent.
 *
 * Initialise le tableau passe en parametre avec la contribution
 *     de l'operateur.
 *
 * @param (DoubleTab& inconnue) tableau contenant les donnees sur lesquelles on applique l'operateur.
 * @param (DoubleTab& resu) tableau dans lequel stocke la contribution de l'operateur
 * @return (DoubleTab&) le tableau contenant le resultat
 */
inline const Champ_base& Op_Diff_K_Omega::diffusivite() const
{
  return la_diffusivite.valeur();
}


/*! @brief Associe la diffusivite a l'operateur.
 *
 * @param (Champ_Don& nu) le champ representant la diffusivite
 * @return (Champ_Don&) le champ representant la diffusivite
 */
inline Champ_base& Op_Diff_K_Omega::associer_diffusivite(const Champ_base& nu)
{
  la_diffusivite=nu;
  return la_diffusivite.valeur();
}




#endif
