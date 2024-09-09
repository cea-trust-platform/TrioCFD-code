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
// File:        Op_Diff_Flux_Chaleur_Turb_Base.h
// Directory:   $TURBULENCE_ROOT/src/Kernel/Operateurs
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Op_Diff_Flux_Chaleur_Turb_Base_included
#define Op_Diff_Flux_Chaleur_Turb_Base_included

#include <Operateur_Diff_base.h>
#include <Operateur.h>
#include <Operateur_negligeable.h>
#include <Matrice_Morse.h>
#include <Op_Diff_Flux_Chaleur_Turb_Base.h>
#include <TRUST_Deriv.h>
#include <Champ_Fonc.h>
#include <Champ_Don.h>
#include <TRUSTTabs_forward.h>

/*! @brief class Op_Diff_Flux_Chaleur_Turb_Base Sert a modeliser le terme diffusif dans l'equation de transport
 *
 *  des fluctuations thermiques
 *  On traite les deux equations de transport en une seule
 *  equation vectorielle.
 *  La viscosite_turbulente est a diviser par la constante Prdt_teta
 *  pour la diffusion de K et par la constante Prdt_epsteta pour la
 *  diffusion de Eps
 *
 */
class Op_Diff_Flux_Chaleur_Turb_Base : public Operateur_base
{
  Declare_base(Op_Diff_Flux_Chaleur_Turb_Base);

public:
  virtual void associer_diffusivite_turbulente() =0;

};




class   Op_Diff_Flux_Chaleur_Turb_negligeable : public Operateur_negligeable, public Op_Diff_Flux_Chaleur_Turb_Base

{
  Declare_instanciable(Op_Diff_Flux_Chaleur_Turb_negligeable);

public:

  inline void associer(const Domaine_dis&, const Domaine_Cl_dis&, const Champ_Inc& ) override;
  inline DoubleTab& ajouter(const DoubleTab& ,  DoubleTab& ) const override;
  inline DoubleTab& calculer(const DoubleTab& , DoubleTab& ) const override;
  inline void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const override;
  inline void contribuer_au_second_membre(DoubleTab& ) const override;
  inline void modifier_pour_Cl(Matrice_Morse&, DoubleTab&) const override;
  inline void dimensionner(Matrice_Morse& ) const override;
  inline void associer_diffusivite_turbulente() override;
};


//////////////////////////////////////////////////////////////////////////////
//
// CLASS: Op_Diff_Flux_Chaleur_Turb
//
//////////////////////////////////////////////////////////////////////////////


class Op_Diff_Flux_Chaleur_Turb :  public Operateur,
  public DERIV(Op_Diff_Flux_Chaleur_Turb_Base)
{

  Declare_instanciable(Op_Diff_Flux_Chaleur_Turb);

public:

  inline Operateur_base& l_op_base() override;
  inline const Operateur_base& l_op_base() const override;
  inline void associer_diffusivite_turbulente();
  inline DoubleTab& ajouter(const DoubleTab& , DoubleTab& ) const override;
  inline DoubleTab& calculer(const DoubleTab& , DoubleTab& ) const override;
  void typer() override;
  inline int op_non_nul() const override;


};

///////////////////////////////////////////////////////////////
//   Fonctions inline de la classe Op_Diff_Flux_Chaleur_Turb_negligeable
///////////////////////////////////////////////////////////////

inline void Op_Diff_Flux_Chaleur_Turb_negligeable::associer(const Domaine_dis& domaine_dis,
                                                            const Domaine_Cl_dis& domaine_cl_dis,
                                                            const Champ_Inc& inco)
{
  Operateur_negligeable::associer(domaine_dis,domaine_cl_dis,inco);
}

inline DoubleTab& Op_Diff_Flux_Chaleur_Turb_negligeable::ajouter(const DoubleTab& x,  DoubleTab& y) const
{
  return Operateur_negligeable::ajouter(x,y);
}

inline DoubleTab& Op_Diff_Flux_Chaleur_Turb_negligeable::calculer(const DoubleTab& x, DoubleTab& y) const
{
  return Operateur_negligeable::calculer(x,y);
}

/*! @brief on assemble la matrice.
 *
 */
inline void Op_Diff_Flux_Chaleur_Turb_negligeable::contribuer_a_avec(const DoubleTab& inco,
                                                                     Matrice_Morse& matrice) const
{
  ;
}

/*! @brief on ajoute la contribution du second membre.
 *
 */
inline void Op_Diff_Flux_Chaleur_Turb_negligeable::contribuer_au_second_membre(DoubleTab& resu) const
{
  ;
}

// Modification des Cl
inline void  Op_Diff_Flux_Chaleur_Turb_negligeable::modifier_pour_Cl(Matrice_Morse& matrice, DoubleTab& resu) const
{
  ;
}

inline void  Op_Diff_Flux_Chaleur_Turb_negligeable::dimensionner(Matrice_Morse& matrice) const
{
  ;
}

inline void Op_Diff_Flux_Chaleur_Turb_negligeable::associer_diffusivite_turbulente()
{
  ;
}

/*! @brief
 *
 */
inline Operateur_base& Op_Diff_Flux_Chaleur_Turb::l_op_base()
{
  if(!non_nul())
    Cerr << "Op_Diff_Flux_Chaleur_Turb n'a pas ete typer" << finl;
  return valeur();
}

/*! @brief
 *
 */
inline const Operateur_base& Op_Diff_Flux_Chaleur_Turb::l_op_base() const
{
  if(!non_nul())
    Cerr << "Op_Diff_Flux_Chaleur_Turb n'a pas ete typer" << finl;
  return valeur();
}

inline void Op_Diff_Flux_Chaleur_Turb::associer_diffusivite_turbulente()
{
  valeur().associer_diffusivite_turbulente();
}

inline DoubleTab& Op_Diff_Flux_Chaleur_Turb::ajouter(const DoubleTab& inconnue, DoubleTab& resu) const
{
  return valeur().ajouter(inconnue, resu);
}


inline DoubleTab& Op_Diff_Flux_Chaleur_Turb::calculer(const DoubleTab& inconnue, DoubleTab& resu) const
{
  return valeur().calculer(inconnue, resu);
}

inline int Op_Diff_Flux_Chaleur_Turb::op_non_nul() const
{
  if (non_nul())
    return 1;
  else
    return 0;
}

#endif

