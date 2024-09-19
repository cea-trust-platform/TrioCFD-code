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
// File:        Modele_rayo_semi_transp.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src
// Version:     /main/16
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modele_rayo_semi_transp_included
#define Modele_rayo_semi_transp_included

#include <Equation_rayonnement_base.h>
#include <Probleme_base.h>
#include <Champ_front.h>


class Nom;
class Probleme_base;

/*! @brief Le Modele_rayo_semi_transp est un Probleme_base qui a 4 particularites : * Son equation doit etre typee en fonction de la dicretisation.
 *
 *     Cela impose de differer certaines initialisations jusqu'a
 *     connaitre la discretisation utilisee.
 *   * Il partage son domaine avec un probleme de type hydraulique
 *     En l'etat actuel, toute la geometrie est clonee (donc postraitee 2 fois!)
 *     Apres le travail de B. Mathieu sur la geometrie, ce ne sera plus necessaire.
 *   * Il n'y a qu'une seule valeur temporelle (futur=present).
 *     Il faudrait en faire un probleme independant du temps.
 *   * Il conserve une ref sur le probleme hydraulique. Cette ref est utilisee de
 *     maniere intensive.
 *
 *
 * @sa Pb_Couple_rayo_semi_transp Equation_rayonnement_base
 */
class Modele_rayo_semi_transp: public Probleme_base
{
  Declare_instanciable(Modele_rayo_semi_transp);

public:

  ///////////////////////////////////////////////
  //                                           //
  // Implementation de l'interface de Problem  //
  //                                           //
  ///////////////////////////////////////////////

  void terminate() override  {  finir(); }

  double computeTimeStep(bool& stop) const override
  {
    stop=false;
    return DMAXFLOAT;
  }
  bool initTimeStep(double dt) override;
  bool iterateTimeStep(bool& converged) override;
  void validateTimeStep() override;

  ////////////////////////////////////////////////////////
  //                                                    //
  // Fin de l'implementation de l'interface de Problem  //
  //                                                    //
  ////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////
  //                                                 //
  // Implementation de l'interface de Probleme_base  //
  //                                                 //
  /////////////////////////////////////////////////////

  void completer() override { }
  int nombre_d_equations() const override   { return 1; }
  const Equation_base& equation(int i) const override;
  Equation_base& equation(int i) override;
  const Equation_base& get_equation_by_name(const Nom&) const override;
  Equation_base& getset_equation_by_name(const Nom&) override;
  double calculer_pas_de_temps() const override  {  return DMAXFLOAT;  }

  // Cette methode ne doivent pas servir : on passe par l'interface de Problem
  void mettre_a_jour(double temps) override
  {
    exit();
  }

  bool is_pb_rayo() override { return true ; }

  void preparer_calcul() override;
  void discretiser(Discretisation_base&) override;
  void associer_sch_tps_base(const Schema_Temps_base&) override;

  //////////////////////////////////////////////////////////////
  //                                                          //
  // Fin de l'implementation de l'interface de Probleme_base  //
  //                                                          //
  //////////////////////////////////////////////////////////////

  inline void associer_probleme(Probleme_base& Pb);
  Champ_Inc_base& put_irradience();
  inline Probleme_base& probleme();
  inline const Probleme_base& probleme() const;
  inline Equation_rayonnement_base& eq_rayo();
  inline const Equation_rayonnement_base& eq_rayo() const;
  inline const double& valeur_sigma() const;
  const Champ_front& flux_radiatif(const Nom& nom_bord) const;
  void calculer_flux_radiatif();

protected :
  REF(Probleme_base) mon_probleme_;
  OWN_PTR(Equation_rayonnement_base) Eq_rayo_;
  static const double sigma;
};


inline Equation_rayonnement_base& Modele_rayo_semi_transp::eq_rayo()
{
  assert(Eq_rayo_.non_nul());
  return Eq_rayo_.valeur();
}

inline const Equation_rayonnement_base& Modele_rayo_semi_transp::eq_rayo() const
{
  assert(Eq_rayo_.non_nul());
  return Eq_rayo_.valeur();
}

inline const Equation_base& Modele_rayo_semi_transp::equation(int i) const
{
  assert(i==0);
  return Eq_rayo_;
}

inline Equation_base& Modele_rayo_semi_transp::equation(int i)
{
  assert(i==0);
  return Eq_rayo_;
}

inline const Equation_base& Modele_rayo_semi_transp::get_equation_by_name(const Nom& un_nom) const
{
  assert(Motcle(un_nom)==Motcle("Eq_rayo_semi_transp"));
  return Eq_rayo_;
}

inline Equation_base& Modele_rayo_semi_transp::getset_equation_by_name(const Nom& un_nom)
{
  assert(Motcle(un_nom)==Motcle("Eq_rayo_semi_transp"));
  return Eq_rayo_;
}

inline const double& Modele_rayo_semi_transp::valeur_sigma() const
{
  return sigma;
}

inline void Modele_rayo_semi_transp::associer_probleme(Probleme_base& Pb)
{
  mon_probleme_ = Pb;
}

inline Probleme_base& Modele_rayo_semi_transp::probleme()
{
  return mon_probleme_.valeur();
}

inline const Probleme_base& Modele_rayo_semi_transp::probleme() const
{
  return mon_probleme_.valeur();
}

#endif
