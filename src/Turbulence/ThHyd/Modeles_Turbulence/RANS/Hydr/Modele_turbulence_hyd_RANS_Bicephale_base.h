/****************************************************************************
* Copyright (c) 2017, CEA
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
// File:        Modele_turbulence_hyd_RANS_Bicephale_base.h
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/RANS/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modele_turbulence_hyd_RANS_Bicephale_base_included
#define Modele_turbulence_hyd_RANS_Bicephale_base_included

#include <Modele_turbulence_hyd_2_eq_base.h>

class Transport_K_ou_Eps_base;

/*! @brief Classe Modele_turbulence_hyd_RANS_Bicephale_base Classe de base des modeles de type RANS en formulation bicephale : les equations de k et epsilon sont gerees separement
 *
 * @sa Modele_turbulence_hyd_2_eq_base
 */
class Modele_turbulence_hyd_RANS_Bicephale_base: public Modele_turbulence_hyd_2_eq_base
{
  Declare_base_sans_constructeur(Modele_turbulence_hyd_RANS_Bicephale_base);
public:

  Modele_turbulence_hyd_RANS_Bicephale_base()
  {
    EPS_MIN_ = 1.e-10;
    K_MIN_ = 1.e-10;
  }

  void set_param(Param& param) override;
  void completer() override;
  int sauvegarder(Sortie& os) const override;
  int reprendre(Entree& is) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;

  const Champ_base& get_champ(const Motcle& nom) const override;
  void get_noms_champs_postraitables(Noms& nom, Option opt = NONE) const override;

  int nombre_d_equations() const override { return 2; }

  virtual Transport_K_ou_Eps_base& eqn_transp_K()=0;
  virtual Transport_K_ou_Eps_base& eqn_transp_Eps()=0;
  virtual const Transport_K_ou_Eps_base& eqn_transp_K() const =0;
  virtual const Transport_K_ou_Eps_base& eqn_transp_Eps() const =0;
  virtual const Equation_base& equation_k_eps(int) const =0;

  void associer_seconde_eqn(const Equation_base&);

  inline Equation_base& seconde_equation();
  inline const Equation_base& seconde_equation() const;

protected:
  REF(Equation_base) ma_seconde_equation_;
};

/*! @brief Renvoie la seconde equation associee au modele de turbulence en formulation bicephale.
 *
 * (c'est une equation du type Equation_base)
 *
 * @return (Equation_base&) la seconde equation associee au modele de turbulence en formulation bicephale
 */
inline Equation_base& Modele_turbulence_hyd_RANS_Bicephale_base::seconde_equation()
{
  if (ma_seconde_equation_.non_nul() == 0)
    {
      Cerr << "\nError in Modele_turbulence_hyd_RANS_Bicephale_base::seconde_equation() : The equation is unknown !" << finl;
      Process::exit();
    }
  return ma_seconde_equation_.valeur();
}

inline const Equation_base& Modele_turbulence_hyd_RANS_Bicephale_base::seconde_equation() const
{
  if (ma_seconde_equation_.non_nul() == 0)
    {
      Cerr << "\nError in Modele_turbulence_hyd_RANS_Bicephale_base::seconde_equation() : The equation is unknown !" << finl;
      Process::exit();
    }
  return ma_seconde_equation_.valeur();
}

#endif /* Modele_turbulence_hyd_RANS_Bicephale_base_included */
