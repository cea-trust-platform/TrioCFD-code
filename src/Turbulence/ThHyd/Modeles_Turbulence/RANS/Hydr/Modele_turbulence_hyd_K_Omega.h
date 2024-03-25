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
// File:        Modele_turbulence_hyd_K_Omega.h
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/RANS/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modele_turbulence_hyd_K_Omega_included
#define Modele_turbulence_hyd_K_Omega_included

#include <TRUSTTabs_forward.h>
#include <Transport_K_Omega.h>


/*! @brief Classe Modele_turbulence_hyd_K_Omega Cette classe represente le modele de turbulence (k, omega) pour les
 *
 *     equations de Navier-Stokes.
 *
 * @sa Modele_turbulence_hyd_base Modele_turbulence_hyd_LES_base
 */
class Modele_turbulence_hyd_K_Omega: public Modele_turbulence_hyd_RANS_K_Omega_base, public Modele_turbulence_hyd_RANS_Gen<Modele_turbulence_hyd_K_Omega>
{
  Declare_instanciable(Modele_turbulence_hyd_K_Omega);
public:
  void set_param(Param& param) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  int preparer_calcul() override;
  void verifie_loi_paroi() override;
  bool initTimeStep(double dt) override;
  void mettre_a_jour(double) override;
  virtual inline Champ_Inc& K_Omega();
  virtual inline const Champ_Inc& K_Omega() const;

  inline Transport_K_Omega_base& eqn_transp_K_Omega() override;
  inline const Transport_K_Omega_base& eqn_transp_K_Omega() const override;

  void init_F1_F2_enstrophy();

  inline const DoubleTab& get_blenderF1() const;
  inline DoubleTab& get_blenderF1();
  inline const DoubleTab& get_fieldF2() const;
  inline DoubleTab& get_fieldF2();
  inline const DoubleTab& get_enstrophy() const;
  inline DoubleTab& get_enstrophy();

  const Equation_base& equation_k_omega(int i) const override
  {
    assert((i == 0));
    return eqn_transport_K_Omega_;
  }

  void controler() { eqn_transport_K_Omega_.controler_K_Omega(); }
  virtual Champ_Fonc& calculer_viscosite_turbulente(double );

protected:
  Transport_K_Omega eqn_transport_K_Omega_;
  void fill_turbulent_viscosity_tab(const int n, const DoubleTab& , DoubleTab& );
  DoubleTab blenderF1_; // Blending field for SST model
  DoubleTab fieldF2_; // for the turbulent viscosity in the SST model
  DoubleTab enstrophy_; // for the turbulent viscosity in the SST model
};

/*! @brief Renvoie le champ inconnue du modele de turbulence i.
 *
 * e. : (K, Omega). Cette inconnue est portee
 *     par l'equation de transport K-Omega porte par le modele.
 *     (version const)
 *
 * @return (Champ_Inc&) le champ inconnue (K, Omega)
 */
inline const Champ_Inc& Modele_turbulence_hyd_K_Omega::K_Omega() const
{
  return eqn_transport_K_Omega_.inconnue();
}


/*! @brief Renvoie le champ inconnue du modele de turbulence i.
 *
 * e. : (K, Omega). Cette inconnue est portee
 *     par l'equation de transport K-Omega porte par le modele.
 *
 * @return (Champ_Inc&) le champ inconnue (K, Omega)
 */
inline Champ_Inc& Modele_turbulence_hyd_K_Omega::K_Omega()
{
  return eqn_transport_K_Omega_.inconnue();
}

/*! @brief Renvoie l'equation du modele de turbulence i.
 *
 * e. : (K, Omega).
 *
 * @return (Transport_K_Omega&) equation (K, Omega)
 */
inline Transport_K_Omega_base& Modele_turbulence_hyd_K_Omega::eqn_transp_K_Omega()
{
  return eqn_transport_K_Omega_;
}

/*! @brief Renvoie l'equation du modele de turbulence i.
 *
 * e. : (K, Omega).
 *     (version const)
 *
 * @return (Transport_K_Omega&) equation (K, Omega)
 */
inline const Transport_K_Omega_base& Modele_turbulence_hyd_K_Omega::eqn_transp_K_Omega() const
{
  return eqn_transport_K_Omega_;
}

/*! @brief Returns the blending field F1.
 *
 *     (version const)
 *
 * @return (DoubleTab&) blenderF1
 */
inline const DoubleTab& Modele_turbulence_hyd_K_Omega::get_blenderF1() const
{
  return blenderF1_;
}

/*! @brief Returns the blending field F1.
 *
 * @return (DoubleTab&) blenderF1
 */
inline DoubleTab& Modele_turbulence_hyd_K_Omega::get_blenderF1()
{
  return blenderF1_;
}

/*! @brief Returns the field F2.
 *
 *     (version const)
 *
 * @return (DoubleTab&) blenderF2
 */
inline const DoubleTab& Modele_turbulence_hyd_K_Omega::get_fieldF2() const
{
  return fieldF2_;
}

/*! @brief Returns the field F2.
 *
 * @return (DoubleTab&) fieldF2
 */
inline DoubleTab& Modele_turbulence_hyd_K_Omega::get_fieldF2()
{
  return fieldF2_;
}

/*! @brief Returns the field enstrophy.
 *
 *     (version const)
 *
 * @return (DoubleTab&) enstrophy
 */
inline const DoubleTab& Modele_turbulence_hyd_K_Omega::get_enstrophy() const
{
  return enstrophy_;
}

/*! @brief Returns the field F2.
 *
 * @return (DoubleTab&) fieldF2
 */
inline DoubleTab& Modele_turbulence_hyd_K_Omega::get_enstrophy()
{
  return enstrophy_;
}

#endif
