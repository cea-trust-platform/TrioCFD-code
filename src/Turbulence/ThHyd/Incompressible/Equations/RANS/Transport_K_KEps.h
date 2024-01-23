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
// File:        Transport_K_KEps.h
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Incompressible/Equations/RANS
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Transport_K_KEps_included
#define Transport_K_KEps_included

#include <Transport_K_Eps_non_std.h>
#include <Loi_2couches.h>
#include <Operateur_Conv.h>
#include <TRUST_Ref.h>

class Modele_turbulence_hyd_K_Eps_2_Couches;
class Motcle;

/*! @brief classe Transport_K_KEps Cette classe represente l'equation de transport de l'energie cinetique
 *
 *     turbulente K et du taux de dissipation epsilon (eps) associee au modele
 *     de turbulence (k,eps) a deux couches.
 *     On traite en une seule equation le transport des deux
 *     grandeurs turbulentes. Il s'agit donc d'une equation vectorielle, dont
 *     le champ inconnue possede 2 composantes: K et epsilon.
 *     On instanciera un objet de cette classe uniquement si
 *     on utilise le modele a deux couches  pour traiter la turbulence
 *
 * @sa Equation_base
 */
class Transport_K_KEps : public Transport_K_Eps_non_std
{

  Declare_instanciable(Transport_K_KEps);

public :

  void set_param(Param& titi) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  void associer_milieu_base(const Milieu_base&) override;
  void associer_modele_turbulence(const Mod_turb_hyd_RANS_keps&) override;
  inline int get_nbcouches() const;
  inline int get_yswitch() const;
  inline int get_switch() const;
  inline int get_nutswitch() const;
  inline int get_impr() const;
  inline const Champ_Inc& vitesse_transportante();

  const Motcle& domaine_application() const override;
  Loi_2couches_base& loi2couches()
  {
    return loi_2couches.valeur();
  }

protected :

  Loi_2couches loi_2couches;

private:
  int nb_couches;        //le nombre de couches de mailles ou on applique le modele 1-equation
  int ystar_switch;        //la valeur de y* qui delimite les deux couches
  int type_switch;                // indique si on choisit un switch par y* ou par nu_t
  int nut_switch;        //la valeur de nut/nu qui delimite les deux couches.
  int impr;                //Indique si on doit afficher la zon,e des 2 couches dans le .out
};


/*! @brief Renvoie le nombre de couches utilisees pour le modele 1-equation.
 *
 * (version const)
 *
 * @return (int) le nombre de couches
 */
inline int Transport_K_KEps::get_nbcouches() const
{
  return nb_couches;
}

/*! @brief Renvoie le y* de switch entre les deux couches pour le modele a deux couches.
 *
 * (version const)
 *
 * @return (int) le nombre de couches
 */
inline int Transport_K_KEps::get_yswitch() const
{
  return ystar_switch;
}



/*! @brief Renvoie 0 si on choisit le switch par y*, 0 par nu_t.
 *
 * (version const)
 *
 * @return (int)
 */
inline int Transport_K_KEps::get_switch() const
{
  return type_switch;
}

/*! @brief Renvoie la valeur de nut/nu qui delimite les deux couches.
 *
 * (version const)
 *
 * @return (int) valeur limite de nu_t/nu
 */
inline int Transport_K_KEps::get_nutswitch() const
{
  return nut_switch;
}

/*! @brief indique si on doit ecrire le domaine des 2 couches dans le .
 *
 * out. (version const)
 *
 * @return (int) 1 si on imprime , 0 sinon.
 */
inline int Transport_K_KEps::get_impr() const
{
  return impr;
}


#endif


