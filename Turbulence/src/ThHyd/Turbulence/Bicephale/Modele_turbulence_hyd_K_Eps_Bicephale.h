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
// File:        Modele_turbulence_hyd_K_Eps_Bicephale.h
// Directory:   $TRUST_ROOT/src/ThHyd/Turbulence
// Version:     /main/20
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modele_turbulence_hyd_K_Eps_Bicephale_included
#define Modele_turbulence_hyd_K_Eps_Bicephale_included

#include <Transport_K_ou_Eps.h>
#include <Modele_Fonc_Bas_Reynolds.h>
//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//    Classe Modele_turbulence_hyd_K_Eps_Bicephale
//    Cette classe represente le modele de turbulence (k,eps) pour les
//    equations de Navier-Stokes ou les 2 equations de k et eps sont gerees separement du point de vue informatique.
// .SECTION voir aussi
//    Mod_turb_hyd_base Mod_turb_hyd_ss_maille
//////////////////////////////////////////////////////////////////////////////

class Modele_turbulence_hyd_K_Eps_Bicephale : public Mod_turb_hyd_RANS_Bicephale
{

  Declare_instanciable(Modele_turbulence_hyd_K_Eps_Bicephale);

public:

  void set_param(Param& param);
  int lire_motcle_non_standard(const Motcle&, Entree&);
  int preparer_calcul();
  void verifie_loi_paroi();
  virtual bool initTimeStep(double dt);
  void mettre_a_jour(double );
  virtual inline Champ_Inc& K();
  virtual inline const Champ_Inc& K() const;
  virtual inline Champ_Inc& Eps();
  virtual inline const Champ_Inc& Eps() const;

  inline int nombre_d_equations() const;
  inline Transport_K_ou_Eps_base& eqn_transp_K();
  inline const Transport_K_ou_Eps_base& eqn_transp_K() const;
  inline Transport_K_ou_Eps_base& eqn_transp_Eps();
  inline const Transport_K_ou_Eps_base& eqn_transp_Eps() const;
  const Equation_base& equation_k_eps(int) const ;


  inline Modele_Fonc_Bas_Reynolds& associe_modele_fonction();
  inline const Modele_Fonc_Bas_Reynolds& associe_modele_fonction() const;

  virtual const Champ_base& get_champ(const Motcle& nom) const;
  virtual void get_noms_champs_postraitables(Noms& nom,Option opt=NONE) const;
protected:
  Modele_Fonc_Bas_Reynolds mon_modele_fonc;
  Transport_K_ou_Eps  eqn_transport_K;
  Transport_K_ou_Eps  eqn_transport_Eps;
  virtual Champ_Fonc& calculer_viscosite_turbulente(double temps);

};


// Description:
//    Renvoie le champ inconnue K du modele de turbulence
//    Cette inconnue est portee
//    par l equation de transport K porte par le modele.
//    (version const)
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Champ_Inc&
//    Signification: le champ inconnue (K)
//    Contraintes: reference constante
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l objet
inline const Champ_Inc& Modele_turbulence_hyd_K_Eps_Bicephale::K() const
{
  return eqn_transport_K.inconnue();
}


// Description:
//    Renvoie le champ inconnue K du modele de turbulence
//    Cette inconnue est portee
//    par l equation de transport K porte par le modele.
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Champ_Inc&
//    Signification: le champ inconnue (K)
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l objet
inline Champ_Inc& Modele_turbulence_hyd_K_Eps_Bicephale::K()
{
  return eqn_transport_K.inconnue();
}

// Description:
//    Renvoie le champ inconnue epsilon du modele de turbulence
//    Cette inconnue est portee
//    par l equation de transport epsilon porte par le modele.
//    (version const)
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Champ_Inc&
//    Signification: le champ inconnue (epsilon)
//    Contraintes: reference constante
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l objet
inline const Champ_Inc& Modele_turbulence_hyd_K_Eps_Bicephale::Eps() const
{
  return eqn_transport_Eps.inconnue();
}


// Description:
//    Renvoie le champ inconnue epsilon du modele de turbulence
//    Cette inconnue est portee
//    par l equation de transport epsilon porte par le modele.
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Champ_Inc&
//    Signification: le champ inconnue (epsilon)
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l objet
inline Champ_Inc& Modele_turbulence_hyd_K_Eps_Bicephale::Eps()
{
  return eqn_transport_Eps.inconnue();
}

// Description:
//    Renvoie l equation d evolution de K du modele de turbulence
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Transport_K_ou_Eps&
//    Signification: equation (K)
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l objet
inline Transport_K_ou_Eps_base& Modele_turbulence_hyd_K_Eps_Bicephale::eqn_transp_K()
{
  return eqn_transport_K;
}

// Description:
//    Renvoie l equation d evolution de K du modele de turbulence
//    (version const)
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Transport_K_ou_Eps&
//    Signification: equation (K)
//    Contraintes: reference constante
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l objet
inline const Transport_K_ou_Eps_base& Modele_turbulence_hyd_K_Eps_Bicephale::eqn_transp_K() const
{
  return eqn_transport_K;
}

// Description:
//    Renvoie l equation d evolution de epsilon du modele de turbulence
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Transport_K_ou_Eps&
//    Signification: equation (epsilon)
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l objet
inline Transport_K_ou_Eps_base& Modele_turbulence_hyd_K_Eps_Bicephale::eqn_transp_Eps()
{
  return eqn_transport_Eps;
}

// Description:
//    Renvoie l equation d evolution de epsilon du modele de turbulence
//    (version const)
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Transport_K_ou_Eps&
//    Signification: equation (epsilon)
//    Contraintes: reference constante
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l objet
inline const Transport_K_ou_Eps_base& Modele_turbulence_hyd_K_Eps_Bicephale::eqn_transp_Eps() const
{
  return eqn_transport_Eps;
}



inline Modele_Fonc_Bas_Reynolds& Modele_turbulence_hyd_K_Eps_Bicephale::associe_modele_fonction()
{
  return mon_modele_fonc;
}

inline const Modele_Fonc_Bas_Reynolds& Modele_turbulence_hyd_K_Eps_Bicephale::associe_modele_fonction() const
{
  return  mon_modele_fonc;
}
inline int Modele_turbulence_hyd_K_Eps_Bicephale::nombre_d_equations() const
{
  return 2;
}
#endif
