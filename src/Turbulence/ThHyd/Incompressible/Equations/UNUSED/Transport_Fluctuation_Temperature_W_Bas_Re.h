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
// File:        Transport_Fluctuation_Temperature_W_Bas_Re.h
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Incompressible/Equations/UNUSED
//
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

#ifndef Transport_Fluctuation_Temperature_W_Bas_Re_included
#define Transport_Fluctuation_Temperature_W_Bas_Re_included

#include <Fluide_Incompressible.h>
#include <Transport_Fluctuation_Temperature_W.h>
#include <TRUST_Ref.h>

class Modele_turbulence_scal_Fluctuation_Temperature_W_Bas_Re;

class Motcle;

class Transport_Fluctuation_Temperature_W_Bas_Re : public Transport_Fluctuation_Temperature_W
{

  Declare_instanciable(Transport_Fluctuation_Temperature_W_Bas_Re);

public :

  void set_param(Param& titi) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  void associer_modele_turbulence(const Modele_turbulence_scal_Fluctuation_Temperature_W_Bas_Re& );
  inline const Milieu_base& milieu() const override;
  inline Milieu_base& milieu() override;
  inline Modele_turbulence_scal_Fluctuation_Temperature_W_Bas_Re& modele_turbulence();
  inline const Modele_turbulence_scal_Fluctuation_Temperature_W_Bas_Re& modele_turbulence() const;
  const Operateur& operateur(int) const override;
  Operateur& operateur(int) override;
  // int a_pour_Champ_Inc(const Motcle&, REF(Champ_base)& ) const;

protected :
  REF(Modele_turbulence_scal_Fluctuation_Temperature_W_Bas_Re) mon_modele_Fluctu_Temp;
};


// Fonctions inline:


inline const Milieu_base& Transport_Fluctuation_Temperature_W_Bas_Re::milieu() const
{
  if(!le_fluide.non_nul())
    {
      Cerr << "Il semble que vous n'ayez pas associe "
           << "Transport_Fluctuation_Temperature_W_Bas_Re avec un fluide" << finl;
      exit();
    }
  return le_fluide.valeur();
}

inline Milieu_base& Transport_Fluctuation_Temperature_W_Bas_Re::milieu()
{
  if(!le_fluide.non_nul())
    {
      Cerr << "Il semble que vous n'ayez pas associe "
           << "Transport_Fluctuation_Temperature_W_Bas_Re avec un fluide" << finl;
      exit();
    }
  return le_fluide.valeur();
}



inline Modele_turbulence_scal_Fluctuation_Temperature_W_Bas_Re& Transport_Fluctuation_Temperature_W_Bas_Re::modele_turbulence()
{
  assert(mon_modele_Fluctu_Temp.non_nul());
  return mon_modele_Fluctu_Temp.valeur();
}

inline const Modele_turbulence_scal_Fluctuation_Temperature_W_Bas_Re& Transport_Fluctuation_Temperature_W_Bas_Re::modele_turbulence() const
{

  assert(mon_modele_Fluctu_Temp.non_nul());
  return mon_modele_Fluctu_Temp.valeur();
}


#endif


