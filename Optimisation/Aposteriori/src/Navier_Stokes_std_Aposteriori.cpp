/****************************************************************************
* Copyright (c) 2022, CEA
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

#include <Navier_Stokes_std_Aposteriori.h>
#include <Schema_Temps_base.h>
#include <Discret_Thyd.h>

Implemente_instanciable_sans_constructeur(Navier_Stokes_std_Aposteriori,"Navier_Stokes_std_Aposteriori",Navier_Stokes_std);

Sortie& Navier_Stokes_std_Aposteriori::printOn(Sortie& is) const { return Navier_Stokes_std::printOn(is); }
Entree& Navier_Stokes_std_Aposteriori::readOn(Entree& is) { return Navier_Stokes_std::readOn(is); }

void Navier_Stokes_std_Aposteriori::creer_champ(const Motcle& motlu)
{
  Navier_Stokes_std::creer_champ(motlu);

  if (motlu == "estimateur_aposteriori")
    {
      if (!estimateur_aposteriori.non_nul())
        {
          const Discret_Thyd& dis=ref_cast(Discret_Thyd, discretisation());
          dis.estimateur_aposteriori(zone_dis(), zone_Cl_dis(), la_vitesse, la_pression, diffusivite_pour_transport(), estimateur_aposteriori);
          champs_compris_.ajoute_champ(estimateur_aposteriori);
        }
    }
}

const Champ_base& Navier_Stokes_std_Aposteriori::get_champ(const Motcle& nom) const
{
  const double temps_init = schema_temps().temps_init();

  if (nom == "estimateur_aposteriori")
    {
      Champ_Fonc_base& ch = ref_cast_non_const(Champ_Fonc_base, estimateur_aposteriori.valeur());
      if (((ch.temps() != la_vitesse->temps()) || (ch.temps() == temps_init)) && (la_vitesse->mon_equation_non_nul()))
        ch.mettre_a_jour(la_vitesse->temps());
      return champs_compris_.get_champ(nom);
    }
  else
    return Navier_Stokes_std::get_champ(nom);
}
