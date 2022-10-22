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

#include <Schema_Euler_explicite_Phase_field.h>
#include <Mass_Redistribution_Phase_Field.h>
#include <Convection_Diffusion_Phase_field.h>
#include <Source_Con_Phase_field.h>
#include <Equation.h>
#include <Debog.h>

Implemente_instanciable(Schema_Euler_explicite_Phase_field,"Schema_Euler_explicite_Phase_field",Schema_Euler_explicite);

Sortie& Schema_Euler_explicite_Phase_field::printOn(Sortie& s) const { return Schema_Euler_explicite::printOn(s); }
Entree& Schema_Euler_explicite_Phase_field::readOn(Entree& s) { return Schema_Euler_explicite::readOn(s) ; }

int Schema_Euler_explicite_Phase_field::faire_un_pas_de_temps_eqn_base(Equation_base& eqn)
{
  DoubleTab& present = eqn.inconnue().valeurs(); // Un
  DoubleTab& futur   = eqn.inconnue().futur();   // Un+1
  DoubleTab dudt(futur);

  // Boundary conditions applied on Un+1:
  eqn.zone_Cl_dis()->imposer_cond_lim(eqn.inconnue(),temps_courant()+pas_de_temps());

  const Zone_VDF& zvdf = ref_cast(Zone_VDF, eqn.zone_dis().valeur());
  const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field, mon_probleme.valeur().equation(1));
  Sources list_sources = eq_c.sources();
  Source_Con_Phase_field& source_pf = ref_cast(Source_Con_Phase_field, list_sources(0).valeur());
  DoubleVect minnX = source_pf.get_minX();
  DoubleVect maxxX = source_pf.get_maxX();

  Mass_Redistribution_Phase_Field::impose_mass_redistribution(zvdf, present, minnX, maxxX);


  // On tourne la roue pour que les operateurs utilisent les champs au temps futur
  eqn.inconnue().avancer();
  eqn.derivee_en_temps_inco(dudt);
  eqn.inconnue().reculer();

  // Un+1=Un+dt_*dU/dt
  futur=dudt;
  futur*=dt_;
  futur+=present;

  Mass_Redistribution_Phase_Field::impose_mass_redistribution(zvdf, futur, minnX, maxxX);

  eqn.zone_Cl_dis()->imposer_cond_lim(eqn.inconnue(),temps_courant()+pas_de_temps());
  update_critere_statio(dudt, eqn);

  return 1;
}
