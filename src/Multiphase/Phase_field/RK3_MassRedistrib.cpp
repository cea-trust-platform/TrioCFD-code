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

#include <Mass_Redistribution_Phase_Field.h>
#include <Convection_Diffusion_Phase_field.h>
#include <Source_Con_Phase_field.h>
#include <RK3_MassRedistrib.h>
#include <Equation_base.h>

Implemente_instanciable(RK3_MassRedistrib,"Runge_Kutta_ordre_3_MassRedistrib",RK3);

Sortie& RK3_MassRedistrib::printOn(Sortie& s) const { return  RK3::printOn(s); }
Entree& RK3_MassRedistrib::readOn(Entree& s) { return RK3::readOn(s) ; }

int RK3_MassRedistrib::faire_un_pas_de_temps_eqn_base(Equation_base& eqn)
{
  if ( nb_pas_dt()>=0 && nb_pas_dt()<=100 && facsec_ == 1 ) print_warning(100);

  const double a2 = -5. / 9., a3 = -153. / 128.;
  const double b1 = 1. / 3., b2 = 15. / 16., b3 = 8. / 15.;

  DoubleTab& xi=eqn.inconnue().valeurs();
  DoubleTab& xipls1=eqn.inconnue().futur();
  DoubleTab qi(xi) ;
  DoubleTab present(xi);

  const Domaine_VDF& zvdf = ref_cast(Domaine_VDF, eqn.domaine_dis());
  const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field, mon_probleme->equation(1));
  Sources& list_sources = ref_cast_non_const(Sources, eq_c.sources());
  Source_Con_Phase_field& source_pf = ref_cast(Source_Con_Phase_field, list_sources(0).valeur());
  DoubleVect minnX = source_pf.get_minX();
  DoubleVect maxxX = source_pf.get_maxX();

  Mass_Redistribution_Phase_Field::impose_mass_redistribution(zvdf, xi, minnX, maxxX);

  // q1 = f(x0)
  eqn.derivee_en_temps_inco(qi);

  // q1 = h f(x0)
  qi *= dt_;

  // x1 = x0 + b1 q1
  xi.ajoute(b1, qi);

  // xipls1 = f(x1)

  //Mass_Redistribution_Phase_Field::impose_mass_redistribution(zvdf, xi, minnX, maxxX);

  eqn.derivee_en_temps_inco(xipls1);

  // q2 = h f(x1) + a2 q1
  qi *= a2;
  qi.ajoute(dt_,xipls1);

  // x2 = x1 + b2 q2
  xi.ajoute(b2, qi);

  // xipls1 = f(x2)

  //Mass_Redistribution_Phase_Field::impose_mass_redistribution(zvdf, xi, minnX, maxxX);

  eqn.derivee_en_temps_inco(xipls1);

  // q3 = h f(x2) + a3 q2
  qi *= a3;
  qi.ajoute(dt_,xipls1);

  // x3 = x2 + b3 q3
  xi.ajoute(b3, qi);

  //Mass_Redistribution_Phase_Field::impose_mass_redistribution(zvdf, xi, minnX, maxxX);

  // futur = x3
  xipls1 = xi;

  Mass_Redistribution_Phase_Field::impose_mass_redistribution(zvdf, xipls1, minnX, maxxX);

  // xi = (x3-x0)/dt
  xi -= present;
  xi /= dt_;
  update_critere_statio(xi, eqn);

  //Mass_Redistribution_Phase_Field::impose_mass_redistribution(zvdf, xi, minnX, maxxX);

  // Update boundary condition on futur:
  eqn.domaine_Cl_dis().imposer_cond_lim(eqn.inconnue(),temps_courant()+pas_de_temps());
  xipls1.echange_espace_virtuel();

  /*  Mass_Redistribution_Phase_Field::impose_mass_redistribution(zvdf, xi, minnX, maxxX);*/

  // xi = x0;
  xi = present;
  return 1;
}
