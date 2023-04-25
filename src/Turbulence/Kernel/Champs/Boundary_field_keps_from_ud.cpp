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
// File:        Boundary_field_keps_from_ud.cpp
// Directory:   $TURBULENCE_ROOT/src/Kernel/Champs
//
//////////////////////////////////////////////////////////////////////////////


#include <Boundary_field_keps_from_ud.h>
#include <Boundary_field_uniform_keps_from_ud.h>
#include <Param.h>
#include <Front_VF.h>
#include <Champ_front_uniforme.h>
#include <Champ_Inc_base.h>
#include <Equation_base.h>
#include <Schema_Temps_base.h>
#include <Champ_front_calc.h>

Implemente_instanciable(Boundary_field_keps_from_ud,"Boundary_field_keps_from_ud",Ch_front_var_instationnaire_dep);

Sortie& Boundary_field_keps_from_ud::printOn(Sortie& os) const
{
  return (os);
}

// XD Boundary_field_keps_from_ud front_field_base Boundary_field_keps_from_ud 1 To specify a K-Eps inlet field with hydraulic diameter, speed, and turbulence intensity (VDF only)
Entree& Boundary_field_keps_from_ud::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("U", &vitesse_, Param::REQUIRED); 	// XD_ADD_P front_field_base U 0 Initial velocity magnitude
  param.ajouter("D", &d, Param::REQUIRED); 		// XD_ADD_P double Hydraulic diameter
  param.ajouter("I", &I, Param::REQUIRED); 		// XD_ADD_P double Turbulence intensity [%]
  param.lire_avec_accolades(is);
  if (I<=0 || I>=1)
    {
      Cerr << "Turbulence intensity is a number between 0 and 1. For instance, 0.05 means 5%." << finl;
      Process::exit();
    }
  fixer_nb_comp(2);
  return (is);
}

int Boundary_field_keps_from_ud::initialiser(double tps, const Champ_Inc_base& inco)
{
  const Front_VF& le_bord = ref_cast(Front_VF,frontiere_dis());
  int nb_faces = le_bord.nb_faces();
  DoubleTab& tab=les_valeurs->valeurs();
  tab.resize(nb_faces, 2);

  if (!Ch_front_var_instationnaire_dep::initialiser(tps,inco)) return 0;

  //sch = inco.equation().schema_temps();
  int nb_cases=inco.equation().schema_temps().nb_valeurs_temporelles();
  // Define the boundary field for the velocity:
  Champ_front_base& vitesse_frontiere = vitesse_.valeur();
  vitesse_frontiere.associer_fr_dis_base(frontiere_dis());
  vitesse_frontiere.completer();
  vitesse_frontiere.fixer_nb_valeurs_temporelles(nb_cases);
  vitesse_frontiere.initialiser(tps, inco);

  // pour utiliser la trace ...
  if (sub_type(Champ_front_calc,vitesse_frontiere))
    ref_cast(Champ_front_calc,vitesse_frontiere).set_distant(0);

  vitesse_frontiere.mettre_a_jour(tps);

  mettre_a_jour(tps);
  return 1;
}

void Boundary_field_keps_from_ud::mettre_a_jour(double tps)
{
  // Update speed:
  Champ_front_base& vitesse_frontiere = vitesse_.valeur();
  vitesse_frontiere.changer_temps_futur(tps, 1);
  vitesse_frontiere.mettre_a_jour(tps); // OK

  // Compute KEps:
  //const DoubleTab& vit = vitesse_frontiere.valeurs(); // Se base sur temps_defaut qui reste a 0 !!
  const DoubleTab& vit = vitesse_frontiere.valeurs_au_temps(tps);
  DoubleTab& keps = les_valeurs->valeurs();
  int nb_faces = ref_cast(Front_VF,frontiere_dis()).nb_faces();
  for (int face=0; face<nb_faces; face++)
    {
      int i = (vit.dimension(0)==1 ? 0 : face); // Uniform field
      double u = (vit.line_size() == 1) ? std::fabs(vit(i, 0)) : sqrt(vit(i, 0)*vit(i, 0)+vit(i, 1)*vit(i, 1));

      pair val = k_eps_from_udi(u, d, I, dimension);
      keps(face, 0) = val.k;
      keps(face, 1) = val.eps;
      if (face==0)
        {
          Cerr << "Provisoire t=" << tps << " " << " U(t)=" << 2+tps << " u(t)=" << u << " k=" << val.k << " eps=" << val.eps << finl;
        }
    }
}


