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
// File:        Neumann_paroi_rayo_semi_transp_VEF.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src/VEF
// Version:     /main/15
//
//////////////////////////////////////////////////////////////////////////////

#include <Neumann_paroi_rayo_semi_transp_VEF.h>
#include <Modele_rayo_semi_transp.h>
#include <Milieu_base.h>
#include <Champ_Uniforme.h>
#include <Front_VF.h>

Implemente_instanciable(Neumann_paroi_rayo_semi_transp_VEF,"Paroi_flux_impose_rayo_semi_transp_VEF",Neumann_paroi);


Sortie& Neumann_paroi_rayo_semi_transp_VEF::printOn(Sortie& os) const
{
  return os;
}

Entree& Neumann_paroi_rayo_semi_transp_VEF::readOn(Entree& is)
{
  is >>  le_champ_front;
  return is;
}

void Neumann_paroi_rayo_semi_transp_VEF::mettre_a_jour(double temps)
{
  calculer_temperature_bord(temps);
}

double Neumann_paroi_rayo_semi_transp_VEF::flux_impose(int i) const
{
  const DoubleTab& flux_radiatif = modele().flux_radiatif(frontiere_dis().le_nom()).valeurs();

  if (le_champ_front.valeurs().size()==1)
    return le_champ_front(0,0)-flux_radiatif(i,0);
  else if (le_champ_front.valeurs().dimension(1)==1)
    return le_champ_front(i,0)-flux_radiatif(i,0);
  else
    Cerr << "Neumann_paroi_rayo_semi_transp_VEF::flux_impose erreur" << finl;

  return 0;
}


double Neumann_paroi_rayo_semi_transp_VEF::flux_impose(int i,int j) const
{
  const DoubleTab& flux_radiatif = modele().flux_radiatif(frontiere_dis().le_nom()).valeurs();

  if (le_champ_front.valeurs().dimension(0)==1)
    return le_champ_front(0,j)-flux_radiatif(i);
  else
    return le_champ_front(i,j)-flux_radiatif(i);
}


const Cond_lim_base& Neumann_paroi_rayo_semi_transp_VEF::la_cl() const
{
  return (*this);
}


void Neumann_paroi_rayo_semi_transp_VEF::calculer_temperature_bord(double temps)
{
  // Cerr<<"Neumann_paroi_rayo_semi_transp_VEF::calculer_temperature_bord() : Debut"<<finl;
  // Pour une discretisation VEF, la temperature est codee sur les faces
  // des elements, il suffit donc ici  de recuperer les valeurs des temperatures
  // sur les faces de bord.
  DoubleTab& temperature = domaine_Cl_dis().equation().inconnue().valeurs();
  DoubleTab& tab = temperature_bord_->valeurs_au_temps(temps);

  const Front_VF& front_vf = ref_cast(Front_VF,frontiere_dis());
  int ndeb = front_vf.num_premiere_face();
  int nfin = ndeb + front_vf.nb_faces();
  int face=0;
  for(face=ndeb; face<nfin; face++)
    tab(face-ndeb,0) = temperature(face);
  //  Cerr<<"Neumann_paroi_rayo_semi_transp_VEF::calculer_temperature_bord() : Fin"<<finl;
}


int Neumann_paroi_rayo_semi_transp_VEF::compatible_avec_eqn(const Equation_base& eqn) const
{
  Motcle dom_app=eqn.domaine_application();
  Motcle Thermique="Thermique";
  Motcle indetermine="indetermine";
  if ( (dom_app==Thermique) || (dom_app==indetermine) )
    return 1;
  else
    {
      err_pas_compatible(eqn);
      return 0;
    }
}

void Neumann_paroi_rayo_semi_transp_VEF::completer()
{
  Neumann_paroi::completer();

  // On type et on dimmensionne le champ_front temperature_bord_
  const Front_VF& front_vf = ref_cast(Front_VF,frontiere_dis());
  int nb_comp = 1;

  temperature_bord_.typer("Champ_front_fonc");
  temperature_bord_->fixer_nb_comp(nb_comp);
  DoubleTab& tab= temperature_bord_->valeurs();
  tab.resize(front_vf.nb_faces(),nb_comp);
}
