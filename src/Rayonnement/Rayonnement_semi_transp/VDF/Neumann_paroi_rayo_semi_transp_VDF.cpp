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
// File:        Neumann_paroi_rayo_semi_transp_VDF.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src/VDF
// Version:     /main/20
//
//////////////////////////////////////////////////////////////////////////////

#include <Neumann_paroi_rayo_semi_transp_VDF.h>
#include <Champ_front_uniforme.h>
#include <Modele_rayo_semi_transp.h>
#include <Schema_Temps_base.h>
#include <Fluide_Incompressible.h>
#include <Champ_Uniforme.h>
#include <Domaine_VDF.h>
#include <communications.h>
Implemente_instanciable(Neumann_paroi_rayo_semi_transp_VDF,"Paroi_flux_impose_rayo_semi_transp_VDF",Neumann_paroi);


/*! @brief
 *
 * @param (Sortie& os) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Neumann_paroi_rayo_semi_transp_VDF::printOn(Sortie& os) const
{
  return os;
}

/*! @brief Lecture des parametres de la condition Neumann_paroi Lecture de l'emissivite de la paroi
 *
 *     Lecture du coefficient A
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Neumann_paroi_rayo_semi_transp_VDF::readOn(Entree& is)
{
  is >>  le_champ_front;
  return is;
}


double Neumann_paroi_rayo_semi_transp_VDF::flux_impose(int i) const
{
  const Domaine_VDF& zvdf = ref_cast(Domaine_VDF,domaine_Cl_dis().domaine_dis());
  const IntTab& face_voisins = zvdf.face_voisins();
  const Front_VF& front_vf = ref_cast(Front_VF,frontiere_dis());
  int ndeb = front_vf.num_premiere_face();

  int elem = face_voisins(ndeb+i,0);
  if (elem == -1)
    elem = face_voisins(ndeb+i,1);

  double signe = 1;
  const DoubleTab& flux_radiatif = modele().flux_radiatif(frontiere_dis().le_nom()).valeurs();

  if (le_champ_front->valeurs().size()==1)
    return signe * le_champ_front->valeurs()(0,0)-flux_radiatif(i,0);
  else if (le_champ_front->valeurs().dimension(1)==1)
    return signe * (le_champ_front->valeurs()(i,0) - flux_radiatif(i,0));
  else
    Cerr << "Neumann_paroi_rayo_semi_transp_VDF::flux_impose erreur" << finl;

  return 0;
}


double Neumann_paroi_rayo_semi_transp_VDF::flux_impose(int i,int j) const
{
  const Domaine_VDF& zvdf = ref_cast(Domaine_VDF,domaine_Cl_dis().domaine_dis());
  const IntTab& face_voisins = zvdf.face_voisins();
  const Front_VF& front_vf = ref_cast(Front_VF,frontiere_dis());
  int ndeb = front_vf.num_premiere_face();

  int elem = face_voisins(ndeb+i,0);
  if (elem == -1)
    elem = face_voisins(ndeb+i,1);

  const DoubleTab& flux_radiatif = modele().flux_radiatif(frontiere_dis().le_nom()).valeurs();

  if (le_champ_front->valeurs().dimension(0)==1)
    return le_champ_front->valeurs()(0,j)-flux_radiatif(i);
  else
    return le_champ_front->valeurs()(i,j)-flux_radiatif(i);
}


const Cond_lim_base& Neumann_paroi_rayo_semi_transp_VDF::la_cl() const
{
  return (*this);
}

void Neumann_paroi_rayo_semi_transp_VDF::calculer_temperature_bord(double temps)
{
  const Domaine_VDF& zvdf = ref_cast(Domaine_VDF,domaine_Cl_dis().domaine_dis());
  const IntTab& face_voisins=zvdf.face_voisins();
  const Milieu_base& le_milieu = mon_dom_cl_dis->equation().milieu();
  ////const Champ_Uniforme& Lambda = ref_cast(Champ_Uniforme,le_milieu.conductivite().valeur());
  const Front_VF& front_vf = ref_cast(Front_VF,frontiere_dis());
  int nb_faces = front_vf.nb_faces();
  int ndeb = front_vf.num_premiere_face();
  const DoubleTab& rho = le_milieu.masse_volumique()->valeurs();
  const DoubleTab& Cp = le_milieu.capacite_calorifique()->valeurs();
  const DoubleTab& Lambda= le_milieu.conductivite()->valeurs();

  assert(le_milieu.capacite_calorifique()->nb_comp() == 1);
  assert(le_milieu.masse_volumique()->nb_comp() == 1);
  double d_Cp, d_rho, d_Lambda;
  if (Cp.get_md_vector().non_nul() || rho.get_md_vector().non_nul() || Lambda.get_md_vector().non_nul())
    {
      // L'un des champs n'est pas uniforme
      ArrOfDouble tmp(3);
      tmp[0] = local_max_vect(Cp);
      tmp[1] = local_max_vect(rho);
      tmp[2] = local_max_vect(Lambda);
      mp_max_for_each_item(tmp);
      d_Cp = tmp[0];
      d_rho = tmp[1];
      d_Lambda = tmp[2];
    }
  else
    {
      // Tous les champs sont uniforme. Raccourci car mp_max penalisant
      d_Cp = Cp(0,0);
      d_rho = rho(0,0);
      d_Lambda = Lambda(0,0);
    }

  Schema_Temps_base& sch = mon_dom_cl_dis->equation().probleme().schema_temps();
  double dt= sch.pas_de_temps() ;


  //  Cerr<<"Nom du bord : "<<front_vf.le_nom()<<finl;
  ////double d_Lambda = Lambda(0,0);

  const DoubleTab& T_f = mon_dom_cl_dis->equation().inconnue().valeurs();
  DoubleTab& T_b = temperature_bord_->valeurs_au_temps(temps);

  const DoubleTab& flux_radiatif = modele().flux_radiatif(frontiere_dis().le_nom()).valeurs();

  int face=0;
  int num_face;
  for (face=0; face<nb_faces; face++)
    {
      num_face = face+ndeb;
      double eF = zvdf.dist_norm_bord(num_face);
      int elem = face_voisins(num_face,0);
      if(elem < 0)
        elem = face_voisins(num_face,1);
      // Pour eviter division par 0 au 1er pas de temps:
      // double omega=1./(1.+eF*eF*d_rho*d_Cp/d_Lambda/dt);
      double omega=d_Lambda*dt/(d_Lambda*dt+eF*eF*d_rho*d_Cp);
      //      Cerr << "omega=" << omega << finl;
      if (sub_type(Champ_front_uniforme,le_champ_front.valeur()))
        {
          T_b(face,0) = omega*(T_f(elem) + (le_champ_front->valeurs()(0,0) - flux_radiatif(face,0))/(d_Lambda/eF))
                        + (1-omega) * T_b(face,0);
        }
      else
        T_b(face,0) = omega*(T_f(elem) + (le_champ_front->valeurs()(face,0) - flux_radiatif(face,0))/(d_Lambda/eF))
                      + (1-omega) * T_b(face,0);
    }
}


int Neumann_paroi_rayo_semi_transp_VDF::compatible_avec_eqn(const Equation_base& eqn) const
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


void Neumann_paroi_rayo_semi_transp_VDF::completer()
{
  Neumann_paroi::completer();

  // On type et on dimmensionne le champ_front temperature_bord_
  //  const Milieu_base& mil=mon_dom_cl_dis->equation().milieu();
  const Front_VF& front_vf = ref_cast(Front_VF,frontiere_dis());
  int nb_comp = 1;
  temperature_bord_.typer("Champ_front_fonc");
  temperature_bord_->fixer_nb_comp(nb_comp);
  DoubleTab& tab= temperature_bord_->valeurs();
  tab.resize(front_vf.nb_faces(),nb_comp);

  // On initialise le tableau des temperatures de bord egale a la
  // temperature initiale du milieu courant
  const Domaine_VDF& zvdf = ref_cast(Domaine_VDF,domaine_Cl_dis().domaine_dis());
  int ndeb = front_vf.num_premiere_face();
  const IntTab& face_voisins = zvdf.face_voisins();
  const DoubleTab& T = mon_dom_cl_dis->equation().inconnue().valeurs();
  int face=0;
  //
  // Debut de la boucle sur les faces de bord
  //

  for(face=0; face<front_vf.nb_faces(); face++)
    {
      int elem = face_voisins(face+ndeb,0);
      if (elem<0)
        elem = face_voisins(face+ndeb,1);
      tab(face,0) = T(elem);
    }
}
