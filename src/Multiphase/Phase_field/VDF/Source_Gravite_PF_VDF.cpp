/****************************************************************************
* Copyright (c) 2021, CEA
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
// File:        Source_Gravite_PF_VDF.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Multiphase/Phase_field/src/VDF
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Gravite_PF_VDF.h>
#include <Domaine_VDF.h>
#include <Domaine_Cl_VDF.h>
#include <Periodique.h>
#include <Dirichlet.h>
#include <Dirichlet_homogene.h>

#include <Navier_Stokes_phase_field.h>
#include <Convection_Diffusion_Phase_field.h>
#include <Probleme_base.h>
#include <Fluide_Incompressible.h>
#include <Source_Con_Phase_field.h>


Implemente_instanciable(Source_Gravite_PF_VDF,"Source_Gravite_PF_VDF",Source_base);


/*! @brief Imprime la source sur un flot de sortie.
 *
 * @param (Sortie& os) le flot de sortie pour l'impression
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Source_Gravite_PF_VDF::printOn(Sortie& os) const
{
  os <<que_suis_je()<< finl;
  return os;
}

/*! @brief Lecture de la source sur un flot d'entree.
 *
 * @param (Entree& is) le flot d'entree pour la lecture des parametres
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Source_Gravite_PF_VDF::readOn(Entree& is)
{
  return is;
}

/*! @brief Remplit le tableau volumes
 *
 * @param (Entree& is) le flot d'entree pour la lecture des parametres
 * @return le flot d'entree modifie
 */
void Source_Gravite_PF_VDF::associer_domaines(const Domaine_dis_base& dds,const Domaine_Cl_dis_base& domaine_cl)
{
  le_dom = ref_cast(Domaine_VDF,dds);
  le_dom_Cl = ref_cast(Domaine_Cl_VDF,domaine_cl);
}



/*! @brief Ajoute les termes sources
 *
 * @return (DoubleTab&)
 */
DoubleTab& Source_Gravite_PF_VDF::ajouter(DoubleTab& resu) const
{
  int face, nb_faces = le_dom->nb_faces();
  int premiere_face_interne = le_dom->premiere_face_int();

  const IntVect& orientation = le_dom->orientation();
  const DoubleVect& volumes_entrelaces = le_dom->volumes_entrelaces();
  DoubleVect porosite_surf ;          // porosites surfaciques
  porosite_surf.ref(equation().milieu().porosite_face());
  Navier_Stokes_phase_field& eq_NS_PF = ref_cast_non_const(Navier_Stokes_phase_field, equation());
  const DoubleVect& g = eq_NS_PF.get_g_();

  int num_cl;

  const int boussi = eq_NS_PF.get_boussi_();

  const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme->equation(1));
  Sources& list_sources = ref_cast_non_const(Sources, eq_c.sources());
  Source_Con_Phase_field& source_pf = ref_cast(Source_Con_Phase_field, list_sources(0).valeur());
  int type_systeme_naire = source_pf.get_type_systeme_naire(); //si type_systeme_naire = 0 ou 1

  if (type_systeme_naire==0)
    {
      if (boussi == 0)
        {
          for (num_cl=0 ; num_cl<le_dom->nb_front_Cl() ; num_cl++)
            {
              const Cond_lim& la_cl = le_dom_Cl->les_conditions_limites(num_cl);
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
              int ndeb = le_bord.num_premiere_face();
              int nfin = ndeb + le_bord.nb_faces();
              if (sub_type(Dirichlet,la_cl.valeur()) || sub_type(Dirichlet_homogene,la_cl.valeur()))
                ;
              else
                for (face=ndeb ; face<nfin ; face++)
                  {
                    resu(face) += g(orientation(face)) * porosite_surf(face) * volumes_entrelaces(face);
                  }
            }

          for (face=premiere_face_interne ; face<nb_faces; face++)
            {
              resu(face) += g(orientation(face)) * porosite_surf(face) * volumes_entrelaces(face);
            }
        }
      else if (boussi == 1)
        {
          const DoubleTab& rho = eq_NS_PF.rho().valeurs();
          double rho0 = eq_NS_PF.rho0();
          const IntTab& face_voisins = le_dom->face_voisins();
          const DoubleVect& volumes = le_dom->volumes();

          for (num_cl=0 ; num_cl<le_dom->nb_front_Cl() ; num_cl++)
            {
              const Cond_lim& la_cl = le_dom_Cl->les_conditions_limites(num_cl);
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
              int ndeb = le_bord.num_premiere_face();
              int nfin = ndeb + le_bord.nb_faces();
              if (sub_type(Dirichlet,la_cl.valeur()) || sub_type(Dirichlet_homogene,la_cl.valeur()))
                ;
              else
                for (face=ndeb ; face<nfin ; face++)
                  {
                    int el0=face_voisins(face,0);
                    int el1=face_voisins(face,1);
                    double vol0=volumes(el0);
                    double vol1=volumes(el1);
                    double cbetaface=(vol0*(rho(el0)/rho0-1.0)+vol1*(rho(el1)/rho0-1.0))/(vol0+vol1);
                    resu(face) += cbetaface * g(orientation(face)) * porosite_surf(face) * volumes_entrelaces(face);
                  }
            }

          for (face=premiere_face_interne ; face<nb_faces; face++)
            {
              int el0=face_voisins(face,0);
              int el1=face_voisins(face,1);
              double vol0=volumes(el0);
              double vol1=volumes(el1);
              double cbetaface=(vol0*(rho(el0)/rho0-1.0)+vol1*(rho(el1)/rho0-1.0))/(vol0+vol1);
              resu(face) += cbetaface * g(orientation(face)) * porosite_surf(face) * volumes_entrelaces(face);
            }
        }
    }
  else if (type_systeme_naire==1)
    {
      if (boussi == 0)
        {
          for (num_cl=0 ; num_cl<le_dom->nb_front_Cl() ; num_cl++)
            {
              const Cond_lim& la_cl = le_dom_Cl->les_conditions_limites(num_cl);
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
              int ndeb = le_bord.num_premiere_face();
              int nfin = ndeb + le_bord.nb_faces();
              if (sub_type(Dirichlet,la_cl.valeur()) || sub_type(Dirichlet_homogene,la_cl.valeur()))
                ;
              else
                for (face=ndeb ; face<nfin ; face++)
                  {
                    resu(face) += g(orientation(face)) * porosite_surf(face) * volumes_entrelaces(face);
                  }
            }

          for (face=premiere_face_interne ; face<nb_faces; face++)
            {
              resu(face) += g(orientation(face)) * porosite_surf(face) * volumes_entrelaces(face);
            }
        }
      else if (boussi == 1)
        {
          const DoubleTab& rho = eq_NS_PF.rho().valeurs();
          double rho0 = eq_NS_PF.rho0();
          const IntTab& face_voisins = le_dom->face_voisins();
          const DoubleVect& volumes = le_dom->volumes();

          for (num_cl=0 ; num_cl<le_dom->nb_front_Cl() ; num_cl++)
            {
              const Cond_lim& la_cl = le_dom_Cl->les_conditions_limites(num_cl);
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
              int ndeb = le_bord.num_premiere_face();
              int nfin = ndeb + le_bord.nb_faces();
              if (sub_type(Dirichlet,la_cl.valeur()) || sub_type(Dirichlet_homogene,la_cl.valeur()))
                ;
              else
                for (face=ndeb ; face<nfin ; face++)
                  {
                    int el0=face_voisins(face,0);
                    int el1=face_voisins(face,1);
                    double vol0=volumes(el0);
                    double vol1=volumes(el1);
                    double cbetaface=(vol0*(rho(el0)/rho0-1.0)+vol1*(rho(el1)/rho0-1.0))/(vol0+vol1);
                    resu(face) += cbetaface * g(orientation(face)) * porosite_surf(face) * volumes_entrelaces(face);
                  }
            }

          for (face=premiere_face_interne ; face<nb_faces; face++)
            {
              int el0=face_voisins(face,0);
              int el1=face_voisins(face,1);
              double vol0=volumes(el0);
              double vol1=volumes(el1);
              double cbetaface=(vol0*(rho(el0)/rho0-1.0)+vol1*(rho(el1)/rho0-1.0))/(vol0+vol1);
              resu(face) += cbetaface * g(orientation(face)) * porosite_surf(face) * volumes_entrelaces(face);
            }
        }
    }



  return resu;
}


/*! @brief Calcule la contribution de cette source
 *
 * @param (DoubleTab& resu) flux
 * @return (DoubleTab&) le flux
 */
DoubleTab& Source_Gravite_PF_VDF::calculer(DoubleTab& resu) const
{
  resu = 0.;
  return ajouter(resu);
}


