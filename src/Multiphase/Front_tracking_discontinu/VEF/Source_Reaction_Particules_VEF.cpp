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
// File:        Source_Reaction_Particules_VEF.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src/VEF
// Version:     /main/12
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Reaction_Particules_VEF.h>
#include <TRUSTTab.h>
#include <Transport_Marqueur_FT.h>
#include <Domaine.h>
#include <Domaine_VF.h>
#include <Champ_implementation_P1.h>
#include <Dirichlet.h>
#include <Dirichlet_homogene.h>
#include <Milieu_base.h>
#include <Scatter.h>
#include <Echange_EV_Options.h>
#include <MD_Vector_tools.h>
#include <Probleme_base.h>

Implemente_instanciable(Source_Reaction_Particules_VEF,"Two_Way_Coupling_VEF",Source_Reaction_Particules);


Sortie& Source_Reaction_Particules_VEF::printOn(Sortie& os) const
{
  return os;
}

Entree& Source_Reaction_Particules_VEF::readOn(Entree& is)
{
  return is;
}

//-Interpolation de source_stockage_ aux sommets (par ponderation par coordonnes barycentriques)
//-Dans le cas parallele une sommation aux sommets items communs est faite
// pour tenir compte des contributions de chaque processeur
//-Calcul de la source aux faces en sommant la condribution de chacun des sommets de la face
//Rq : Le terme source est divisee par rho si le milieu n est pas quasi-compressible
DoubleTab& Source_Reaction_Particules_VEF::ajouter(DoubleTab& resu) const
{
  const DoubleTab& source_action = eq_marqueur_->source_stockage();
  const Domaine_dis_base& zdis = equation().domaine_dis();
  const Domaine_VF& zvf = ref_cast(Domaine_VF,zdis);
  const Domaine_Cl_dis_base& zcldis = equation().domaine_Cl_dis();
  const Domaine& domaine_geom = zdis.domaine();
  const DoubleTab& coord = domaine_geom.coord_sommets();
  const IntTab& sommet_poly = domaine_geom.les_elems();
  const IntTab& face_som = zvf.face_sommets();
  const Maillage_FT_Disc& maillage = eq_marqueur_->maillage_interface();
  const DoubleTab& positions = maillage.sommets();
  const ArrOfInt& som_elem = maillage.sommet_elem();
  int dim = Objet_U::dimension;
  int nb_faces = zvf.nb_faces();
  int premiere_face_interne = zvf.premiere_face_int();
  int nb_pos = positions.dimension(0);
  int nb_som_elem = domaine_geom.nb_som_elem();
  int nb_som_face = domaine_geom.type_elem()->nb_som_face();
  int is_QC=0;
  int is_FT=0;
  double coord_b=0.;
  double x_pos,y_pos;
  double z_pos=0.;
  int elem,som_glob;
  const int nb_elem_reels = domaine_geom.nb_elem();


  //creation d un espace virtuel pour source_som
  DoubleTab source_som(0, dim);
  domaine_geom.creer_tableau_sommets(source_som);

  //is_QC=0 le terme source est homogne a dv/dt sinon homogene a d(rho*v)/dt
  is_QC = equation().probleme().is_dilatable();
  if (sub_type(Navier_Stokes_FT_Disc,equation()))
    is_FT = 1;

  const Champ_base& champ_rho = (is_FT==1?ref_cast(Navier_Stokes_FT_Disc,equation()).champ_rho_faces():equation().milieu().masse_volumique().valeur());
  const DoubleTab& rho_faces = champ_rho.valeurs();

  //Remplissage de source_som (interpolation de source_stockage aux sommets)

  //verifier que 0<=coord_b<=1
  for (int pos=0; pos<nb_pos; pos++)
    {
      x_pos = positions(pos,0);
      y_pos = positions(pos,1);
      if (dim==3)
        z_pos = positions(pos,2);

      elem = som_elem[pos];
      // On va sommer les contributions de tous les sommets communs,
      // ne pas traiter les elements virtuels, on aurait plusieurs fois la contribution du meme point
      if (elem < nb_elem_reels)
        {
          for (int som=0; som<nb_som_elem; som++)
            {
              som_glob = sommet_poly(elem,som);
              if (dim==2)
                coord_b = coord_barycentrique_P1_triangle(sommet_poly,coord,x_pos,y_pos,elem,som) ;
              else if (dim==3)
                coord_b = coord_barycentrique_P1_tetraedre(sommet_poly,coord,x_pos,y_pos,z_pos,elem,som) ;

              for (int k=0; k<dim; k++)
                source_som(som_glob,k) -= source_action(pos,k)*coord_b;
            }
        }
    }

  // Sommation des contributions de tous les sommets partages entre plusieurs processeurs,
  // puis echange espace virtuel.
  MD_Vector_tools::echange_espace_virtuel(source_som, MD_Vector_tools::EV_SOMME_ECHANGE);

  //Remplissage de resu (Interpolation de source_som aux faces)

  //Traitement des conditions limites
  int num_cl;
  for (num_cl=0 ; num_cl<zdis.nb_front_Cl() ; num_cl++)
    {
      const Cond_lim& la_cl = zcldis.les_conditions_limites(num_cl);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
      int ndeb = le_bord.num_premiere_face();
      int nfin = ndeb + le_bord.nb_faces();
      if (sub_type(Dirichlet,la_cl.valeur()) || sub_type(Dirichlet_homogene,la_cl.valeur()))
        ;
      else
        for (int face=ndeb ; face<nfin ; face++)
          {
            double rho = 1.;
            if (!is_QC)
              rho = (is_FT==1?rho_faces(face):rho_faces(0,0));
            for (int nbsf=0; nbsf<nb_som_face; nbsf++)
              {
                som_glob = face_som(face,nbsf);
                for (int k=0; k<dim; k++)
                  resu(face,k) += (1./nb_som_face)*source_som(som_glob,k)/rho;
              }
          }
    }

  //Traitement des faces internes
  for (int face=premiere_face_interne; face<nb_faces; face++)
    {
      double rho = 1.;
      if (!is_QC)
        rho = (is_FT==1?rho_faces(face):rho_faces(0,0));

      for (int nbsf=0; nbsf<nb_som_face; nbsf++)
        {
          som_glob = face_som(face,nbsf);
          for (int k=0; k<dim; k++)
            resu(face,k) += (1./nb_som_face)*source_som(som_glob,k)/rho;
        }
    }

  return resu;
}


