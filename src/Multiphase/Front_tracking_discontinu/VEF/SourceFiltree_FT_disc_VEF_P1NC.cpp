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
// File:        SourceFiltree_FT_disc_VEF_P1NC.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src/VEF
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#include <SourceFiltree_FT_disc_VEF_P1NC.h>
#include <Probleme_FT_Disc_gen.h>
#include <Transport_Interfaces_FT_Disc.h>
#include <Parser.h>
#include <Domaine_VEF.h>
#include <Domaine_Cl_VEF.h>
#include <Neumann_sortie_libre.h>
#include <Symetrie.h>
#include <Periodique.h>
#include <Dirichlet.h>
#include <Dirichlet_homogene.h>

Implemente_instanciable(SourceFiltree_FT_disc_VEF_P1NC,"SourceFiltree_FT_disc_VEF_P1NC",Source_base);


Entree& SourceFiltree_FT_disc_VEF_P1NC::readOn(Entree& is)
{
  return lire(is);
}
Sortie& SourceFiltree_FT_disc_VEF_P1NC::printOn(Sortie& os) const
{
  SourceFiltree_FT_disc_base::ecrire_donnees(os);
  return os << que_suis_je();
}
/*! @brief Lit les parametres du terme source a partir d'un flot d'entree.
 *
 * @param (Entree& is) un flot d'entree
 */
Entree& SourceFiltree_FT_disc_VEF_P1NC::lire(Entree& is)
{
  return SourceFiltree_FT_disc_base::lire_donnees(is);
}

DoubleTab& SourceFiltree_FT_disc_VEF_P1NC::ajouter(DoubleTab& resu) const
{
  const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
  const Domaine_Cl_VEF& domaine_Cl_VEF = le_dom_Cl_VEF.valeur();
  const IntTab& face_voisins = domaine_VEF.face_voisins();
  const DoubleTab& xv = domaine_VEF.xv();
  const DoubleTab& Indicatrice = Indic_->valeur().valeurs();
  const DoubleVect& volumes_entrelaces = domaine_VEF.volumes_entrelaces();
  const DoubleVect& porosite_surf = equation().milieu().porosite_face();
  const int nb_front_Cl = domaine_VEF.nb_front_Cl();

  int n_bord,face, elem1,elem2,k, ndeb,nfin, offset;
  double indic, x,y,z = 0;
  // Boucle sur les conditions limites pour traiter les faces de bord
  for (n_bord=0 ; n_bord<nb_front_Cl ; n_bord++)
    {
      // pour chaque Condition Limite on regarde son type
      // Si face de Dirichlet on ne fait rien
      // Si face de Neumann, Periodique ou de Symetrie on calcule la contribution au terme source
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Neumann_sortie_libre,la_cl.valeur())||sub_type(Symetrie,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          ndeb = le_bord.num_premiere_face();
          nfin = ndeb + le_bord.nb_faces();
          for (face=ndeb; face<nfin; face++)
            {
              elem1 = face_voisins(face,0);
              indic = Indicatrice(elem1);
              x = xv(face,0);
              y = xv(face,1);
              if (dimension==3)
                {
                  z = xv(face,2);
                }
              if (indic<0.5)
                {
                  offset = 0;
                }
              else
                {
                  offset = dimension_;
                }
              for (k=0 ; k<dimension ; k++)
                {
                  fI_xyz_t[k]->setVar("x",x);
                  fI_xyz_t[k]->setVar("y",y);
                  fI_xyz_t[k]->setVar("z",z);
                  fI_xyz_t[k]->setVar("t",temps_);
                  resu(face,k) += fI_xyz_t[offset+k]->eval() * volumes_entrelaces(face) * porosite_surf(face);
                }
            }
        }
      else if (sub_type(Periodique,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          ndeb = le_bord.num_premiere_face();
          nfin = ndeb + le_bord.nb_faces();
          for (face=ndeb; face<nfin; face++)
            {
              elem1 = face_voisins(face,0);
              elem2 = face_voisins(face,1);
              indic = (Indicatrice(elem1) + Indicatrice(elem2)) /2.;
              x = xv(face,0);
              y = xv(face,1);
              if (dimension==3)
                {
                  z = xv(face,2);
                }
              if (indic<0.5)
                {
                  offset = 0;
                }
              else
                {
                  offset = dimension_;
                }
              for (k=0 ; k<dimension ; k++)
                {
                  fI_xyz_t[k]->setVar("x",x);
                  fI_xyz_t[k]->setVar("y",y);
                  fI_xyz_t[k]->setVar("z",z);
                  fI_xyz_t[k]->setVar("t",temps_);
                  resu(face,k) += fI_xyz_t[offset+k]->eval() * volumes_entrelaces(face) * porosite_surf(face);
                }
            }
        }
      else if ( (sub_type(Dirichlet,la_cl.valeur())) || (sub_type(Dirichlet_homogene,la_cl.valeur())) )
        {
          ;
        }
    }
  //faces internes
  for (face =domaine_VEF.premiere_face_int(); face<domaine_VEF.nb_faces(); face++)
    {
      elem1 = face_voisins(face,0);
      elem2 = face_voisins(face,1);
      indic = (Indicatrice(elem1) + Indicatrice(elem2)) /2.;
      x = xv(face,0);
      y = xv(face,1);
      if (dimension==3)
        {
          z = xv(face,2);
        }
      if (indic<0.5)
        {
          offset = 0;
        }
      else
        {
          offset = dimension_;
        }
      for (k=0 ; k<dimension ; k++)
        {
          fI_xyz_t[k]->setVar("x",x);
          fI_xyz_t[k]->setVar("y",y);
          fI_xyz_t[k]->setVar("z",z);
          fI_xyz_t[k]->setVar("t",temps_);
          resu(face,k) += fI_xyz_t[offset+k]->eval() * volumes_entrelaces(face) * porosite_surf(face);
        }
    }

  return resu;
}
DoubleTab& SourceFiltree_FT_disc_VEF_P1NC::calculer(DoubleTab& resu) const
{
  resu = 0.;
  return ajouter(resu);
}
void SourceFiltree_FT_disc_VEF_P1NC::mettre_a_jour(double temps)
{
  temps_ = temps;
}
void SourceFiltree_FT_disc_VEF_P1NC::completer()
{
  Source_base::completer();
}

void SourceFiltree_FT_disc_VEF_P1NC::associer_domaines(const Domaine_dis_base& domaine_dis,const Domaine_Cl_dis_base& domaine_Cl_dis)
{
  le_dom_VEF = ref_cast(Domaine_VEF, domaine_dis);
  le_dom_Cl_VEF = ref_cast(Domaine_Cl_VEF, domaine_Cl_dis);
}
void SourceFiltree_FT_disc_VEF_P1NC::associer_pb(const Probleme_base& pb)
{
  if (sub_type(Probleme_FT_Disc_gen,pb))
    {
      const Probleme_FT_Disc_gen& pbFT = ref_cast(Probleme_FT_Disc_gen,pb);
      const Transport_Interfaces_FT_Disc& eq_transport = pbFT.equation_interfaces(nom_eq_transport_);
      Indic_ = eq_transport.inconnue();
    }
  else
    {
      Cerr<<"Un probleme "<<pb.que_suis_je()<<" ne peut etre associe a une source SourceFiltree_FT_disc_VEF_P1NC"<<finl;
      Cerr<<"Il faut utiliser un pb Pb_Front_Tracking_base"<<finl;
    }
}
