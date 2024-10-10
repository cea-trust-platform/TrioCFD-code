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
// File:        Modele_Jones_Launder_Thermique_VDF.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Modeles_Turbulence/RANS/Fonc
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_Jones_Launder_Thermique_VDF.h>
#include <Domaine_VDF.h>
#include <Domaine_Cl_VDF.h>
#include <TRUSTTrav.h>
#include <Probleme_base.h>
#include <Dirichlet.h>
#include <Dirichlet_entree_fluide_leaves.h>
#include <Symetrie.h>
#include <Neumann.h>
#include <Neumann_homogene.h>
#include <TRUSTTrav.h>

Implemente_instanciable(Modele_Jones_Launder_Thermique_VDF,"Modele_Jones_Launder_Thermique_VDF",Modele_Fonc_Bas_Reynolds_Thermique_Base);

// printOn et readOn

Sortie& Modele_Jones_Launder_Thermique_VDF::printOn(Sortie& s ) const
{
  return s;
}

Entree& Modele_Jones_Launder_Thermique_VDF::readOn(Entree& is )
{
  Motcle motlu, accolade_fermee="}", accolade_ouverte="{";
  is >> motlu;
  if (motlu==accolade_ouverte)
    {
      is >> motlu;
      if (motlu != accolade_fermee)
        {
          Cerr << "Erreur a la lecture du Modele fonc bas reynolds Jones et Launder pour la thermique" << finl;
          Cerr << "On attendait } a la place de " << motlu << finl;
          exit();
        }
    }
  else
    {
      Cerr << "Erreur a la lecture du Modele fonc bas reynolds Jones et Launder pour la thermique" << finl;
      Cerr << "On attendait { a la place de " << motlu << finl;
      exit();
    }
  return is;
}

Entree& Modele_Jones_Launder_Thermique_VDF::lire(const Motcle& , Entree& is)
{
  return is;
}
///////////////////////////////////////////////////////////////
//   Implementation des fonctions de la classe
///////////////////////////////////////////////////////////////

void  Modele_Jones_Launder_Thermique_VDF::associer(const Domaine_dis_base& domaine_dis,
                                                   const Domaine_Cl_dis_base& domaine_Cl_dis)
{
  //const Domaine_VDF& le_dom = ref_cast(Domaine_VDF,domaine_dis);
  //  const Domaine_Cl_VDF& le_dom_Cl = ref_cast(Domaine_Cl_VDF,domaine_Cl_dis);
}

void Modele_Jones_Launder_Thermique_VDF::associer_pb(const Probleme_base& pb )
{
  eq_transport_Fluctu_Temp_Bas_Re = ref_cast(Transport_Fluctuation_Temperature_W_Bas_Re,equation());
}

DoubleTab& Modele_Jones_Launder_Thermique_VDF::Calcul_D(DoubleTab& D,const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis,
                                                        const DoubleTab& vitesse,const DoubleTab& Fluctu_Temp, double diffu ) const
{
  const Domaine_VDF& le_dom = ref_cast(Domaine_VDF,domaine_dis);
  const Domaine_Cl_VDF& le_dom_Cl = ref_cast(Domaine_Cl_VDF,domaine_Cl_dis);
  D = 0;
  //  const DoubleVect& volumes = le_dom.volumes();
  const DoubleVect& porosite_surf = domaine_Cl_dis.equation().milieu().porosite_face();
  const DoubleVect& volume_entrelaces = le_dom.volumes_entrelaces();
  //  int nb_elem = le_dom.nb_elem();
  //  int nb_elem_tot = le_dom.nb_elem_tot();
  //const Domaine& domaine=le_dom.domaine();

  //int nb_faces_elem = domaine.nb_faces_elem();
  //IntTrav numfa(nb_faces_elem);
  double coef;
  //  const IntTab& elem_faces = le_dom.elem_faces();
  const IntTab& face_voisins = le_dom.face_voisins();
  int nb_faces = le_dom.nb_faces();

  DoubleTab gradth(nb_faces);
  int num_face,poly1,poly2,ori, ndeb, nfin;

  // Calcul de Gradient de racine de theta^2.

  // Boucle sur les bords pour traiter les conditions aux limites
  for (int n_bord=0; n_bord<le_dom.nb_front_Cl(); n_bord++)

    {
      const Cond_lim& la_cl = le_dom_Cl.les_conditions_limites(n_bord);

      if ( sub_type(Dirichlet,la_cl.valeur()) )
        {
          const Dirichlet_entree_fluide& la_cl_typee = ref_cast(Dirichlet_entree_fluide,la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          ndeb = le_bord.num_premiere_face();
          nfin = ndeb + le_bord.nb_faces();

          for (num_face=ndeb; num_face<nfin; num_face++)
            {
              gradth = 0;
              poly1 =  face_voisins(num_face,0);
              poly2=-1;
              // poiur faire planter si on en avait besoin
              if (poly1 != -1)
                {
                  coef = volume_entrelaces(num_face)*porosite_surf(num_face)*0.5;
                  if ( (Fluctu_Temp(poly1,0)>0) && (Fluctu_Temp(poly2,0)>0) )
                    gradth(num_face) += (coef*(la_cl_typee.val_imp(num_face-ndeb,0) - sqrt(Fluctu_Temp(poly1,0))))/le_dom.dist_norm_bord(num_face);
                  D[poly1] += 2*diffu*(gradth(num_face)*gradth(num_face));
                }
              else
                {
                  poly2 = face_voisins(num_face,1);
                  coef = volume_entrelaces(num_face)*porosite_surf(num_face)*0.5;
                  if ( (Fluctu_Temp(poly1,0)>0) && (Fluctu_Temp(poly2,0)>0) )
                    gradth(num_face) += (coef*(sqrt(Fluctu_Temp(poly2,0)) - la_cl_typee.val_imp(num_face-ndeb,0)))/le_dom.dist_norm_bord(num_face);
                  D[poly2] += 2*diffu*(gradth(num_face)*gradth(num_face));
                }
            }

        }
      else if (sub_type(Symetrie,la_cl.valeur()))
        ;
      else if ( (sub_type(Neumann,la_cl.valeur()))
                ||
                (sub_type(Neumann_homogene,la_cl.valeur()))
              )
        {
          // do nothing
          ;
        }
    }
  // Traitement des faces internes
  for (num_face=le_dom.premiere_face_int(); num_face<nb_faces; num_face++)
    {
      poly1 = face_voisins(num_face,0);
      poly2 = face_voisins(num_face,1);
      ori = le_dom.orientation(num_face);
      coef = volume_entrelaces(num_face)*porosite_surf(num_face);

      if ( (Fluctu_Temp(poly1,0)>0) && (Fluctu_Temp(poly2,0)>0) )
        gradth(num_face) += coef*(sqrt(Fluctu_Temp(poly1,0))-sqrt(Fluctu_Temp(poly2,0)))/(le_dom.xp(poly2,ori)- le_dom.xp(poly1,ori));

      D[poly1] += 2*diffu*(gradth(num_face)*gradth(num_face));
      D[poly2] += 2*diffu*(gradth(num_face)*gradth(num_face));
    }
  return D;
}

DoubleTab& Modele_Jones_Launder_Thermique_VDF::Calcul_E(DoubleTab& E,const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis, const DoubleTab& temp,const DoubleTab& Fluctu_Temp,double diffu, const DoubleTab& diffu_turb ) const
{
  const Domaine_VDF& le_dom = ref_cast(Domaine_VDF,domaine_dis);
  const Domaine_Cl_VDF& le_dom_Cl = ref_cast(Domaine_Cl_VDF,domaine_Cl_dis);
  E = 0;
  const DoubleVect& volumes = le_dom.volumes();
  //  const DoubleVect& porosite_vol = la_equation().milieu().porosite_elem();

  int nb_elem = le_dom.nb_elem();
  int nb_elem_tot = le_dom.nb_elem_tot();
  const IntTab& elem_faces = le_dom.elem_faces();
  int nb_faces = le_dom.nb_faces();
  const IntTab& face_voisins = le_dom.face_voisins();
  //  const IntTab& Qdm = le_dom.Qdm();
  //  const IntVect& orientation = le_dom.orientation();
  int ndeb,nfin,poly1, poly2, num_face, ori;
  const Domaine& domaine=le_dom.domaine();
  int nb_faces_elem = domaine.nb_faces_elem();

  DoubleTrav dT_dy(nb_faces);
  DoubleTrav dT_dz(nb_faces);
  DoubleTrav dT_dx(nb_faces);
  DoubleTrav d2T_dy2(nb_elem_tot);
  DoubleTrav d2T_dz2(nb_elem_tot);
  DoubleTrav d2T_dx2(nb_elem_tot);

  //Calcul des derives de T

  //Boucle sur les frontieres pour traiter les cds aux limites.
  for (int n_bord=0; n_bord<le_dom.nb_front_Cl(); n_bord++)       //boucle sur les frontieres
    {
      const Cond_lim& la_cl = le_dom_Cl.les_conditions_limites(n_bord);

      if ( sub_type(Dirichlet,la_cl.valeur()) )
        {
          const Dirichlet_entree_fluide& la_cl_typee = ref_cast(Dirichlet_entree_fluide,la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          ndeb = le_bord.num_premiere_face();
          nfin = ndeb + le_bord.nb_faces();

          for (num_face=ndeb; num_face<nfin; num_face++)      //boucle sur les faces des frontieres
            {
              //         dT_dx = 0;
              //         dT_dy = 0;
              //         dT_dz = 0;
              poly1 = face_voisins(num_face,0);
              ori = le_dom.orientation(num_face);
              if (poly1 != -1)
                {
                  if (ori == 0)
                    {
                      dT_dx(num_face) = (la_cl_typee.val_imp(num_face-ndeb,0) - temp(poly1))/le_dom.dist_norm_bord(num_face);
                    }
                  if (ori == 1)
                    {
                      dT_dy(num_face) = (la_cl_typee.val_imp(num_face-ndeb,0) - temp(poly1))/le_dom.dist_norm_bord(num_face);
                    }
                  if (ori == 2)
                    {
                      dT_dz(num_face) = (la_cl_typee.val_imp(num_face-ndeb,0) - temp(poly1))/le_dom.dist_norm_bord(num_face);
                    }
                }
              else
                {
                  poly2 = face_voisins(num_face,1);
                  if (ori == 0)
                    {
                      dT_dx(num_face) = (temp(poly2) - la_cl_typee.val_imp(num_face-ndeb,0))/le_dom.dist_norm_bord(num_face);
                    }
                  if (ori == 1)
                    {
                      dT_dy(num_face) = (temp(poly2) - la_cl_typee.val_imp(num_face-ndeb,0))/le_dom.dist_norm_bord(num_face);
                    }
                  if (ori == 2)
                    {
                      dT_dz(num_face) = (temp(poly2) - la_cl_typee.val_imp(num_face-ndeb,0))/le_dom.dist_norm_bord(num_face);
                    }
                }
            }
        }
    }


  //Traitement des faces internes toujours pour le calcul des derivees de T
  for (num_face=le_dom.premiere_face_int(); num_face<nb_faces; num_face++)      //boucle sur les faces internes
    {
      poly1 = face_voisins(num_face,0);
      poly2 = face_voisins(num_face,1);
      ori = le_dom.orientation(num_face);
      if (ori == 0)
        {
          dT_dx(num_face) = ( temp(poly2) - temp(poly1) ) / (le_dom.xp(poly2,ori) - le_dom.xp(poly1,ori));
        }
      if (ori == 1)
        {
          dT_dy(num_face) = ( temp(poly2) - temp(poly1) ) / (le_dom.xp(poly2,ori) - le_dom.xp(poly1,ori));
        }
      if (ori == 2)
        {
          dT_dz(num_face) = ( temp(poly2) - temp(poly1) ) / (le_dom.xp(poly2,ori) - le_dom.xp(poly1,ori));
        }
    }


  //on passe maintenant au calcul des derivees secondes de T et du terme E
  IntTrav numfa(nb_faces_elem);
  for (int elem=0; elem<nb_elem; elem++)      //boucle sur les elements du domaine
    {
      for (int i=0; i<nb_faces_elem; i++)
        numfa[i] = elem_faces(elem,i);
      if ( dimension == 2)
        {
          d2T_dx2(elem) = ( dT_dx(numfa[2]) - dT_dx(numfa[0]) ) * le_dom.dist_face(numfa[1],numfa[3],0);
          d2T_dy2(elem) = ( dT_dy(numfa[1]) - dT_dx(numfa[3]) ) * le_dom.dist_face(numfa[0],numfa[2],1);
          E[elem] = 2 * diffu * diffu_turb(elem) * ( d2T_dx2(elem) + d2T_dy2(elem) ) * ( d2T_dx2(elem) + d2T_dy2(elem) );
        }
      else if (dimension == 3)
        {
          d2T_dx2(elem) = ( dT_dx(numfa[3]) - dT_dx(numfa[0]) ) * volumes(numfa[3]);
          d2T_dy2(elem) = ( dT_dy(numfa[1]) - dT_dx(numfa[4]) ) * volumes(numfa[4]);
          d2T_dz2(elem) = ( dT_dy(numfa[5]) - dT_dx(numfa[2]) ) * volumes(numfa[5]);
          E[elem] = 2 * diffu * diffu_turb(elem) * (d2T_dx2(elem)+d2T_dy2(elem)+d2T_dz2(elem)) * (d2T_dx2(elem)+d2T_dy2(elem)+d2T_dz2(elem));
        }

    }

  return E;
}

DoubleTab& Modele_Jones_Launder_Thermique_VDF::Calcul_F1( DoubleTab& F1, const Domaine_dis_base& domaine_dis,const DoubleTab& K_Eps_Bas_Re,const DoubleTab& FluctuTemp_Bas_Re,double visco,double diffu) const
{
  const Domaine_VDF& le_dom = ref_cast(Domaine_VDF,domaine_dis);
  int nb_elem = le_dom.nb_elem();
  for (int elem=0; elem <nb_elem; elem ++ )
    F1[elem] = 1.;
  return F1;
}

DoubleTab& Modele_Jones_Launder_Thermique_VDF::Calcul_F2( DoubleTab& F2, const Domaine_dis_base& domaine_dis,const DoubleTab& K_Eps_Bas_Re,const DoubleTab& FluctuTemp_Bas_Re,double visco,double diffu) const
{
  const Domaine_VDF& le_dom = ref_cast(Domaine_VDF,domaine_dis);
  int nb_elem = le_dom.nb_elem();
  /* DoubleTab Re(nb_elem);
     int elem;

     for (elem=0; elem< nb_elem ; elem++)
     {
     Re(elem) = (K_eps_Bas_Re(elem,0)*K_eps_Bas_Re(elem,0))/(visco*K_eps_Bas_Re(elem,1));
     F2[elem] = 1. - (0.3*exp(-1*carre(Re(elem))));
     }*/
  for (int elem=0; elem <nb_elem; elem ++ )
    F2[elem] = 1.;
  return F2;
}

DoubleTab& Modele_Jones_Launder_Thermique_VDF::Calcul_F3( DoubleTab& F3, const Domaine_dis_base& domaine_dis,const DoubleTab& K_Eps_Bas_Re,const DoubleTab& FluctuTemp_Bas_Re,double visco,double diffu) const
{
  const Domaine_VDF& le_dom = ref_cast(Domaine_VDF,domaine_dis);
  int nb_elem = le_dom.nb_elem();
  for (int elem=0; elem <nb_elem; elem ++ )
    F3[elem] = 1.;
  return F3;
}
DoubleTab& Modele_Jones_Launder_Thermique_VDF::Calcul_F4( DoubleTab& F4, const Domaine_dis_base& domaine_dis,const DoubleTab& K_Eps_Bas_Re,const DoubleTab& FluctuTemp_Bas_Re,double visco,double diffu) const
{
  const Domaine_VDF& le_dom = ref_cast(Domaine_VDF,domaine_dis);
  int nb_elem = le_dom.nb_elem();
  for (int elem=0; elem <nb_elem; elem ++ )
    F4[elem] = 1.;
  return F4;
}


DoubleTab&  Modele_Jones_Launder_Thermique_VDF::Calcul_Flambda( DoubleTab& Flambda,const Domaine_dis_base& domaine_dis,const DoubleTab& K_Eps_Bas_Re, const DoubleTab& FluctuTemp_Bas_Re, double visco, double diffu) const
{
  const Domaine_VDF& le_dom = ref_cast(Domaine_VDF,domaine_dis);
  Flambda = 0;
  int nb_elem = le_dom.nb_elem();
  DoubleTab Rt(nb_elem);
  int elem;
  for (elem=0; elem< nb_elem ; elem++)
    {
      if ( (K_Eps_Bas_Re(elem,1)*FluctuTemp_Bas_Re(elem,1) >= 1.e-6) && (K_Eps_Bas_Re(elem,0)*FluctuTemp_Bas_Re(elem,0) >= 1.e-6))
        {
          Rt(elem) = K_Eps_Bas_Re(elem,0)/sqrt(visco*diffu)*sqrt( (K_Eps_Bas_Re(elem,0)*FluctuTemp_Bas_Re(elem,0))/2/(K_Eps_Bas_Re(elem,1)*FluctuTemp_Bas_Re(elem,1)) );
          //     Re(elem) = (K_Eps_Bas_Re(elem,0)*K_Eps_Bas_Re(elem,0))/(visco*K_Eps_Bas_Re(elem,1));
          Flambda[elem] = exp(-2.5/(1.+Rt(elem)/50.));
        }
      else
        Flambda[elem] = 1.;
    }
  return Flambda;
}

void  Modele_Jones_Launder_Thermique_VDF::mettre_a_jour(double temps)
{
  ;
}
