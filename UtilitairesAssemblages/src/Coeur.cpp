/****************************************************************************
* Copyright (c) 2015, CEA
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
// File:        Coeur.cpp
// Directory:   $TRUST_ROOT/src/UtilitairesAssemblages
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////

#include <Coeur.h>
#include <Interprete.h>
#include <Probleme_base.h>
#include <Equation_base.h>
#include <Motcle.h>
#include <Static_Int_Lists.h>
#include <Param.h>
#include <Zone_VEF.h>
#include <Zone_Cl_VEF.h>
#include <UtilitairesPrepro.h>
#include <Domaine.h>
#include <SFichier.h>

Implemente_instanciable(Coeur,"Coeur",Objet_U);


// Description:
//    Simple appel a: Interprete::printOn(Sortie&)
// Precondition:
// Parametre: Sortie& os
//    Signification: un flot de sortie
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: Sortie&
//    Signification: le flot de sortie modifie
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
Sortie& Coeur::printOn(Sortie& os) const
{
  return printOn(os);
}


// Description:
//    Simple appel a: Interprete::readOn(Entree&)
// Precondition:
// Parametre: Entree& is
//    Signification: un flot d'entree
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: Entree&
//    Signification: le flot d'entree modifie
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
Entree& Coeur::readOn(Entree& is)
{
  test=0;
  Param param(que_suis_je());
  param.ajouter("probleme",&nom_pb,Param::REQUIRED);
  param.ajouter("type_probleme",&jeux_ou_assemblages,Param::REQUIRED);
  param.ajouter("nom_bord",&nom_bord,Param::REQUIRED);
  param.ajouter("entreplat",&ep,Param::REQUIRED);
  param.ajouter("epaisseur_jeu",&ep_jeux,Param::REQUIRED);
  param.ajouter("nb_couronnes",&nb_couronnes,Param::REQUIRED);
  param.ajouter("origine_numerotation",&M0,Param::REQUIRED);
  param.ajouter("test",&test);
  param.lire_avec_accolades_depuis(is);


  calculer();
  //calculer_geom_2D();

  return is;
}

// Description:
// Precondition:
// Parametre: Zone& zone
//    Signification: la zone dont on veut raffiner les elements
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour:
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
void Coeur::calculer()
{
  // La connectivite est ici traitee sur le probleme discretise
  // suivant la numerotation GLOBALE des faces de bord pour la frontiere consideree

  // A voir + tard si on conserve cette methode (sert-elle ?)


  Probleme_base& pb = ref_cast(Probleme_base, interprete().objet(nom_pb));

  Zone_dis_base& zdis = pb.equation(0).zone_dis().valeur();

  if (!sub_type(Zone_VEF, zdis))
    {
      Cerr << "L'interprete CalculerConnectFacesHexaRNR ne marche qu'en VEF-3D pour une geometrie de type coeur RNR a assemblages hexaedriques" << finl;
      Process::exit();
    }

  Zone_VEF& zvef = ref_cast(Zone_VEF, zdis);
  Zone_Cl_VEF& zclvef = ref_cast(Zone_Cl_VEF, pb.equation(0).zone_Cl_dis().valeur());
  DoubleTab& xv = zvef.xv();

  Conds_lim& les_cl = zclvef.les_conditions_limites();
  int nb_front=zvef.nb_front_Cl();

  int icl=-1;
  for(int i=0; i< nb_front; i++)
    {
      Frontiere& fr = les_cl[i].valeur().frontiere_dis().frontiere();
      if(fr.le_nom() == nom_bord)
        icl = i;
    }

  if (icl<0)
    {
      Cerr << nom_bord << " n'est pas une CL du probleme " <<pb << finl;
      Process::exit();
    }

  const Cond_lim& la_cl = zclvef.les_conditions_limites(icl);
  const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
  int nb_faces_tot=le_bord.nb_faces_tot();

  IntVect type_face(nb_faces_tot); // vaut 1, 2, ...6 selon la position relative de la face sur son assemblage d'appartenance
  DoubleVect face_hexa_x(nb_faces_tot); // fournit pour chaque face du bord  la coordonnee x du CG de l'assemblage la contenant
  DoubleVect face_hexa_y(nb_faces_tot); // fournit pour chaque face du bord  la coordonnee y du CG de l'assemblage la contenant
  IntVect face_hexa_M(nb_faces_tot);
  IntVect face_hexa_N(nb_faces_tot);

  int nb_assemblages=(2*nb_couronnes+1)*(2*nb_couronnes+1); // sur-dimensionne pour l'instant
  vect_assemblages.dimensionner(nb_assemblages);

  pas_y = ep+ep_jeux;
  pas_x = pas_y*sqrt(3.);

  double ep_test_int = ep + 0.5*ep_jeux;

  for (int iface=0; iface<nb_faces_tot; iface++)
    {
      int num_face = le_bord.num_face(iface);

      double x = xv(num_face,0);
      double y = xv(num_face,1);

      double sgn_x = 1.;
      if(x<0.) sgn_x = -1.;

      double sgn_y = 1.;
      if(y<0.) sgn_y = -1.;


      int i = int((sgn_x*x)/(pas_x));
      int j = int((sgn_y*y)/(pas_y));

      double xj = sgn_x*(2*i+1)*0.5*pas_x;
      double yj = sgn_y*(2*j+1)*0.5*pas_y;

      int k = int((sgn_x*x+0.5*pas_x)/(pas_x));
      int l = int((sgn_y*y+0.5*pas_y)/(pas_y));

      double xi = sgn_x*k*pas_x;
      double yi = sgn_y*l*pas_y;


      if  ( UtilitairesPrepro::hexa(x-xi, y-yi, ep_test_int)==0 )
        {
          face_hexa_x(iface) = xi ;
          face_hexa_y(iface) = yi ;
          face_hexa_M(iface) = UtilitairesPrepro::coord_assemblage_M(xi, yi, M0, pas_y) ;
          face_hexa_N(iface) = UtilitairesPrepro::coord_assemblage_N(xi, yi, M0, pas_y) ;

          int iass = (face_hexa_M(iface)-M0+nb_couronnes)*(2*nb_couronnes+1)+(face_hexa_N(iface)-M0+nb_couronnes);

          type_face(iface) = UtilitairesPrepro::hexa(xv(num_face,0)-xi, xv(num_face,1)-yi, 0.9*ep);
          assert(type_face(iface)!=0);

          vect_assemblages[iass].addFace(type_face(iface), num_face);
        }
      if  ( UtilitairesPrepro::hexa(x-xj, y-yj, ep_test_int)==0 )
        {
          face_hexa_x(iface) = xj ;
          face_hexa_y(iface) = yj ;
          face_hexa_M(iface) = UtilitairesPrepro::coord_assemblage_M(xj, yj, M0, pas_y) ;
          face_hexa_N(iface) = UtilitairesPrepro::coord_assemblage_N(xj, yj, M0, pas_y) ;

          int iass = (face_hexa_M(iface)-M0+nb_couronnes)*(2*nb_couronnes+1)+(face_hexa_N(iface)-M0+nb_couronnes);

          type_face(iface) = UtilitairesPrepro::hexa(xv(num_face,0)-xj, xv(num_face,1)-yj, 0.9*ep);
          assert(type_face(iface)!=0);

          vect_assemblages[iass].addFace(type_face(iface), num_face);
        }

    }

  for (int iass=0; iass<nb_assemblages; iass++)
    vect_assemblages[iass].reordonner_faces(xv);


  if (test)
    {
      SFichier fic("assemblages.dat");
      //for (int iass=0;iass<nb_assemblages;iass++)
      int iass=12;
      {
        fic << "Assemblage " << iass << finl;
        int nb = vect_assemblages[iass].getFaces(1).size();
        fic << "nb_faces " << nb << finl;
        assert (nb==vect_assemblages[iass].getFaces(2).size());
        assert (nb==vect_assemblages[iass].getFaces(3).size());
        assert (nb==vect_assemblages[iass].getFaces(4).size());
        assert (nb==vect_assemblages[iass].getFaces(5).size());
        assert (nb==vect_assemblages[iass].getFaces(6).size());
        for (int j=0; j<nb; j++)
          {
            fic << vect_assemblages[iass].getFaces(1)[j] << " "
                << vect_assemblages[iass].getFaces(2)[j]  << " "
                << vect_assemblages[iass].getFaces(3)[j]  << " "
                << vect_assemblages[iass].getFaces(4)[j]  << " "
                << vect_assemblages[iass].getFaces(5)[j]  << " "
                << vect_assemblages[iass].getFaces(6)[j]  << " "
                << finl;
          }
      }
    }
}



void Coeur::calculer_geom()
{
  // La connectivite est ici traitee sur le probleme NON discretise
  // suivant la numerotation LOCALE des faces de a frontiere consideree
  // (pour agir directement au niveau du fichier .geom)


  const Domaine& dom =  ref_cast(Domaine, interprete().objet(nom_dom));
  const Zone& zone = dom.zone(0);
  const DoubleTab& coord_sommets = dom.coord_sommets();

  int icl = zone.rang_frontiere(nom_bord);
  const Frontiere& fr = zone.frontiere(icl);
  const Faces& les_faces_bord=fr.faces();
  int nb_faces_bord = les_faces_bord.nb_faces();

  int nb_assemblages=(2*nb_couronnes+1)*(2*nb_couronnes+1); // sur-dimensionne pour l'instant
  vect_assemblages.dimensionner(nb_assemblages);

  pas_y = ep+ep_jeux;
  pas_x = pas_y*sqrt(3.);

  double x,y;

  double ep_test_int = ep + 0.5*ep_jeux ;
  for (int iface=0; iface<nb_faces_bord; iface++)
    {

      int i0=les_faces_bord.sommet(iface,0);
      int i1=les_faces_bord.sommet(iface,1);

      if(dimension==2)
        {
          x = (coord_sommets(i0,0)+coord_sommets(i1,0))/2.;
          y = (coord_sommets(i0,1)+coord_sommets(i1,1))/2.;
        }
      else
        {
          int i2=les_faces_bord.sommet(iface,2);

          x = (coord_sommets(i0,0)+coord_sommets(i1,0)+coord_sommets(i2,0))/3.;
          y = (coord_sommets(i0,1)+coord_sommets(i1,1)+coord_sommets(i2,1))/3.;
        }


      double sgn_x = 1.;
      if(x<0.) sgn_x = -1.;

      double sgn_y = 1.;
      if(y<0.) sgn_y = -1.;

      int i = int((sgn_x*x)/(pas_x));
      int j = int((sgn_y*y)/(pas_y));

      double xj = sgn_x*(2*i+1)*0.5*pas_x;
      double yj = sgn_y*(2*j+1)*0.5*pas_y;

      int k = int((sgn_x*x+0.5*pas_x)/(pas_x));
      int l = int((sgn_y*y+0.5*pas_y)/(pas_y));

      double xi = sgn_x*k*pas_x;
      double yi = sgn_y*l*pas_y;


      if  ( UtilitairesPrepro::hexa(x-xi, y-yi, ep_test_int)==0 )
        {
          int face_hexa_M = UtilitairesPrepro::coord_assemblage_M(xi, yi, M0, pas_y) ;
          int face_hexa_N = UtilitairesPrepro::coord_assemblage_N(xi, yi, M0, pas_y) ;
          int type_face   = UtilitairesPrepro::hexa(x-xi, y-yi, 0.9*ep);

          int iass = (face_hexa_M-M0+nb_couronnes)*(2*nb_couronnes+1)+(face_hexa_N-M0+nb_couronnes);

          vect_assemblages[iass].addFace(type_face, iface);
        }
      if  ( UtilitairesPrepro::hexa(x-xj, y-yj, ep_test_int)==0 )
        {
          int face_hexa_M = UtilitairesPrepro::coord_assemblage_M(xj, yj, M0, pas_y) ;
          int face_hexa_N = UtilitairesPrepro::coord_assemblage_N(xj, yj, M0, pas_y) ;
          int type_face   = UtilitairesPrepro::hexa(x-xj, y-yj, 0.9*ep);

          int iass = (face_hexa_M-M0+nb_couronnes)*(2*nb_couronnes+1)+(face_hexa_N-M0+nb_couronnes);

          vect_assemblages[iass].addFace(type_face, iface);
        }

    }
}
