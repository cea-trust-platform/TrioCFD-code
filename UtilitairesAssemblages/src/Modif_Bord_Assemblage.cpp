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
// File:        Modif_Bord_Assemblage.cpp
// Directory:   $TRUST_ROOT/src/UtilitairesAssemblages
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////

#include <Modif_Bord_Assemblage.h>
#include <Domaine.h>
#include <Motcle.h>
#include <EFichier.h>
#include <UtilitairesPrepro.h>
#include <Coeur.h>
#include <Param.h>

Implemente_instanciable(Modif_Bord_Assemblage,"Modif_Bord_Assemblage",Interprete_geometrique_base);


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
Sortie& Modif_Bord_Assemblage::printOn(Sortie& os) const
{
  return Interprete::printOn(os);
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
Entree& Modif_Bord_Assemblage::readOn(Entree& is)
{
  return Interprete::readOn(is);
}


Entree& Modif_Bord_Assemblage::interpreter_(Entree& is)
{
  Nom nom_dom_loc;
  Param param(que_suis_je());
  param.ajouter("domaine",&nom_dom_loc,Param::REQUIRED);
  param.ajouter("nom_cl",&nom_cl,Param::REQUIRED);
  param.ajouter("M0",&M0,Param::REQUIRED);
  param.ajouter("entreplat",&entreplat,Param::REQUIRED);
  param.ajouter("epaisseur_jeu",&epaisseur_jeu,Param::REQUIRED);
  param.ajouter("nom_fichier",&nom_fichier,Param::REQUIRED);
  param.lire_avec_accolades_depuis(is);
  associer_domaine(nom_dom_loc);
  nb_couronnes = M0;
  differencier_faces_assemblages(domaine());
  return is;
}



void Modif_Bord_Assemblage::differencier_faces_assemblages(Domaine& dom)
{
  Zone& zone = dom.zone(0);

  EFichier fic(nom_fichier);

  if(fic.good())
    {

      int icl = zone.rang_frontiere(nom_cl);
      Frontiere& fr = zone.frontiere(icl);
      Faces& les_faces_bord=fr.faces();
      int nb_faces = les_faces_bord.nb_faces();

      IntVect type_face(nb_faces);
      type_face=0;


      // definition des parametres d'entree pour le calcul de correspondance
      //                   assemblage - faces de bord
      //////////////////////////////////////////////////////////////////////

      Coeur coeur;
      coeur.setEntreplat(entreplat);
      coeur.setEpaisseurJeu(epaisseur_jeu);
      coeur.setNbCouronnes(nb_couronnes); // dimensionne a M0 - a augmenter si necessaire
      coeur.setOrigineNum(M0);
      coeur.setBord(nom_cl);
      coeur.setDomaine(nom_dom);

      coeur.calculer_geom();


      // Lecture des assemblages a differencier au niveau des C-L
      //////////////////////////////////////////////////////////////////////

      Noms nom_type_assemblage;
      int type_assemblage,M_lu,N_lu;

      fic >> nom_type_assemblage;
      fic >> type_assemblage >> M_lu >> N_lu ;

      while(type_assemblage!=-1)
        {
          for(int i=1; i<7; i++)
            {
              int nb = coeur.getAssemblage_size(M_lu, N_lu, i);

              for (int j=0; j<nb; j++)
                type_face(coeur.getAssemblage_face(M_lu, N_lu, i, j)) = type_assemblage;

            }
          fic >> type_assemblage >> M_lu >> N_lu ;
        }


      // reprise de l'ancien bord asssemblage "allege"
      ////////////////////////////////////////////////

      IntTab som_old_bord(nb_faces,dimension);

      int k=0;

      for (int j=0; j<nb_faces; j++)
        {
          if (type_face(j)==0)
            {
              som_old_bord(k,0) = les_faces_bord.sommet(j,0);
              som_old_bord(k,1) = les_faces_bord.sommet(j,1);
              if(dimension==3) som_old_bord(k,2) = les_faces_bord.sommet(j,2);
              k++;
            }
        }


      // creation des nouveaux bords
      ///////////////////////////////

      int nsz = nom_type_assemblage.size();

      for (int i=0; i<nsz; i++)
        {
          Bord& new_bord = zone.faces_bord().add(Bord());
          new_bord.nommer(nom_type_assemblage[nsz-i-1]);     // pour contrer la logique VECT de lecture de Noms
          Cerr << i << " " << nom_type_assemblage[nsz-i-1] << finl;
          Faces& les_faces_new_bord=new_bord.faces();

          if(dimension==2) les_faces_new_bord.typer(segment_2D);
          else                   les_faces_new_bord.typer(triangle_3D);

          IntTab som_new_bord(nb_faces,dimension);

          int kk=0;

          for (int j=0; j<nb_faces; j++)
            {
              if (type_face(j)==i+1)
                {
                  som_new_bord(kk,0) = les_faces_bord.sommet(j,0);
                  som_new_bord(kk,1) = les_faces_bord.sommet(j,1);
                  if(dimension==3) som_new_bord(kk,2) = les_faces_bord.sommet(j,2);
                  kk++;
                }
            }

          som_new_bord.resize(kk,dimension);
          les_faces_new_bord.voisins().resize(kk,dimension);
          les_faces_new_bord.voisins()=-1;
          les_faces_new_bord.les_sommets().ref(som_new_bord);
        }

      som_old_bord.resize(k,dimension);
      les_faces_bord.voisins().resize(k,dimension);
      les_faces_bord.voisins()=-1;
      les_faces_bord.les_sommets().ref(som_old_bord);

    }//if(fic.good())
}
