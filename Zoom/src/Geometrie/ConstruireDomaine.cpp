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
// File:        ConstruireDomaine.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Geometrie
// Version:     /main/14
//
//////////////////////////////////////////////////////////////////////////////

#include <ConstruireDomaine.h>
#include <Domaine.h>
#include <Sous_Zone.h>
#include <Deriv_Zone.h>
#include <Static_Int_Lists.h>
#include <Connectivite_som_elem.h>
#include <Faces_builder.h>

Implemente_instanciable(ConstruireDomaine,"ConstruireDomaine",Interprete_geometrique_base);


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
Sortie& ConstruireDomaine::printOn(Sortie& os) const
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
Entree& ConstruireDomaine::readOn(Entree& is)
{
  return Interprete::readOn(is);
}


// Description:
//    Fonction principale de l'interprete Raffiner
//    Triangule 1 a 1 toutes les zones du domaine
//    specifie par la directive.
//    On triangule la zone grace a la methode:
//      void Raffiner::raffiner(Zone& zone) const
//    Raffiner signifie ici transformer en triangle des
//    elements geometrique d'une zone.
// Precondition: on doit etre en 2D (dimension d'espace=2)
// Parametre: Entree& is
//    Signification: un flot d'entree
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: Entree&
//    Signification: le flot d'entree
//    Contraintes:
// Exception: l'objet a mailler n'est pas du type Domaine
// Effets de bord:
// Postcondition:
Entree& ConstruireDomaine::interpreter_(Entree& is)
{

  Nom nom_dom, nom_sszone;
  is >> nom_dom >> nom_sszone;
  associer_domaine(nom_dom);
  Domaine& dom=domaine();
  Sous_Zone& sszone=ref_cast(Sous_Zone, objet(nom_sszone));





  DERIV(Zone) zone_raf;
  zone_raf.typer("Zone");
  Zone& zoneraf = dom.add(zone_raf.valeur());
  zoneraf.associer_domaine(dom);



  Cout << "Type elem = " << sszone.zone().type_elem()->que_suis_je() << finl;

  zoneraf.typer("Rectangle");

  IntTab correspo_new_som;
  creer_sommet_et_elem(dom,sszone,correspo_new_som);

  creer_bords(zoneraf,sszone,correspo_new_som);



  //Domaine& dom0 = zone0.domaine();
  // Cerr << "Zone 0 nb faces bord = " << zone0.nb_faces_bord() << finl;


  return is;
}



void ConstruireDomaine::creer_sommet_et_elem(Domaine& dom, Sous_Zone& ssz,IntTab& correspo_newsom)
{
  Zone& zone0 = ssz.zone();
  Domaine& dom0 = zone0.domaine();

  const int nb_poly = ssz.nb_elem_tot();
  const int nb_som_poly = zone0.nb_som_elem();
  int i,j;
  const DoubleTab& les_coord = dom0.coord_sommets();
  DoubleTab newsom(nb_som_poly*nb_poly,Objet_U::dimension);
  int compteur = 0;
  IntTab new_elems(nb_poly, nb_som_poly); // les nouveaux elements
  correspo_newsom.resize(zone0.nb_som_tot());
  correspo_newsom = -1;
  IntTab& les_elems_= dom.zone(0).les_elems();
  les_elems_.ref(new_elems);

  for (i=0; i<nb_poly; i++)
    {
      IntTab indice(nb_som_poly);
      //Cout << "poly = " << i << " " << ssz(i) << finl;
      for (j=0; j<nb_som_poly; j++)
        {
          DoubleTab un_som(1, 3);
          //Cout << "som = " << j << " " << zone0.sommet_elem(i,j) << finl;

          un_som(0,0) = les_coord(zone0.sommet_elem(ssz(i),j),0);
          un_som(0,1) = les_coord(zone0.sommet_elem(ssz(i),j),1);
          //Cout << "coord   = " << un_som(0,0) << " " << un_som(0,1) << finl;

          if (dimension==3)  un_som(0,2) = les_coord(zone0.sommet_elem(ssz(i),j),2);
          indice(j) = ajouterSom(dom,un_som,newsom,compteur);
          new_elems(i,j) = indice(j);
          correspo_newsom(zone0.sommet_elem(ssz(i),j)) = indice(j);
        }
    }

  Cerr << "Correspondance  = " << correspo_newsom << finl;

}

int ConstruireDomaine::ajouterSom(Domaine& dom, const DoubleTab& un_som, DoubleTab& newsom, int& compteur)
{
  const DoubleTab& les_coord = dom.coord_sommets();
  const int nbsom = dom.nb_som();
  int i;
  int trouve = -1;
  //Cout << "nb som = " << nbsom << finl;


  if (dimension == 3)
    {
      for (i=0; i<nbsom; i++)
        {
          if ( (les_coord(i,0) == un_som(0,0))  && (les_coord(i,1) == un_som(0,1)) && (les_coord(i,2) == un_som(0,2)))
            {
              trouve = i;
              break;
            }
        }
      if (trouve == -1)
        {
          newsom(compteur,0) = un_som(0,0);
          newsom(compteur,1) = un_som(0,1);
          newsom(compteur,2) = un_som(0,2);
          trouve = compteur;
          compteur++;
          dom.ajouter(un_som);
        }
    }
  else
    {
      for (i=0; i<nbsom; i++)
        {
          if ( (les_coord(i,0) == un_som(0,0)) && (les_coord(i,1) == un_som(0,1)))
            {
              trouve = i;
              break;
            }
        }
      if (trouve == -1)
        {
          newsom(compteur,0) = un_som(0,0);
          newsom(compteur,1) = un_som(0,1);
          trouve = compteur;
          compteur++;
          dom.ajouter(un_som);
        }
    }
  return trouve;
}


void ConstruireDomaine::creer_bords(Zone& zoneraf, Sous_Zone& ssz, IntTab& correspo_ns)
{
  Zone& zone0 = ssz.zone();


  Bords& mes_faces_bord = zoneraf.faces_bord();
  Elem_geom& elem = zoneraf.type_elem();
  LIST_CURSEUR(Bord) curseur(zone0.faces_bord());
  while(curseur)
    {
      Faces& les_faces=curseur->faces();
      Bord& bord=mes_faces_bord.add(Bord());
      Cout << "Dans ConstruireDomaine  Nom bord = " << curseur->le_nom() << finl;
      //Cout << "Dans ConstruireDomaine  type face  = " << les_faces.type_face() << finl;
      bord.nommer(curseur->le_nom());
      bord.typer_faces(les_faces.type_face());
      bord.associer_zone(zoneraf);

      int nb_faces = les_faces.nb_faces_tot();
      int nb_som_faces = les_faces.nb_som_faces();
      IntTab& les_sommets = les_faces.les_sommets();
      IntTab uneface(1,nb_som_faces);

      Cout << "Dans ConstruireDomaine  nb face & som/face = " << nb_faces  << " " << nb_som_faces << finl;
      int i;
      bord.dimensionner(0);
      for (i=0; i<nb_faces; i++)
        {
          Cout << "sommet (i) = " << i << " " << les_sommets(i,0) << finl;
          Cout << "sommet (i) = " << i << " " << les_sommets(i,1) << finl;
          if (correspo_ns(les_sommets(i,0)) !=-1 && correspo_ns(les_sommets(i,1)) !=-1)
            {
              //        Cout << "   new sommet (i) = " << i << " " << correspo_ns(les_sommets(i,0)) << finl;
              //        Cout << "   new sommet (i) = " << i << " " << correspo_ns(les_sommets(i,1)) << finl;
              // la face doit etre prise en compte pour le nouveau domaine
              uneface(0,0) = correspo_ns(les_sommets(i,0));
              uneface(0,1) = correspo_ns(les_sommets(i,1));

              bord.ajouter_faces(uneface);
            }
        }

      //Cout << " faces = " << les_faces << finl;
      ++curseur;
    }

  Bord& bord=mes_faces_bord.add(Bord());
  bord.nommer("interface");
  bord.typer_faces(elem.type_face());
  bord.associer_zone(zoneraf);
  bord.dimensionner(0);

  Faces mes_faces;
  {
    // bloc a factoriser avec Zone_VF.cpp :
    Type_Face type_face = zone0.type_elem().type_face(0);
    mes_faces.typer(type_face);
    mes_faces.associer_zone(zone0);

    Static_Int_Lists connectivite_som_elem;
    const int     nb_sommets_tot = zone0.domaine().nb_som_tot();
    const IntTab&    elements       = zone0.les_elems();

    construire_connectivite_som_elem(nb_sommets_tot,
                                     elements,
                                     connectivite_som_elem,
                                     1 /* include virtual elements */);

    Faces_builder faces_builder;
    IntTab elem_faces; // Tableau dont on aura pas besoin
    faces_builder.creer_faces_reeles(zone0,
                                     connectivite_som_elem,
                                     mes_faces,
                                     elem_faces);
  }

  Cerr << "Dans ConstruireDomaine  La zone construite nb elem ; nb som = " << zoneraf.nb_elem() << " " << zoneraf.nb_som() << finl;
  Cerr << "Dans ConstruireDomaine  La zone construite nb bords = " << zoneraf.nb_faces_bord() << finl;
  Cerr << "Dans ConstruireDomaine  La zone construite  bords = " << zoneraf.faces_bord() << finl;
  //Cerr << "Dans ConstruireDomaine  La zone0  bords = " << zone0.faces_bord() << finl;

}
