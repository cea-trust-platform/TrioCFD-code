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

#include <Connectivite_som_elem.h>
#include <ConstruireDomaine.h>
#include <Static_Int_Lists.h>
#include <Faces_builder.h>
#include <TRUST_Deriv.h>
#include <Sous_Domaine.h>
#include <Domaine.h>

Implemente_instanciable(ConstruireDomaine,"ConstruireDomaine",Interprete_geometrique_base);

Sortie& ConstruireDomaine::printOn(Sortie& os) const
{
  return Interprete::printOn(os);
}

Entree& ConstruireDomaine::readOn(Entree& is)
{
  return Interprete::readOn(is);
}


/*! @brief Fonction principale de l'interprete Raffiner Triangule 1 a 1 tous les domaines
 *
 *     specifie par la directive.
 *     On triangule le domaine grace a la methode:
 *       void Raffiner::raffiner(Domaine& domaine) const
 *     Raffiner signifie ici transformer en triangle des
 *     elements geometrique d'un domaine.
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree
 * @throws l'objet a mailler n'est pas du type Domaine
 */
Entree& ConstruireDomaine::interpreter_(Entree& is)
{
  Nom nom_dom, nom_ssdomaine;
  is >> nom_dom >> nom_ssdomaine;
  associer_domaine(nom_dom);
  Domaine& dom=domaine();
  Sous_Domaine& ssdomaine=ref_cast(Sous_Domaine, objet(nom_ssdomaine));

  Cout << "Type elem = " << ssdomaine.domaine().type_elem()->que_suis_je() << finl;

  dom.typer("Rectangle");

  IntTab correspo_new_som;
  creer_sommet_et_elem(dom,ssdomaine,correspo_new_som);

  creer_bords(dom,ssdomaine,correspo_new_som);

  return is;
}



void ConstruireDomaine::creer_sommet_et_elem(Domaine& dom, Sous_Domaine& ssz,IntTab& correspo_newsom)
{
  Domaine& domaine0 = ssz.domaine();
  Domaine& dom0 = domaine0;

  const int nb_poly = ssz.nb_elem_tot();
  const int nb_som_poly = domaine0.nb_som_elem();
  int i,j;
  const DoubleTab& les_coord = dom0.coord_sommets();
  DoubleTab newsom(nb_som_poly*nb_poly,Objet_U::dimension);
  int compteur = 0;
  IntTab new_elems(nb_poly, nb_som_poly); // les nouveaux elements
  correspo_newsom.resize(domaine0.nb_som_tot());
  correspo_newsom = -1;
  IntTab& les_elems_= dom.les_elems();
  les_elems_.ref(new_elems);

  for (i=0; i<nb_poly; i++)
    {
      IntTab indice(nb_som_poly);
      //Cout << "poly = " << i << " " << ssz(i) << finl;
      for (j=0; j<nb_som_poly; j++)
        {
          DoubleTab un_som(1, 3);
          //Cout << "som = " << j << " " << domaine0.sommet_elem(i,j) << finl;

          un_som(0,0) = les_coord(domaine0.sommet_elem(ssz(i),j),0);
          un_som(0,1) = les_coord(domaine0.sommet_elem(ssz(i),j),1);
          //Cout << "coord   = " << un_som(0,0) << " " << un_som(0,1) << finl;

          if (dimension==3)  un_som(0,2) = les_coord(domaine0.sommet_elem(ssz(i),j),2);
          indice(j) = ajouterSom(dom,un_som,newsom,compteur);
          new_elems(i,j) = indice(j);
          correspo_newsom(domaine0.sommet_elem(ssz(i),j)) = indice(j);
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


void ConstruireDomaine::creer_bords(Domaine& domaineraf, Sous_Domaine& ssz, IntTab& correspo_ns)
{
  Domaine& domaine0 = ssz.domaine();


  Bords& mes_faces_bord = domaineraf.faces_bord();
  Elem_geom& elem = domaineraf.type_elem();
  for (auto& itr : domaine0.faces_bord())
    {
      Faces& les_faces=itr.faces();
      Bord& bord=mes_faces_bord.add(Bord());
      Cout << "Dans ConstruireDomaine  Nom bord = " << itr.le_nom() << finl;
      //Cout << "Dans ConstruireDomaine  type face  = " << les_faces.type_face() << finl;
      bord.nommer(itr.le_nom());
      bord.typer_faces(les_faces.type_face());
      bord.associer_domaine(domaineraf);

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
    }

  Bord& bord=mes_faces_bord.add(Bord());
  bord.nommer("interface");
  bord.typer_faces(elem->type_face());
  bord.associer_domaine(domaineraf);
  bord.dimensionner(0);

  Faces mes_faces;
  {
    // bloc a factoriser avec Domaine_VF.cpp :
    Type_Face type_face = domaine0.type_elem()->type_face(0);
    mes_faces.typer(type_face);
    mes_faces.associer_domaine(domaine0);

    Static_Int_Lists connectivite_som_elem;
    const int     nb_sommets_tot = domaine0.nb_som_tot();
    const IntTab&    elements       = domaine0.les_elems();

    construire_connectivite_som_elem(nb_sommets_tot,
                                     elements,
                                     connectivite_som_elem,
                                     1 /* include virtual elements */);

    Faces_builder faces_builder;
    IntTab elem_faces; // Tableau dont on aura pas besoin
    faces_builder.creer_faces_reeles(domaine0,
                                     connectivite_som_elem,
                                     mes_faces,
                                     elem_faces);
  }

  Cerr << "Dans ConstruireDomaine  Le domaine construit nb elem ; nb som = " << domaineraf.nb_elem() << " " << domaineraf.nb_som() << finl;
  Cerr << "Dans ConstruireDomaine  Le domaine construit nb bords = " << domaineraf.nb_faces_bord() << finl;
  Cerr << "Dans ConstruireDomaine  Le domaine construit  bords = " << domaineraf.faces_bord() << finl;
  //Cerr << "Dans ConstruireDomaine  Le domaine0  bords = " << domaine0.faces_bord() << finl;

}
