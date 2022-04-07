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
// File:        Prolongement_elem_elem_FMG.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Operateurs
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////

#include <Prolongement_elem_elem_FMG.h>
#include <Domaine.h>
#include <Zone_VF.h>

Implemente_instanciable(Prolongement_elem_elem_FMG,"Prolongement_elem_elem_FMG",Prolongement_base);

//// printOn
//

Sortie& Prolongement_elem_elem_FMG::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

//// readOn
//

Entree& Prolongement_elem_elem_FMG::readOn(Entree& s )
{
  return s ;
}




// Description:
//    Prolongement du centre de gravite des faces grossieres
//    au centre de gravite de TOUTES les faces fines
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour:
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:

//ATTENTION
//NE FONCTIONNE QU'EN 2D VDF - VDF!!!!!!!!!!!!!!
void Prolongement_elem_elem_FMG::prolonger(Zone_VF& zone_VFG,
                                           Zone_VF& zone_VFF,
                                           const Frontiere& frontF,
                                           IntVect& connect,
                                           const DoubleTab& valG, DoubleTab& tab,
                                           int nb_comp)
{
  int nb_elemF = connect.size_array();
  int nbelemsF;
  int num_elemG;
  double epsilonP = 0;
  int nbfaceG;
  Zone& zoneg = zone_VFG.zone();
  int nb_faceG = zoneg.nb_faces_elem();
  Zone& zonef = zone_VFF.zone();
  Elem_geom& elem_geomF = zonef.type_elem();

  int nb_sommets_elemG = zoneg.nb_som_elem();
  IntVect numerosG1;
  DoubleVect numerosG1_dist;
  IntVect numerosG2;
  DoubleVect numerosG2_dist;
  numerosG1.resize(nb_faceG+1);
  numerosG2.resize(nb_faceG+1);
  numerosG1_dist.resize(nb_faceG+1);
  numerosG2_dist.resize(nb_faceG+1);
  int numG, num_faceG;
  IntTab& num_face_elemG = zone_VFG.elem_faces();

  //coord des centres de gravite des elems
  DoubleTab& cgF = zone_VFF.xp();
  DoubleTab& cgG = zone_VFG.xp();
  double dist;
  tab = 0.;

  Domaine& domG = zoneg.domaine();
  const DoubleTab& coordG = domG.coord_sommets();
  DoubleVect coord_somG;
  coord_somG.resize(dimension);
  int dim;
  ArrOfInt elem;
  int compo;
  int numG2 = -1;
  IntTab& sommet_elemG = zoneg.les_elems();

  //elemG voisins de chaque face grossiere
  IntTab& voisins = zone_VFG.face_voisins();

  //Pour chaque elemF
  for (nbelemsF=0; nbelemsF<nb_elemF; nbelemsF++)
    {
      // ELEM G CORRESPONDANT
      // on recupere le numero de l'elemG correspondant
      num_elemG = connect(nbelemsF);
      //on ajoute la valeur valG(num_elemG)
      if (nb_comp == 1)
        {
          tab(nbelemsF) = tab(nbelemsF) + (1-1.5*epsilonP)*valG(num_elemG);
        }
      else
        {
          for(compo = 0; compo<nb_comp; compo++)
            tab(nbelemsF, compo) = tab(nbelemsF, compo) + (1-1.5*epsilonP)*valG(num_elemG, compo);
        }


      //ELEM VOISIN G CORRESPONDANT
      //Pour chaque faceG
      for (nbfaceG=0; nbfaceG<nb_faceG; nbfaceG++)
        {
          num_faceG = num_face_elemG(num_elemG, nbfaceG);

          //on trouve l'elemG voisin de num_elemG
          numG = voisins(num_faceG, 0);

          if (numG == num_elemG )
            {
              numG = voisins(num_faceG, 1);   // element exterieur
            }
          //  if (numG != num_elemG)
          numerosG1(nbfaceG) = numG;
          double som_dist = 0.;
          dist = 0.;

          for(dim=0; dim<dimension; dim++)
            {
              if (numG != -1)
                som_dist = som_dist + (cgG(numG, dim)-cgF(nbelemsF, dim))*(cgG(numG, dim)- cgF(nbelemsF, dim));
            }
          dist = sqrt(som_dist);
          numerosG1_dist(nbfaceG) = dist;
        }

      // On classe les elemG du plus proche au moins proche de l'elemF
      // on copie la version  G1
      for (nbfaceG=0; nbfaceG<nb_faceG; nbfaceG++)
        {
          numerosG2(nbfaceG) = numerosG1(nbfaceG);
          numerosG2_dist(nbfaceG) = numerosG1_dist(nbfaceG);
        }

      for (nbfaceG=1; nbfaceG<numerosG2.size()-1; nbfaceG++)
        {
          double distbis = numerosG2_dist(nbfaceG);
          int elem2 = numerosG2(nbfaceG);
          int i = nbfaceG -1;
          // dist doit etre insere dans le tableau ordonne 0..j-1
          while  ( i>=0 && numerosG2_dist[i]>distbis )
            {
              numerosG2_dist(i+1) = numerosG2_dist(i);
              numerosG2(i+1)      = numerosG2(i);
              i = i-1;
            }
          numerosG2_dist(i+1) = distbis;
          numerosG2(i+1)      =  elem2;
        }

      //somme
      // les deux premiers du tableau ont le coef B (1/32 +1/4epsilon)
      // les deux derniers le coeff C (-3/32 +1/4epsilon)

      for (nbfaceG=0; nbfaceG<2; nbfaceG++)
        {
          if (nb_comp == 1)
            {
              if (numerosG2(nbfaceG) != -1)
                tab(nbelemsF) = tab(nbelemsF) + (1./32. + 0.25*epsilonP)*valG(numerosG2(nbfaceG));
            }
          else
            {
              if (numerosG2(nbfaceG) != -1)
                for(compo = 0; compo<nb_comp; compo++)
                  tab(nbelemsF, compo) = tab(nbelemsF, compo) + (1./32. + 0.25*epsilonP)*valG(numerosG2(nbfaceG), compo);
            }
        }
      for (nbfaceG=2; nbfaceG<nb_faceG; nbfaceG++)
        {
          if (nb_comp == 1)
            {
              if (numerosG2(nbfaceG) != -1)
                tab(nbelemsF) = tab(nbelemsF) + (-3./32. + 0.25*epsilonP)*valG(numerosG2(nbfaceG));
            }
          else
            {
              if (numerosG2(nbfaceG) != -1)
                for(compo = 0; compo<nb_comp; compo++)
                  tab(nbelemsF, compo) = tab(nbelemsF, compo) + (-3./32. + 0.25*epsilonP)*valG(numerosG2(nbfaceG), compo);
            }
        }

      //ELEM G POUR "EQUILIBRE"
      //si l'elemF contient un sommet de num_elemG
      //on trouve l'elemG contenant ce sommet
      //et qui n'a pas encore ete utilise
      int trouve = 0;
      int nb_som = 0;
      while ((trouve == 0) && (nb_som<nb_sommets_elemG))
        {
          //coord du sommetG nb_som
          int num_som = sommet_elemG(num_elemG, nb_som);
          for(dim=0; dim<dimension; dim++)
            {
              coord_somG(dim) = coordG(num_som, dim);
            }
          //on regarde si l'elemF contient le sommet nb_som
          trouve = elem_geomF.contient(coord_somG, nbelemsF);
          if (trouve == 0)
            nb_som = nb_som+1;
        }
      //on trouve les numG contenant ce sommetG
      if (trouve == 1)
        {
          if(dimension == 2)
            zoneg.rang_elems_sommet(elem, coord_somG(0), coord_somG(1), 0);
          else if (dimension == 3)
            zoneg.rang_elems_sommet(elem, coord_somG(0), coord_somG(1), coord_somG(2));
        }

      //on trouve le numG different de ceux deja utilises
      int taille = elem.size_array();
      trouve = 0;
      int n=0;
      int k=0;
      int trouve2;

      while ((n<taille) && (trouve == 0))
        {
          if (elem[n] == num_elemG)
            {
              trouve = 1;
              trouve2 = 1;
            }
          else
            {
              k = 0;
              trouve2 = 0;
            }
          while((trouve2 == 0) && (k<nbfaceG))
            {
              if (elem[n] == numerosG2(k))
                trouve2 = 1;
              else
                k=k+1;
            }
          if(trouve2 == 1)
            {
              trouve = 0;
              n = n+1;
            }
          else
            {
              trouve = 1;
              numG2 = elem[n];
            }
        }
      if (nb_comp == 1)
        {
          if (numG2 != -1)
            tab(nbelemsF) = tab(nbelemsF) + (1./8.)*valG(numG2);
        }
      else
        {
          if (numG2 != -1)
            for(compo = 0; compo<nb_comp; compo++)
              tab(nbelemsF, compo) = tab(nbelemsF, compo) + (1./8.)*valG(numG2, compo);
        }
    }
}




//NE FAIT RIEN
void Prolongement_elem_elem_FMG::calculer(Zone_VF& zonef,
                                          Zone_VF& zoneg,
                                          IntVect& connect_ff)
{
  //ne fait rien mais c'est normal!!!
}
