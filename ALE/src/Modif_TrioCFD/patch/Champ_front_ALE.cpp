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
// File:        Champ_front_ALE.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/ALE/src
// Version:     /main/17
//
//////////////////////////////////////////////////////////////////////////////

#include <Champ_front_ALE.h>
#include <Domaine.h>
#include <Frontiere_dis_base.h>
#include <MD_Vector_tools.h>

Implemente_instanciable(Champ_front_ALE,"Champ_front_ALE",Ch_front_var_instationnaire_dep);
//XD Champ_front_ale front_field_base Champ_front_ale 0 Class to define a boundary condition on a moving boundary of a mesh (only for the Arbitrary Lagrangian-Eulerian framework ).
//XD attr val listchaine val 0 NL2 Example:  2 -y*0.01 x*0.01

/*Champ_front_ALE::Champ_front_ALE()
{
  const Frontiere& front=la_frontiere_dis->frontiere();
  const Zone& zone=front.zone();
  const Domaine& domaine=zone.domaine();
  vit_som_bord_ALE.resize(domaine.nb_som(),nb_comp());
  const MD_Vector& md = zone.domaine().md_vector_sommets();
  MD_Vector_tools::creer_tableau_distribue(md, vit_som_bord_ALE);
}*/

// Description:
//    Impression sur un flot de sortie au format:
//    taille
//    valeur(0) ... valeur(i)  ... valeur(taille-1)
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
Sortie& Champ_front_ALE::printOn(Sortie& os) const
{
  //   const DoubleTab& tab=valeurs();
  //   os << tab.size() << " ";
  //   for(int i=0; i<tab.size(); i++)
  //     os << tab(0,i);
  return os;
}




// Description:
//    Lecture a partir d'un flot d'entree au format:
//    nombre_de_composantes
//    fonction dependant de xyzt pour chaque composante
// Precondition:
// Parametre: Entree& is
//    Signification: un flot d'entree
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: Entree&
//    Signification: le flot d'entree modifie
//    Contraintes:
// Exception: accolade ouvrante attendue
// Exception: mot clef inconnu a cet endroit
// Exception: accolade fermante attendue
// Effets de bord:
// Postcondition:
Entree& Champ_front_ALE::readOn(Entree& is)
{
  int dim;
  is >> dim;
  fixer_nb_comp(dim);

  fxyt.dimensionner(dim);


  //        Cout << "dim = " << dim << finl;
  for (int i = 0; i<dim; i++)
    {
      Nom tmp;
      //Cout << "i = " << i << finl;
      is >> tmp;
      //Cout << "fonc = " << tmp << finl;
      Cerr << "Lecture et interpretation de la fonction " << tmp << finl;
      fxyt[i].setNbVar(3);
      fxyt[i].setString(tmp);
      fxyt[i].addVar("x");
      fxyt[i].addVar("y");
      //         fxyt[i].addVar("z");
      fxyt[i].addVar("t");
      fxyt[i].parseString();
      Cerr << "Interpretation de la fonction " << tmp << " Ok" << finl;
      //        Cout << "end = " << tmp << finl;
    }
  return is;
}


// Description:
//    Pas code !!
// Precondition:
// Parametre: Champ_front_base& ch
//    Signification:
//    Valeurs par defaut:
//    Contraintes: reference constante
//    Acces: NON ACCEDE
// Retour: Champ_front_base&
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
Champ_front_base& Champ_front_ALE::affecter_(const Champ_front_base& ch)
{
  return *this;
}

int Champ_front_ALE::initialiser(double temps, const Champ_Inc_base& inco)
{
  if (!Ch_front_var_instationnaire_dep::initialiser(temps,inco))
    return 0;

  mettre_a_jour(temps);
  return 1;
}

// Description:
//     Mise a jour du champ front et remplie un tableau de dimension
//     egale au nombre total de sommets du domaine. Ce tableau a des
//     valeurs nulles aux sommets qui n'appartiennent pas au bord ALE
//     traite.
// Precondition:
// Parametre: double tps
//    Signification: le temps de mise a jour
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour:
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
void Champ_front_ALE::mettre_a_jour(double temps)
{
  Cerr << "Champ_front_ALE ::mettre_a_jour" << finl;

  const Frontiere& front=la_frontiere_dis->frontiere();
  int nb_faces=front.nb_faces();
  const Faces& faces=front.faces();
  int nbsf=faces.nb_som_faces();
  int i,j,k;
  DoubleTab& tab=valeurs_au_temps(temps);
  tab=0.;

  // PQ : 13/06/13 appel systematique a remplir_vit_som_bord_ALE() sinon pas bon !!
  remplir_vit_som_bord_ALE(temps);

  for( i=0; i<nb_faces; i++)
    {
      for( j=0; j<nb_comp(); j++)
        {
          for( k=0; k<nbsf; k++)
            {
              tab(i,j)+=vit_som_bord_ALE(faces.sommet(i,k),j);
            }
        }
    }
  tab/=nbsf;
  tab.echange_espace_virtuel();
}

void Champ_front_ALE::remplir_vit_som_bord_ALE(double tps)
{
  Cout<<"Champ_front_ALE::remplir_vit_som_bord_ALE"<<finl;
  const Frontiere& front=la_frontiere_dis->frontiere();
  int nb_faces=front.nb_faces();
  const Zone& zone=front.zone();
  const Faces& faces=front.faces();
  const Domaine& domaine=zone.domaine();
  double x,y,z;
  int nbsf=faces.nb_som_faces();
  int i,j,k;
  int nb_som_tot=domaine.nb_som_tot();
  vit_som_bord_ALE.resize(nb_som_tot,nb_comp());
  /*if (vit_som_bord_ALE.dimension(0) != domaine.nb_som())
    {
      vit_som_bord_ALE.resize(domaine.nb_som(),nb_comp());
      const MD_Vector& md = zone.domaine().md_vector_sommets();
      MD_Vector_tools::creer_tableau_distribue(md, vit_som_bord_ALE);
    }*/

  vit_som_bord_ALE=0.;

  for( i=0; i<nb_faces; i++)
    {
      x=y=z=0;
      for( k=0; k<nbsf; k++)
        {
          x=domaine.coord(faces.sommet(i,k),0);
          if(dimension>1)
            y=domaine.coord(faces.sommet(i,k),1);
          if(dimension>2)
            z=domaine.coord(faces.sommet(i,k),2);
          for( j=0; j<nb_comp(); j++)
            {
              fxyt[j].setVar("x",x);
              fxyt[j].setVar("y",y);
              //               fxyt[j].setVar("z",z);
              fxyt[j].setVar("t",tps);
              vit_som_bord_ALE(faces.sommet(i,k),j)=fxyt[j].eval();
              //cout << " x y  " << x << " " << y << " " << z << " " << vit_som_bord_ALE(faces.sommet(i,k),j) << endl;
            }
        }
    }
  //vit_som_bord_ALE.echange_espace_virtuel();
}
