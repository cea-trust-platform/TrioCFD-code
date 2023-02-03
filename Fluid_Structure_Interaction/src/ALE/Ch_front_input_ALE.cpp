/****************************************************************************
 * Copyright (c) 2015 - 2018, CEA
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
// File:        Ch_front_input_ALE.cpp
// Directory:   New class
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////

#include <Ch_front_input_ALE.h>
#include <Probleme_base.h>
#include <Schema_Temps_base.h>
#include <Interprete.h>
#include <Domaine.h>
#include <Frontiere_dis_base.h>
#include <communications.h>
#include <Convert_ICoCoTrioField.h>
#include <MD_Vector_tools.h>
#include <Roue.h>

Implemente_instanciable_sans_constructeur(Ch_front_input_ALE,"Ch_front_input_ALE",Ch_front_input_P1);
// XD Ch_front_input_ALE  front_field_base  Ch_front_input_ALE  -1 Class to define a boundary condition on a moving boundary of a mesh (only for the Arbitrary Lagrangian-Eulerian  framework ) .  NL2  Example:  Ch_front_input_ALE { nb_comp 3 nom VITESSE_IN_ALE probleme pb initial_value 3 1. 0. 0. }

Ch_front_input_ALE::Ch_front_input_ALE():
  alreadyInit_(false)
{
}

Entree& Ch_front_input_ALE::readOn(Entree& is)
{
  return Ch_front_input::readOn(is);
}

Sortie& Ch_front_input_ALE::printOn(Sortie& os) const
{
  return Ch_front_input::printOn(os);
}

void Ch_front_input_ALE::getTemplate(TrioField& afield) const
{
  buildSommetsFaces(); //From Ch_front_input in order to update mesh
  Ch_front_input_P1::getTemplate(afield);
}

void Ch_front_input_ALE::setValue(const TrioField& afield)
{
  int NrOfNodesIn_les_valeurs_som=les_valeurs_som[1].valeurs().size()/nb_comp();
  for(int i=0; i<NrOfNodesIn_les_valeurs_som; i++)
    {
      for(int j=0; j<afield._nb_field_components; j++)
        {
          les_valeurs_som[1].valeurs()(i,j)=afield._field[afield._nb_field_components*i+j]; //Values to nodes.
        }
    }
  les_valeurs_som[1].valeurs().echange_espace_virtuel() ;
}

int Ch_front_input_ALE::initialiser(double temps, const Champ_Inc_base& inco)
{
  Cout << "initialiser Ch_front_input_ALE" << finl ;
  if (alreadyInit_)
    return 1;
  alreadyInit_ = true;
  if (!Ch_front_input_P1::initialiser(temps,inco))
    return 0;

  //ALE part:
  mettre_a_jour(temps);

  return 1;
}

void Ch_front_input_ALE::mettre_a_jour(double temps)
{
  Cerr << "Champ_front_input_ALE ::mettre_a_jour" << finl;

  const Frontiere& front=la_frontiere_dis->frontiere();
  int nb_faces=front.nb_faces();
  const Faces& faces=front.faces();
  int nbsf=faces.nb_som_faces();
  int i,j,k;

  DoubleTab& tab=valeurs_au_temps(temps);
  remplir_vit_som_bord_ALE(temps);
  tab=0.;

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

void Ch_front_input_ALE::remplir_vit_som_bord_ALE(double tps)
{
  Cout<<"Champ_front_input_ALE::remplir_vit_som_bord_ALE"<<finl;

  const Frontiere& front=la_frontiere_dis->frontiere();
  //int nb_faces=front.nb_faces();
  const Domaine& domaine=front.domaine();
  const Faces& faces=front.faces();
  //int nbsf=faces.nb_som_faces();
  int i,j;
  int nb_som_tot=domaine.nb_som_tot();
  vit_som_bord_ALE.resize(nb_som_tot,nb_comp());
  //int nb_som=domaine.nb_som();
  //if (vit_som_bord_ALE.dimension(0) != domaine.nb_som())
  //  {
  //    vit_som_bord_ALE.resize(domaine.nb_som(),nb_comp());
  //    const MD_Vector& md = domaine.domaine().md_vector_sommets();
  //    MD_Vector_tools::creer_tableau_distribue(md, vit_som_bord_ALE);
  //  }


  vit_som_bord_ALE=0.;
  double InitialTimeValue=mon_pb->schema_temps().temps_init();

  // Construction de la liste des sommets
  int nn=0;
  ArrOfInt liste_sommets(nb_som_tot);
  ArrOfInt marqueur(domaine.les_sommets().dimension(0));
  const IntTab& ffaces=faces.les_sommets() ;
  marqueur=-1;
  for (int f=0; f<ffaces.dimension(0); f++)
    {
      for (int s=0; s<ffaces.dimension(1); s++)
        {
          int som=ffaces(f,s);
          if (som >= 0 && marqueur[som]==-1)
            {

              marqueur[som]=nn;
              liste_sommets[nn]=som ;
              nn++;
            }
        }
    }

  /*Cout << "nb_sommets " << sommets_.size() << finl ;
  Cout << "Sommets: " ;
  for (i=0 ; i < nn ; i++)
  { Cout << liste_sommets(i) << ", " ;}
  Cout << finl ;*/

  if(tps==InitialTimeValue)  //Initial values.
    {
      //    for( i=0; i<nb_som; i++)
      for( i=0; i<nn; i++)
        {
          for( j=0; j<nb_comp(); j++)
            {
              vit_som_bord_ALE(liste_sommets[i],j)=initial_value_[j];
            }
        }
    }
  else
    {
      //    for( i=0; i<nb_som; i++)
      for( i=0; i<nn; i++)
        {
          for( j=0; j<nb_comp(); j++)
            {
              vit_som_bord_ALE(liste_sommets[i],j)=les_valeurs_som[1].valeurs()(i,j);
              //std::cout << "i : " << i << ", j: " << j << ",  les_valeurs_som: " << les_valeurs_som[1].valeurs()(i,j)  << std::endl ;
            }
        }
    }
  //vit_som_bord_ALE.echange_espace_virtuel();
}

