/****************************************************************************
* Copyright (c) 2023, CEA
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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Extraire_surface_ALE.cpp
// Directory : $TRIOCFD_ROOT/src/Fluid_Structure_Interaction/ALE
//
/////////////////////////////////////////////////////////////////////////////

#include <Domaine_ALE.h>
#include <Extraire_surface_ALE.h>
#include <Probleme_base.h>
#include <Equation_base.h>
#include <NettoieNoeuds.h>
#include <Parser_U.h>
#include <Domaine.h>
#include <Domaine_VF.h>
#include <Param.h>

Implemente_instanciable( Extraire_surface_ALE, "Extraire_surface_ALE", Extraire_surface ) ;
// XD Extraire_surface_ALE interprete Extraire_surface_ALE -3 Extraire_surface_ALE in order to extract a surface on a mobile boundary (with ALE desciption). NL2 Keyword to specify that the extract surface is done on a mobile domain. The surface mesh is defined by one or two conditions. The first condition is about elements with Condition_elements. For example: Condition_elements x*x+y*y+z*z<1 NL2 Will define a surface mesh with external faces of the mesh elements inside the sphere of radius 1 located at (0,0,0). The second condition Condition_faces is useful to give a restriction.NL2 By default, the faces from the boundaries are not added to the surface mesh excepted if option avec_les_bords is given (all the boundaries are added), or if the option avec_certains_bords is used to add only some boundaries.
// XD attr domaine ref_domaine domaine 0 Domain in which faces are saved
// XD attr probleme ref_Pb_base probleme 0 Problem from which faces should be extracted
// XD attr condition_elements chaine condition_elements 1 not_set
// XD attr condition_faces chaine condition_faces 1 not_set
// XD attr avec_les_bords rien avec_les_bords 1 not_set
// XD attr avec_certains_bords listchaine avec_certains_bords 1 not_set
Sortie& Extraire_surface_ALE::printOn( Sortie& os ) const
{
  return Extraire_surface::printOn( os );
}

Entree& Extraire_surface_ALE::readOn( Entree& is )
{
  return Extraire_surface::readOn( is );
}
/*Entree& Extraire_surface_ALE::interpreter_(Entree& is)
{

  Extraire_surface::interpreter_(is);

  return is;
}*/


Entree& Extraire_surface_ALE::interpreter_(Entree& is)
{
  Nom nom_pb;
  Nom nom_domaine_surfacique;
  Nom expr_elements("1"),expr_faces("1");
  int avec_les_bords;
  Noms noms_des_bords;
  Param param(que_suis_je());
  param.ajouter("Domaine",&nom_domaine_surfacique,Param::REQUIRED);
  param.ajouter("probleme",&nom_pb,Param::REQUIRED);
  param.ajouter("condition_elements",&expr_elements);
  param.ajouter("condition_faces",&expr_faces);
  param.ajouter_flag("avec_les_bords",&avec_les_bords);
  param.ajouter("avec_certains_bords",&noms_des_bords);
  param.lire_avec_accolades_depuis(is);

  associer_domaine(nom_domaine_surfacique);
  Domaine& domaine_surfacique=domaine();

  if (domaine_surfacique.nb_som_tot()!=0)
    {
      Cerr << "Error!" << finl;
      Cerr <<"The domain " << domaine_surfacique.le_nom() << " can't be used by the Extraire_surface keyword." <<finl;
      Cerr <<"The domain for Extraire_surface keyword should be empty and created just before." << finl;
      exit();
    }

  // on recupere le pb
  if(! sub_type(Probleme_base, objet(nom_pb)))
    {
      Cerr << nom_pb << " is of type " << objet(nom_pb).que_suis_je() << finl;
      Cerr << "and not of type Probleme_base" << finl;
      exit();
    }
  Probleme_base& pb=ref_cast(Probleme_base, objet(nom_pb));
  const Domaine_VF& domaine_vf=ref_cast(Domaine_VF,pb.domaine_dis().valeur());
  const Domaine& domaine_volumique = domaine_vf.domaine();

  // Check that actually a Domaine_ALE:
  if (!(sub_type(Domaine_ALE, domaine_surfacique)))
    {
      Cerr << "Extraire_surface_ALE needs a 'Domaine_ALE' object to work with!!"  << finl;
      Process::exit(-1);
    }

  Domaine_ALE& domaine_surfacique_ale = ref_cast(Domaine_ALE, domaine_surfacique);

  Extraire_surface::extraire_surface_without_cleaning(domaine_surfacique_ale, domaine_volumique,nom_domaine_surfacique,domaine_vf,expr_elements,expr_faces,avec_les_bords,noms_des_bords);

  IntTab les_elems_ref= domaine_surfacique.les_elems();
  //We save the elements belonging to the surface; it contains the global numerations of the faces.
  //Because after cleaning the global numerations of the faces will be lost.
  domaine_surfacique_ale.set_les_elems_extrait_surf_reference(les_elems_ref);
  NettoieNoeuds::nettoie(domaine_surfacique_ale);


  return is;
}

