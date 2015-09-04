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
// File:        Modif_Bord_Assemblage.h
// Directory:   $TRUST_ROOT/src/UtilitairesAssemblages
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modif_Bord_Assemblage_included
#define Modif_Bord_Assemblage_included

#include <Interprete_geometrique_base.h>
#include <Zone.h>
class Domaine;
class DoubleTab;

//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//    Classe Modif_Bord_Assemblage
//    Cette classe est un interprete qui sert a lire et executer
//    la directive Modif_Bord_Assemblage :
//        Modif_Bord_Assemblage { nom_domaine nom_cl M0 pas_y nom_fichier }
//    Cette directive est a utiliser en discretisation VEF 2D pour modifier
//    des CL associees a des assemblages de coeur.
// .SECTION voir aussi
//    Interprete Modif_Bord_Assemblage
//    Cette classe est utilisable en 3D
//////////////////////////////////////////////////////////////////////////////
class Modif_Bord_Assemblage : public Interprete_geometrique_base
{
  Declare_instanciable(Modif_Bord_Assemblage);

public :

  Entree& interpreter_(Entree&);

protected :

  void differencier_faces_assemblages(Domaine&) ;

  Nom nom_dom,nom_cl,nom_fichier;
  double entreplat, epaisseur_jeu;
  int M0, nb_couronnes;

};
#endif
