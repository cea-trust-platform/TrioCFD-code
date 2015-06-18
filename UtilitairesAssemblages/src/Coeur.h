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
// File:        Coeur.h
// Directory:   $TRUST_ROOT/src/UtilitairesAssemblages
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Coeur_included
#define Coeur_included


#include <Zone.h>
#include <Vect_Assemblage.h>
class Domaine;
class DoubleTab;
//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
// .SECTION voir aussi
//    Interprete CalculerConnectFacesHexaRNR
//    Cette classe est utilisable en VEF/3D
//////////////////////////////////////////////////////////////////////////////
class Coeur : public Objet_U
{
  Declare_instanciable(Coeur);

public :

  void calculer() ;
  void calculer_geom() ;

  inline Assemblage& getAssemblage(int m, int n);
  inline int getAssemblage_size(int m, int n, int i);
  inline int getAssemblage_face(int m, int n, int i, int j);
  inline void setEntreplat(double e);
  inline void setEpaisseurJeu(double e);
  inline void setNbCouronnes(int nb);
  inline void setOrigineNum(int m0);
  inline void setTypeProbleme(int t);
  inline void setBord(const Nom& n);
  inline void setProbleme(const Nom& n);
  inline void setDomaine(const Nom& n);

protected:

private:

  int test;

  double ep; // entreplat d'un assemblage
  double pas_x, pas_y; //pas entre les centres de deux assemblages
  double ep_jeux; //epaisseur des jeux inter-assemblages
  int nb_couronnes; //nb de couronnes sans compter l'assemblage central
  Nom nom_bord; // nom du bord des assemblages
  Nom nom_pb; // nom du probleme des assemblages
  Nom nom_dom;
  int M0; // origine de la numeration des assemblages. Exemple : mettre 20 pour une origine a (20,20) ou 30 pour une origine a (30x30)

  // PQ : 10/10/08 : Obsolete ?
  int jeux_ou_assemblages; // 0 si le domaine fourni correspond aux jeux, 1 s'il correspond aux assemblages.

  VECT(Assemblage) vect_assemblages;

};

inline Assemblage& Coeur::getAssemblage(int m, int n)
{
  int iass = (m-M0+nb_couronnes)*(2*nb_couronnes+1)+(n-M0+nb_couronnes);
  return vect_assemblages[iass];
}

inline int Coeur::getAssemblage_size(int m, int n, int i)
{
  int iass = (m-M0+nb_couronnes)*(2*nb_couronnes+1)+(n-M0+nb_couronnes);
  return vect_assemblages[iass].getFaces(i).size();
}

inline int Coeur::getAssemblage_face(int m, int n, int i, int j)
{
  int iass = (m-M0+nb_couronnes)*(2*nb_couronnes+1)+(n-M0+nb_couronnes);
  return vect_assemblages[iass].getFaces(i)[j];
}

inline void Coeur::setEntreplat(double e)
{
  ep = e;
}

inline void Coeur::setEpaisseurJeu(double e)
{
  ep_jeux = e;
}

inline void Coeur::setNbCouronnes(int nb)
{
  nb_couronnes = nb;
}

inline void Coeur::setOrigineNum(int m0)
{
  M0 = m0;
}

inline void Coeur::setTypeProbleme(int t)
{
  jeux_ou_assemblages = t;
}

inline void Coeur::setBord(const Nom& n)
{
  nom_bord = n;
}

inline void Coeur::setProbleme(const Nom& n)
{
  nom_pb = n;
}

inline void Coeur::setDomaine(const Nom& n)
{
  nom_dom = n;
}


#endif
