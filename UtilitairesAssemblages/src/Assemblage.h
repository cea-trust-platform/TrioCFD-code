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
// File:        Assemblage.h
// Directory:   $TRUST_ROOT/src/UtilitairesAssemblages
// Version:     /main/7
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Assemblage_included
#define Assemblage_included

#include <Objet_U.h>
#include <IntLists.h>
class DoubleTab;

//
// .DESCRIPTION class Assemblage
//

//
//

class Assemblage  : public Objet_U
{
  Declare_instanciable_sans_constructeur(Assemblage);
public:
  Assemblage();
  inline void setMN(int m, int n);
  inline void addFace(int type, int i);
  inline IntList& getFaces(int type);
  void reordonner_faces(const DoubleTab& xv);

private :

  IntLists liste_face; //liste_face(type)(nface) liste de face par type (1->6 selon sa position dans l'assemblage) et range dans l'ordre des Z croissants.
  int M;
  int N;
  double x0, y0;


};

void Assemblage::setMN(int m, int n)
{
  M = m;
  N = n;
}


void Assemblage::addFace(int type, int f)
{
  liste_face[type-1].add(f);
}

IntList& Assemblage::getFaces(int type)
{
  return liste_face[type-1];
}

#endif
