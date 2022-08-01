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
// File:        SourceFiltree_FT_disc_base.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////

#include <SourceFiltree_FT_disc_base.h>
#include <Parser.h>
#include <Transport_Interfaces_FT_Disc.h>

SourceFiltree_FT_disc_base::SourceFiltree_FT_disc_base()
{
  dimension_ = Objet_U::dimension;
  fI_xyz_t = NULL;
  temps_ = 0.;
}
SourceFiltree_FT_disc_base::~SourceFiltree_FT_disc_base()
{
  int i;
  for (i=0 ; i<2*dimension_ ; i++)
    {
      delete fI_xyz_t[i];
    }
  delete[] fI_xyz_t;
}

Sortie& SourceFiltree_FT_disc_base::ecrire_donnees(Sortie& os) const
{
  int i;
  os<<"phase 0 : Fonctions : { "<<finl;
  for (i=0 ; i<dimension_ ; i++)
    {
      os << (fI_xyz_t[i]->getString()).c_str() ;
      if (i<dimension_-1)
        {
          os << " , ";
        }
    }
  os << " } "<<finl;
  os<<"phase 1 : Fonctions : { "<<finl;
  for (i=0 ; i<dimension_ ; i++)
    {
      os << (fI_xyz_t[dimension_+i]->getString()).c_str() ;
      if (i<dimension_-1)
        {
          os << " , ";
        }
    }
  os << " } "<<finl;

  return os;
}
// Description:
//    Lit les parametres du terme source a partir
//    d'un flot d'entree.
// Precondition:
// Parametre: Entree& is
//    Signification: un flot d'entree
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour:
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
Entree& SourceFiltree_FT_disc_base::lire_donnees(Entree& is)
{
  Motcle accolade_ouverte("{");
  Motcle accolade_fermee("}");
  Motcle motlu;
  Motcles les_motcles(3);
  {
    les_motcles[0] = "phase0";
    les_motcles[1] = "phase1";
    les_motcles[2] = "equation_transport";
  }

  is >> motlu;
  if (motlu != accolade_ouverte)
    {
      Cerr << "Erreur a la lecture des parametres de lire_donnees " << finl;
      Process::exit();
    }
  is >> motlu;
  Nom tmp;
  fI_xyz_t = new Parser*[2*dimension_];
  int i, ph0_lu,ph1_lu;
  ph0_lu = ph1_lu = 0;
  while (motlu != accolade_fermee)
    {
      int rang = les_motcles.search(motlu);
      switch(rang)
        {
        case 0 :
          {
            is >> motlu;
            ph0_lu = 1;
            assert(motlu==accolade_ouverte);
            for (i=0 ; i<dimension_ ; i++)
              {
                is >> tmp;
                const char *s =  tmp.getChar();
                std::string ss(s);
                for (auto & c: ss) c = toupper(c);
                fI_xyz_t[i] = new Parser(ss,4);
                fI_xyz_t[i]->addVar("x");
                fI_xyz_t[i]->addVar("y");
                fI_xyz_t[i]->addVar("z");
                fI_xyz_t[i]->addVar("t");
                fI_xyz_t[i]->parseString();
                Cerr<<"lecture phase0["<<i<<"] : "<<(fI_xyz_t[i]->getString()).c_str()<<finl;
              }
            is >> motlu;
            assert(motlu==accolade_fermee);
            break;
          }
        case 1 :
          {
            is >> motlu;
            ph1_lu = 1;
            assert(motlu==accolade_ouverte);
            for (i=0 ; i<dimension_ ; i++)
              {
                is >> tmp;
                const char *s =  tmp.getChar();
                std::string ss(s);
                for (auto & c: ss) c = toupper(c);
                fI_xyz_t[dimension_+i] = new Parser(ss,4);
                fI_xyz_t[dimension_+i]->addVar("x");
                fI_xyz_t[dimension_+i]->addVar("y");
                fI_xyz_t[dimension_+i]->addVar("z");
                fI_xyz_t[dimension_+i]->addVar("t");
                fI_xyz_t[dimension_+i]->parseString();
                Cerr<<"lecture phase1["<<dimension_+i<<"] : "<<(fI_xyz_t[dimension_+i]->getString()).c_str()<<finl;
              }
            is >> motlu;
            assert(motlu==accolade_fermee);
            break;
          }
        case 2:
          {
            is >> nom_eq_transport_;
            if (Process::je_suis_maitre())
              Cerr << "Equation de transport pour SourceFiltree_FT_Disc : " << nom_eq_transport_ << finl;
            break;
          }
        default :
          {
            Cerr << "Erreur a la lecture des parametres de SourceFiltree_FT_disc " << finl;
            Cerr << "On attendait les " << les_motcles  << "ou } a la place de " << motlu << finl;
            Process::exit();
          }
        }
      is >> motlu;
    }
  assert(motlu==accolade_fermee);
  if (ph0_lu==0)
    {
      for (i=0 ; i<dimension_ ; i++)
        {
          is >> tmp;
          std::string ss("");
          fI_xyz_t[i] = new Parser(ss,4);
          fI_xyz_t[i]->addVar("x");
          fI_xyz_t[i]->addVar("y");
          fI_xyz_t[i]->addVar("z");
          fI_xyz_t[i]->addVar("t");
          fI_xyz_t[i]->parseString();
        }
    }
  if (ph1_lu==0)
    {
      for (i=0 ; i<dimension_ ; i++)
        {
          is >> tmp;
          std::string ss("0");
          fI_xyz_t[dimension_+i] = new Parser(ss,4);
          fI_xyz_t[dimension_+i]->addVar("x");
          fI_xyz_t[dimension_+i]->addVar("y");
          fI_xyz_t[dimension_+i]->addVar("z");
          fI_xyz_t[dimension_+i]->addVar("t");
          fI_xyz_t[dimension_+i]->parseString();
        }
    }

  ecrire_donnees(Cerr);

  return is;
}
