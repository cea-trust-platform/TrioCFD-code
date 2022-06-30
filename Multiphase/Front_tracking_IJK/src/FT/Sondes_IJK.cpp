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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Sondes_IJK.cpp
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#include <Sondes_IJK.h>
#include <stat_counters.h>

Implemente_liste(Sonde_IJK);
Implemente_instanciable(Sondes_IJK,"Sondes_IJK",LIST(Sonde_IJK));

Sortie& Sondes_IJK::printOn( Sortie& os ) const
{
  //  Objet_U::printOn( os );
  return os;
}

Entree& Sondes_IJK::readOn( Entree& s )
{
  //  assert(mon_post.non_nul());

  Motcle motlu;
  Motcle accolade_ouverte("{");
  Motcle accolade_fermee("}");

  s >> motlu;

  if (motlu != accolade_ouverte)
    {
      Cerr << "Error while reading the probes in the postprocessing" << finl;
      Cerr << "We expected { to start to read the probes" << finl;
      exit();
    }
  s >> motlu;

  if (motlu == accolade_fermee)
    {
      Cerr << "Error while reading the probes in the postprocessing" << finl;
      Cerr << "You have not defined any probe" << finl;
      exit();
    }
  while (motlu != accolade_fermee)
    {
      Sonde_IJK une_sonde;
      une_sonde.nommer(motlu);
      // GB: Le constructeur par le nom n'existe pas ?
      //      Sonde_IJK une_sonde(motlu);
      //      une_sonde.associer_post(mon_post.valeur());
      s >> une_sonde;
      add(une_sonde);
      s >> motlu;
    }
  Cerr << "End of reading probes " << finl;
  return s;
}

void Sondes_IJK::postraiter()
{
  statistiques().begin_count(postraitement_counter_);

  LIST_CURSEUR(Sonde_IJK) curseur=*this;
  while(curseur)
    {
      curseur->postraiter();
      ++curseur;
    }
  statistiques().end_count(postraitement_counter_);

}

void Sondes_IJK::mettre_a_jour(double temps, double tinit)
{
  statistiques().begin_count(postraitement_counter_);

  LIST_CURSEUR(Sonde_IJK) curseur=*this;
  while(curseur)
    {
      curseur->mettre_a_jour(temps,tinit);
      ++curseur;
    }
  statistiques().end_count(postraitement_counter_);

}

//void Sondes_IJK::associer_post(const Postraitement& post)
//{
//  mon_post = post;
//}
