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
// File:        Algo_Couple_1.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Algos
// Version:     /main/15
//
//////////////////////////////////////////////////////////////////////////////


#include <Algo_Couple_1.h>
#include <Pb_MG.h>
#include <Schema_Temps_base.h>

Implemente_instanciable(Algo_Couple_1,"Algo_Couple_1",Algo_MG_base);

Sortie& Algo_Couple_1::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

Entree& Algo_Couple_1::readOn(Entree& is )
{
  dt_unif = 0;
  Motcles les_mots(1);
  les_mots[0] = "dt_uniforme";

  Motcle motlu, accolade_fermee="}", accolade_ouverte="{";
  is >> motlu;
  if(motlu!=accolade_ouverte)
    {
      Cerr << "On attendait une { a la lecture d'une " << que_suis_je() << finl;
      Cerr << "et non : " << motlu << finl;
      exit();
    }
  is >> motlu;
  while (motlu != accolade_fermee)
    {
      int rang=les_mots.search(motlu);
      switch(rang)
        {
        case 0 :
          dt_unif = 1;
          break;
        default :
          Cerr << "Mot cle inconnu dans Algo_Couple_1  " << finl;
          exit();
        }
      is >> motlu;
    }
  return is;
}


double Algo_Couple_1::computeTimeStep(bool& stop) const
{

  //Pour ne pas manipuler directement dt_
  double dt;
  ////dt_=DMAXFLOAT;
  dt=DMAXFLOAT;
  const Pb_MG& pbMG=ref_cast(Pb_MG,pb());
  if (dt_unif)
    {
      // dt est le min
      // on s'arrete si quelqu'unveut s'arreter
      if (pbMG.a_un_pb_grossier())
        {
          ////dt_=pbMG.pbG_MG().computeTimeStep(stop);
          dt=pbMG.pbG_MG().computeTimeStep(stop);
          if (stop)
            return 0;
        }
      for (int i=0; i<pbMG.nb_problemes(); i++)
        {
          double dt1=pbMG.probleme(i).computeTimeStep(stop);
          if (stop)
            {
              return 0;
            }
          ///if (dt1<dt_)
          ////dt_=dt1;
          if (dt1<dt)
            dt=dt1;
        }
    }

  else
    {
      ////dt_=pbMG.pbG_MG().computeTimeStep(stop);
      dt=pbMG.pbG_MG().computeTimeStep(stop);
      if (stop)
        return 0;
    }

  ////return dt_;
  return dt;
}


bool Algo_Couple_1::solveTimeStep()
{
  Pb_MG& pbMG=ref_cast(Pb_MG,pb());
  Probleme_base& pbG =pbMG.pbG_MG();

  if (dt_unif == 1)
    {
      pbG.solveTimeStep();
      for(int i=0; i<pbMG.nb_problemes(); i++)
        pbMG.probleme(i).solveTimeStep();
    }
  else
    {
      double t_arrive=pb().presentTime()+dt_;
      pbG.runUntil(t_arrive);
      for(int i=0; i<pbMG.nb_problemes(); i++)
        pbMG.probleme(i).runUntil(t_arrive);
    }

  return true;
}


//Initialisation des prol et restr
//NE FONCTIONNE QUE POUR DES PB DONT
//LES EQUATIONS ONT LE MEME TYPE D'INCONNUES!!!!!!
int Algo_Couple_1::initialiser()
{
  Pb_MG& pbMG=ref_cast(Pb_MG,pb());

  //Pour chaque probleme a 2 grilles
  for(int i=0; i<pbMG.nb_problemes(); i++)
    {
      Pb_2G& pb2G = pbMG.pb_2G(i);
      // Typage des connectivites
      pb2G.typer_Connectivites("Connectivites_faces_couple");
      pb2G.set_nb_prol(0);
      pb2G.set_nb_rest(0);
      // Calcul des connectivites
      pb2G.calculer_connectivites_2G();
    }

  return 1;
}
