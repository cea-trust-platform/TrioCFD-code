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
// File:        Pb_MG.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Kernel
// Version:     /main/23
//
//////////////////////////////////////////////////////////////////////////////

#include <Pb_MG.h>
#include <Champ_front_zoom.h>
#include <Schema_Implicite_base.h>
#include <Equation_base.h>
#include <Algo_MG_base.h>
#include <Postraitement.h>

Implemente_instanciable(Pb_MG,"Pb_MG",Couplage_U);

/*! @brief Calcul des connectivites entre faces fines et faces grossieres, et entre elements fins et elements grossiers
 *
 *     pour chaque probleme a 2 grilles
 *     Le calcul est delegue a l'algo.
 *
 */
void Pb_MG::initialize()
{
  if (a_un_pb_grossier())
    pbG_MG().initialize();
  Couplage_U::initialize();
  mon_algo_->initialiser();
  initialiser_champ_front_zoom();
}

void Pb_MG::terminate()
{
  if (a_un_pb_grossier())
    pbG_MG().terminate();
  Couplage_U::terminate();
}

double Pb_MG::computeTimeStep(bool& stop) const
{
  return mon_algo_->computeTimeStep(stop);
}


bool Pb_MG::initTimeStep(double dt)
{
  bool ok=true;
  if (a_un_pb_grossier())
    ok &= pbG_MG().initTimeStep(dt);
  ok&=Couplage_U::initTimeStep(dt);
  return ok;
}

bool Pb_MG::solveTimeStep()
{
  return mon_algo_->solveTimeStep();
}

void Pb_MG::validateTimeStep()
{
  if (a_un_pb_grossier())
    pbG_MG().validateTimeStep();
  Couplage_U::validateTimeStep();
}

bool Pb_MG::isStationary() const
{
  bool stat=Couplage_U::isStationary();
  if (a_un_pb_grossier())
    stat = stat && pbG_MG().isStationary();
  return stat;
}


int Pb_MG::postraiter(int force)
{
  int ok=1;
  if (a_un_pb_grossier())
    ok &= pbG_MG().postraiter(force);
  ok&=Couplage_U::postraiter(force);
  return ok;
}

bool Pb_MG::updateGivenFields()
{
  bool ok=1;
  if (a_un_pb_grossier())
    ok &= pbG_MG().updateGivenFields();
  ok&=Couplage_U::updateGivenFields();
  return ok;
}


///////////////////////////////////////////////////////////
//                                                       //
// Fin de l'implementation de l'interface de Probleme_U  //
//                                                       //
///////////////////////////////////////////////////////////



// PrintOn et ReadOn

Sortie& Pb_MG::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

Entree& Pb_MG::readOn(Entree& s )
{
  return s ;
}


/*! @brief Associe un algorithme de resolution au probleme.
 *
 * Verifie que l'algorithme est derive de Algo_MG_base.
 *
 * @param (un_algo) l'algorithme a associer
 * @return (0 en cas d'erreur, 1 sinon.)
 */
int Pb_MG::associer_algo(const Algo_MG_base& un_algo)
{
  mon_algo_=un_algo;
  return mon_algo_->associer_pb(*this);
}


//INITIALISATION A ZERO DES CHAMPS FRONT ZOOM
//SANS INCIDENCE SUR LE CALCUL ???????????????
void Pb_MG::initialiser_champ_front_zoom()
{
  int i;
  int j = nb_problemes();
  for(i = 0; i<j; i++)
    {
      Probleme_U& pbF = ref_cast(Probleme_U,probleme(i));

      // cas Pb_MG
      if (sub_type(Pb_MG,pbF))
        ref_cast(Pb_MG,pbF).initialiser_champ_front_zoom();

      // cas Probleme_base
      else if (sub_type(Probleme_base,pbF))
        {
          Probleme_base& pbFb=ref_cast(Probleme_base,pbF);
          int nb_eq;
          int nb_eq_tot = pbFb.nombre_d_equations();

          for(nb_eq=0; nb_eq<nb_eq_tot; nb_eq++)
            {
              //  Cerr<<"num de l'equa : "<<nb_eq<<finl;
              Equation_base& eqF = pbFb.equation(nb_eq);
              Zone_Cl_dis_base& zoneF = eqF.zone_Cl_dis();
              int nbCond = zoneF.nb_cond_lim();
              //Cerr<<"nbCond = "<<nbCond<<finl;
              int nb_cond;
              for(nb_cond=0; nb_cond<nbCond; nb_cond++)
                {
                  Champ_front_base& champFront = zoneF.les_conditions_limites(nb_cond).valeur().champ_front();

                  if(sub_type(Champ_front_zoom, champFront))
                    {
                      DoubleTab& tab = champFront.valeurs();
                      tab = 0.;

                    }
                }
            }
        }

      // Autres types de problemes interdits
      else
        {
          Cerr << "A Pb_MG can not contain a " << pbF.que_suis_je() << finl;
          exit();
        }
    }
}


/*! @brief Association du probleme multigrille avec son probleme grossier
 *
 */
void Pb_MG::associer_pbMG_pbGglobal_(Probleme_base& pb)
{
  // Est-ce que cette ligne servait a quelque chose ???
  associer_pbMG_pbG_(pb);
}

void Pb_MG::associer_pbMG_pbG_(Probleme_base& pb)
{
  probleme_grossier=pb;
}


/*! @brief Ajoute un probleme a la liste des problemes fins Met en place la reference du probleme fin vers this
 *
 */
void Pb_MG::associer_pbMG_pbFin_(Probleme_base& pbF)
{

  addProblem(pbF);

  Pb_2G pb2G;
  associer_pb_2G(pb2G); // Attention : add dans la liste fait une recopie !
  problemes_2G(nb_problemes()-1).associer_probleme_fin(*this,nb_problemes()-1);
}


/*! @brief Association du probleme multigrille global avec un de ses problemes a 2 grilles
 *
 */
void Pb_MG::associer_pb_2G(const Pb_2G& pb2G)
{
  problemes_2G.add(pb2G);
}
