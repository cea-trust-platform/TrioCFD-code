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
// File:        Pb_MG.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Kernel
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Pb_MG_included
#define Pb_MG_included


#include <Couplage_U.h>
//#include <Pb_grossier.h>
#include <List_Pb_2G.h>
#include <Ref_Pb_2G.h>
#include <Ref_Algo_MG_base.h>
#include <Ref_Probleme_base.h>
#include <Probleme_base.h>

class Pb_2G;
class Equation_base;
class Champ_front_zoom;


/*! @brief class Pb_MG
 *
 */

//////////////////////////////////////////////////////////////////////////////////
//
// CLASS: Pb_MG
//
/////////////////////////////////////////////////////////////////////////////////

class Pb_MG : public Couplage_U
{
  Declare_instanciable(Pb_MG);

public:

  //////////////////////////////////////////////////
  //                                              //
  // Implementation de l'interface de Probleme_U  //
  //                                              //
  //////////////////////////////////////////////////

  void initialize() override;
  void terminate() override;
  double computeTimeStep(bool& stop) const override;
  bool initTimeStep(double dt) override;
  bool solveTimeStep() override;
  void validateTimeStep() override;
  bool isStationary() const override;

  int postraiter(int force=1) override;
  bool updateGivenFields() override;

  ///////////////////////////////////////////////////////////
  //                                                       //
  // Fin de l'implementation de l'interface de Probleme_U  //
  //                                                       //
  ///////////////////////////////////////////////////////////

  // Methodes d'acces aux membres :
  virtual int associer_algo(const Algo_MG_base&);

  // Methodes d'acces inline
  inline int nb_pb_2G_() const ;
  inline Pb_2G& pb_2G(int);
  inline const Pb_2G& pb_2G(int) const ;
  inline LIST(Pb_2G)& problemes_2G_();

  // Methodes d'association
  void associer_pbMG_pbGglobal_(Probleme_base&);
  void associer_pbMG_pbG_(Probleme_base&);
  void associer_pbMG_pbFin_(Probleme_base&);
  void associer_pb_2G(const Pb_2G&);

  // Methodes d'acces au probleme grossier
  int a_un_pb_grossier()
  {
    return probleme_grossier.non_nul();
  }
  int a_un_pb_grossier() const
  {
    return probleme_grossier.non_nul();
  }
  inline Probleme_base& pbG_MG();
  inline const Probleme_base& pbG_MG() const ;

  void initialiser_champ_front_zoom();

private:

  REF(Probleme_base) probleme_grossier;
  LIST(Pb_2G) problemes_2G;
  REF(Algo_MG_base) mon_algo_;

};

// Methodes d'acces inline

inline LIST(Pb_2G)& Pb_MG::problemes_2G_()
{
  return problemes_2G;
}

inline int Pb_MG::nb_pb_2G_() const
{
  return problemes_2G.size();
}

inline Pb_2G& Pb_MG::pb_2G(int i)
{
  return (problemes_2G(i));
}

inline const Pb_2G& Pb_MG::pb_2G(int i) const
{
  return (problemes_2G(i));
}


inline Probleme_base& Pb_MG::pbG_MG()
{
  return probleme_grossier;
}

inline const Probleme_base& Pb_MG::pbG_MG() const
{
  return probleme_grossier;
}

#endif
