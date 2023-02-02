/****************************************************************************
* Copyright (c) 2019, CEA
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
// File      : Random_process.h
// Directory : $IJK_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Random_process_included
#define Random_process_included

#include <FixedVector.h>
#include <IJK_Splitting.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <Force_ph.h>
#include <random>
#include <Objet_U.h>
#include <Redistribute_Field.h>
#include <Parser.h>
#include <IJK_Splitting.h>
#include <Multigrille_Adrien.h>
#include <Interprete.h>
#include <TRUST_Ref.h>

class Probleme_base;

// #include <Force_sp.h>
// #include <OpDiffTurbIJK.h>
// #include <OpConvIJKQuickSharp.h>
// #include <OpCentre4IJK.h>
// #include <OpConvIJKAmont.h>
// #include <Linear_algebra_tools.h>
// #include <Boundary_Conditions.h>
// #include <IJK_Interfaces.h>
// #include <IJK_FT_Post.h>


class Random_process : public Objet_U
{
  Declare_instanciable_sans_constructeur_ni_destructeur( Random_process ) ;

public:

  Random_process();
  ~Random_process();
  void initialise(double eps_etoile, double tL,
                  int nl, int nm, int nn, std::string nom_fichier_sortie, Nom nom_sauvegarde);//,int random_fixed_);
  void next_step(double dt, int it);
  void next_step2(double dt, int it);
  void next_step3(ArrOfDouble& advection_velocity, double dt, int it);
  void write(std::string nom_fichier_sortie, double temps);
  void write_separate(std::string nom_fichier_sortie, double t);
  std::vector< std::vector< std:: vector < double > > > get_b();
  // std:: vector < double > get_b_flt();
  ArrOfDouble& get_b_flt();
  std::minstd_rand initialise_gen(int i);
  int get_semi_gen();

private:

  int nl,nm,nn,n_lmn;
  double eps_etoile,tL;
  int kmin,kmax;

  std::minstd_rand gen;
  std::string nom_fichier_;
  Nom nom_sauvegarde_;
//  IntTab semi_gen_et_modulo_reprise_;
  ArrOfInt semi_gen_et_modulo_reprise_;
  std::normal_distribution < double > distribution;
  int moke_gen_;

  std::vector< std::vector< std:: vector < double > > > process;
  ArrOfDouble process_flt;

//  int i_offset;
//  int j_offset;
//  int k_offset;

// je voudrai declarer et me servir de ces attributs
//  std::string const nomFichierReprise_("reprise_gen.sauv");
//  std::string nomFichierReprise_;
//  std::ofstream *gen_write_;
//  std::ofstream gen_write_;
//  std::ifstream *gen_read_;
//  std::ifstream gen_read_;
};


#endif
