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
// File:        Implicit_steady.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Schema_Euler_Implicite_Stationnaire/src/New
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Implicit_steady_included
#define Implicit_steady_included

#include <Simpler.h>
#include <Navier_Stokes_std.h>
#include <Nom.h>

//Description
//Ref. A. P
//De facon identique a l algorithme Implicite (Piso avec avancement_crank_=0) l'algorithme IMPLICITE_STEADY
//se decompose en deux etapes  "prediction" et "correction"
//Le pas de temps dt (un scalaire) est remplace par un vecteur de pas de temps locaux (associe aux faces du maillage)
//Rq: le systeme (B(M/dt)-1Bt) est assemble a chaque iteration

class Implicit_steady : public Simpler
{

  Declare_instanciable_sans_constructeur(Implicit_steady);

public :

  Implicit_steady();
  //Modification of the projection and correction steps for the taking into account of a local time step
  void iterer_NS(Equation_base&, DoubleTab& current, DoubleTab& pression, double, Matrice_Morse&, double, DoubleTrav&,int nb_iter,int& converge, int& ok) override;
  //compute the mass matrix divide by dt_locaux
  void calcul_mat_masse_diviser_par_dt_vef(Navier_Stokes_std& eqnNS, DoubleVect& m_dt, DoubleVect& dt_locaux);
  void calcul_mat_masse_diviser_par_dt_vdf(Navier_Stokes_std& eqnNS, DoubleVect& m_dt, DoubleVect& dt_locaux);
  //void test_periodic_solution(Navier_Stokes_std& eqnNS, DoubleTab& current) const;


  inline const Nom& le_nom() const override
  {
    static Nom nom="Implicit_steady";
    return nom;
  };



protected:


};






#endif

