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
// File:        Source_rayo_semi_transp.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Source_rayo_semi_transp_included
#define Source_rayo_semi_transp_included

#include <Source_rayo_semi_transp_base.h>
#include <TRUST_Deriv.h>
class Modele_rayo_semi_transp;


/*! @brief class Source_rayo_semi_transp Classe generique de la hierarchie des termes sources de l'eqution de
 *
 *     temperature pour les problemes de rayonnement semi transparent
 *     Source_rayo_semi_transp peut referencer n'importe quel
 *     objet derivant de Source_rayo_semi_transp_base.
 *
 * @sa Source_rayo_semi_transp_base
 */
class Source_rayo_semi_transp : public DERIV(Source_rayo_semi_transp_base)
{
  Declare_instanciable(Source_rayo_semi_transp);

public :

  inline Modele_rayo_semi_transp& Modele();
  inline const Modele_rayo_semi_transp& Modele() const;
  inline void associer_modele_rayo(Modele_rayo_semi_transp& m);

protected :

};


inline Modele_rayo_semi_transp& Source_rayo_semi_transp::Modele()
{
  return valeur().Modele();
}

inline const Modele_rayo_semi_transp& Source_rayo_semi_transp::Modele() const
{
  return valeur().Modele();
}

inline void Source_rayo_semi_transp::associer_modele_rayo(Modele_rayo_semi_transp& m)
{
  valeur().associer_modele_rayo(m);
}


#endif
