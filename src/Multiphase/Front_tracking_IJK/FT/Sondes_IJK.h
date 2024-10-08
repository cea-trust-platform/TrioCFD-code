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
// File      : Sondes_IJK.h
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Sondes_IJK_included
#define Sondes_IJK_included

#include <TRUSTList.h>
#include <Sonde_IJK.h>
#include <TRUST_List.h>


/*! @brief : class Sondes_IJK
 *
 *  <Description of class Sondes_IJK>
 *
 *
 *
 */
class Sondes_IJK : public LIST(Sonde_IJK)
{

  Declare_instanciable( Sondes_IJK ) ;

public :
  inline void ouvrir_fichiers();
  inline void fermer_fichiers();
  inline void completer_IJK(const IJK_FT_double& ijk_ft);
  //  void associer_post(const Postraitement&);
  void postraiter();
  void mettre_a_jour(double temps, double tinit);

protected :
  OBS_PTR(IJK_FT_double) ref_ijk_ft_;

};

inline void Sondes_IJK::completer_IJK(const IJK_FT_double& ijk_ft)
{
  for(auto& itr : *this) itr.completer_IJK(ijk_ft);
}

inline void Sondes_IJK::ouvrir_fichiers()
{
  for(auto& itr : *this) itr.ouvrir_fichier();
}
inline void Sondes_IJK::fermer_fichiers()
{
  for(auto& itr : *this) itr.fermer_fichier();
}


#endif /* Sondes_IJK_included */
