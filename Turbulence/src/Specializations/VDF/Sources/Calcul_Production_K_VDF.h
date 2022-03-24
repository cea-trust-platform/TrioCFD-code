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
// File:        Calcul_Production_K_VDF.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Sources
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Calcul_Production_K_VDF_included
#define Calcul_Production_K_VDF_included

/////////////////////////////////////////////////////////////////////////
//
//  CLASS Calcul_Production_K_VDF
//
// Classe qui porte les fonctions de calcul des termes de production
// de l'energie cinetique turbulente et de destruction de cette energie
// Cette classe ne derive d'aucune autre car nous voulons l'utiliser
// pour faire de l'heritage multiple.
//
/////////////////////////////////////////////////////////////////////////

#include <Nom.h>

class DoubleVect;
class DoubleTab;
class Zone_Cl_VDF;
class Champ_Face;
class Zone_VDF;

class Calcul_Production_K_VDF
{

protected:

  Calcul_Production_K_VDF() {}

  DoubleVect& calculer_terme_production_K(const Zone_VDF& ,const Zone_Cl_VDF&, DoubleVect& , const DoubleTab& ,
                                          const DoubleTab& , const Champ_Face& , const DoubleTab& ) const;

  DoubleVect& calculer_terme_production_K_Axi(const Zone_VDF& , const Champ_Face& , DoubleVect& ,
                                              const DoubleTab& , const DoubleTab& ) const;

  DoubleVect& calculer_terme_production_K_BiK(const Zone_VDF& ,const Zone_Cl_VDF&, DoubleVect& , const DoubleTab& ,
                                              const DoubleTab& , const Champ_Face& , const DoubleTab& ) const;

  DoubleVect& calculer_terme_production_K_BiK_Axi(const Zone_VDF& , const Champ_Face& , DoubleVect& ,
                                                  const DoubleTab& , const DoubleTab& ) const;

  DoubleVect& calculer_terme_destruction_K(const Zone_VDF& , const Zone_Cl_VDF& , DoubleVect& ,
                                           const DoubleTab& , const DoubleTab& ,
                                           const DoubleTab& ,const DoubleVect&   ) const;

  DoubleVect& calculer_terme_destruction_K(const Zone_VDF& , const Zone_Cl_VDF& , DoubleVect& ,
                                           const DoubleTab& , const DoubleTab& ,
                                           double ,const DoubleVect&  ) const;

  DoubleVect& calculer_terme_destruction_K(const Zone_VDF& , const Zone_Cl_VDF& , DoubleVect& ,
                                           const DoubleTab& , const DoubleTab& ,
                                           const DoubleVect& ,const DoubleVect&  ,int ) const;
  void mettre_a_jour(double ) { }

protected:


  DoubleTab& calculer_u_teta(const Zone_VDF& , const Zone_Cl_VDF& ,
                             const DoubleTab& , const DoubleTab& ,DoubleTab& ) const;

  DoubleTab& calculer_u_conc(const Zone_VDF& , const Zone_Cl_VDF& ,
                             const DoubleTab& , const DoubleTab& ,DoubleTab& ,int ) const;
};

inline void error_keps(const Nom& source, const Nom& nom)
{
  Cerr << "Error ! You can't use the " << source << " source term with a " << nom << " problem/medium !!" << finl;
  Cerr << "Check the reference manual. It is may be another source term which should be used." << finl;
  Process::exit();
}

template <typename RETURN_TYPE>
RETURN_TYPE not_implemented(const char * nom_funct)
{
  Cerr << "The method " << nom_funct << " should be implemented in a derived class !" << finl;
  throw;
}

#endif /* Calcul_Production_K_VDF_included */
