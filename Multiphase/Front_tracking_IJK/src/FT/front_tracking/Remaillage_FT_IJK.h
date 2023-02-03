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
// File      : Remaillage_FT_IJK.h
// Directory : $IJK_ROOT/src/FT/front_tracking
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Remaillage_FT_IJK_included
#define Remaillage_FT_IJK_included

#include <Remaillage_FT.h>
#include <Maillage_FT_IJK.h>
/*! @brief : class Remaillage_FT_IJK
 *
 *  <Description of class Remaillage_FT_IJK>
 *
 *
 *  Pour initialiser la classe il faut:
 *  - appeler le readOn pour lire les parametres (optionnels)
 *  - appeler associer_domaine()
 *
 */
class Maillage_FT_Disc;

class Remaillage_FT_IJK : public Remaillage_FT
{

  Declare_instanciable(Remaillage_FT_IJK);

public :
  Vecteur3 get_delta_euler(const Maillage_FT_IJK& maillage) const;
  void barycentrer_lisser_systematique_ijk(Maillage_FT_Disc& maillage,
                                           ArrOfDouble& var_volume);
  void barycentrer_lisser_apres_remaillage(Maillage_FT_Disc& maillage, ArrOfDouble& var_volume);

  void remaillage_local_interface(double temps, Maillage_FT_IJK& maillage);
protected :
// N'oublie pas de mettre a jour le tableau compo_connexe_facettes
  int diviser_grandes_aretes(Maillage_FT_IJK& maillage) const;

  double calculer_longueurIdeale2_arete(const Maillage_FT_Disc& maillage, int som0, double x, double y, double z) const override;

};

#endif /* Remaillage_FT_IJK_included */
