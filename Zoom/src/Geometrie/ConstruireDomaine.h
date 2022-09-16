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
// File:        ConstruireDomaine.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Geometrie
// Version:     /main/10
//
//////////////////////////////////////////////////////////////////////////////

#ifndef ConstruireDomaine_included
#define ConstruireDomaine_included


#include <Interprete_geometrique_base.h>
class Domaine;
class Sous_Zone;
#include <TRUSTTabs_forward.h>
class Zone;
#include <TRUSTTabs_forward.h>
/*! @brief Classe ConstruireDomaine Cette classe est un interprete qui sert a lire et executer
 *
 *     la directive ConstruireDomaine :
 *         ConstruireDomaine nom_domaine nom_sszone
 *     Cette directive sert a construire un nouveau domaine a partir
 *     d 'une sous zone
 *
 * @sa Interprete Rectangle
 */
class ConstruireDomaine : public Interprete_geometrique_base
{
  Declare_instanciable(ConstruireDomaine);

public :

  Entree& interpreter_(Entree&) override;

protected :

  static void creer_sommet_et_elem(Domaine&, Sous_Zone&, IntTab&);
  static void creer_bords(Zone&, Sous_Zone&, IntTab&);
  static int ajouterSom(Domaine& , const DoubleTab&, DoubleTab&, int& );

};
#endif
