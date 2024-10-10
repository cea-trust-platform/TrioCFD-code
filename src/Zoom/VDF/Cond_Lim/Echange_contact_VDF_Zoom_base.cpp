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
// File:        Echange_contact_VDF_Zoom_base.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/VDF/Cond_Lim
// Version:     /main/14
//
//////////////////////////////////////////////////////////////////////////////

#include <Echange_contact_VDF_Zoom_base.h>
#include <Motcle.h>

Implemente_base(Echange_contact_VDF_Zoom_base,"Paroi_Echange_contact_VDF_Zoom_base",Echange_externe_impose);


Sortie& Echange_contact_VDF_Zoom_base::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

Entree& Echange_contact_VDF_Zoom_base::readOn(Entree& s )
{
  Motcle motlu;
  Motcles les_motcles(2);
  {
    les_motcles[0] = "h_imp";
    les_motcles[1] = "T_ext";
  }

  int ind = 0;
  while (ind < 2)
    {
      s >> motlu;
      int rang = les_motcles.search(motlu);

      switch(rang)
        {
        case 0:
          {
            s >> h_imp_;
            h_paroi = h_imp_->valeurs()(0,0);
            break;
          }
        case 1:
          {
            s >> le_champ_front;
            break;
          }
        default:
          {
            Cerr << "Erreur a la lecture de la condition aux limites de type Echange_impose " << finl;
            Cerr << "On attendait " << les_motcles << "a la place de " <<  motlu << finl;
            exit();
          }
        }

      ind++;

    }

  return s ;
}

