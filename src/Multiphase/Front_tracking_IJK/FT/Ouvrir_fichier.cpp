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
// File      : Ouvrir_fichier.cpp
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////
#include <Ouvrir_fichier.h>

// reset : mode d'ecriture en append ou en out.
// suffix : fichier nomme nom_du_cas()+suffix
// msg   : message a imprimer
// head  : Si (reset=1), ios::out, on cree une entete : "# "+message+\n
SFichier Ouvrir_fichier(const Nom& suffix,
                        const Nom& head,
                        const int reset,
                        const int prec)
{
  Nom s( Objet_U::nom_du_cas());
  s+= suffix;
  IOS_OPEN_MODE mode = (reset) ? ios::out : ios::app;
  SFichier fic;
  fic.ouvrir(s,mode);
  if (reset)
    fic << "# " << head << finl;
  fic.setf(ios::scientific);
  fic.precision(prec);
  //fic << msg;
  return fic;
}

SFichier Open_file_folder(const Nom& folder,
                          const Nom& suffix,
                          const Nom& head,
                          const int reset,
                          const int prec)
{

  Nom s(folder);
  s+= "/";
  s+= Nom(Objet_U::nom_du_cas());
  s+= suffix;
  IOS_OPEN_MODE mode = (reset) ? ios::out : ios::app;
  SFichier fic;
  fic.ouvrir(s,mode);
  if (reset)
    fic << "# " << head << finl;
  fic.setf(ios::scientific);
  fic.precision(prec);
  //fic << msg;
  return fic;
}
