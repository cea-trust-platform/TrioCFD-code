/****************************************************************************
* Copyright (c) 2015, CEA
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
// File:        Fichier_dx.h
// Directory:   $TRUST_ROOT/src/Front_tracking_discontinu
// Version:     /main/6
//
//////////////////////////////////////////////////////////////////////////////
#ifndef Fichier_dx_included
#define Fichier_dx_included

// .DESCRIPTION        :
//  Cette classe contient un fichier de type EcrFicPartage ou EcrFicPrive,
//  en ASCII ou en BINAIRE selon le format specifie dans le constructeur.
//  Voir la doc du constructeur.

class Fichier_dx
{
public:
  enum Mode {ERASE, APPEND};
  Fichier_dx(const Nom& basename, const Nom& extension,
             Mode mode_append,
             Postraitement_dx::Format format,
             const Comm_Group * comm_group);
  virtual ~Fichier_dx();
  virtual SFichier& get_SFichier();
  virtual const Nom& get_filename() const;
  virtual const Nom& get_filename_sans_chemin() const;
  virtual int is_master() const;
  virtual void syncfile();

  virtual Postraitement_dx::Format get_format() const;

protected:
  Nom filename_;
  Nom filename_sans_ch_;
  SFichier * fichier_;
  const Comm_Group * comm_group_;
  Postraitement_dx::Format format_;
};

// .DESCRIPTION        :
//  Specialisation du Fichier_dx pour le fichier maitre: toujours en ASCII.
//  On peut utiliser une precision differente si on veut
class Fichier_dx_maitre : public Fichier_dx
{
public:
  Fichier_dx_maitre(const Nom& basename, const Nom& extension,
                    Mode mode_append,
                    Postraitement_dx::Format format,
                    const Comm_Group * comm_group);
};

#endif
