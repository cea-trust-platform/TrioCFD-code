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
#ifndef Postraitement_ft_lata_included
#define Postraitement_ft_lata_included

#include <Postraitement.h>
#include <vector>
#include <TRUST_Ref.h>

class Transport_Interfaces_FT_Disc;
class Motcle;
class Maillage_FT_Disc;
class Fichier_Lata;
class Comm_Group;

class Postraitement_ft_lata : public Postraitement
{
  Declare_instanciable(Postraitement_ft_lata);

public:
  void set_param(Param& param) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;

  int write_extra_mesh() override;
  void postprocess_field_values() override;

protected:
  void lire_champ_interface(Entree&);

  int ecrire_maillage_ft_disc();
  int filter_out_virtual_fa7(IntTab& new_fa7);
  void filter_out_array(const DoubleTab& dtab, DoubleTab& new_dtab) const;

  // L'equation de transport contenant les interfaces a postraiter
  OBS_PTR(Transport_Interfaces_FT_Disc) refequation_interfaces;
  // Quels champs d'interface faut-il postraiter ?
  Motcles liste_champs_i_aux_sommets;
  Motcles liste_champs_i_aux_elements;

  Nom id_domaine_;  // INTERFACES or PARTICULES - computed in ecrire_maillage_ft_disc()

  // Renumbering array for interface facettes to keep only real facettes - updated at each time step!
  std::vector<int> renum_;
  bool no_virtuals_ = false;  // whether to exclude virtual elements when writing out interface mesh and fields.
};

#endif
