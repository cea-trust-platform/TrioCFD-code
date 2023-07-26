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
// File      : Statistiques_dns_ijk_FT.h
// Directory : $IJK_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////
#ifndef Statistiques_dns_ijk_FT_H
#define Statistiques_dns_ijk_FT_H
#include <FixedVector.h>
#include <IJK_Field.h>
#include <Statistiques_dns_ijk.h>
#include <TRUSTArrays.h>
#include <TRUST_Ref.h>

class IJK_FT_double;
class IJK_Grid_Geometry;

class Statistiques_dns_ijk_FT : public Statistiques_dns_ijk
{
  Declare_instanciable_sans_constructeur(Statistiques_dns_ijk_FT);
public:
  Statistiques_dns_ijk_FT(); // Je ne sais pas compiler le Vect(Statistiques_dns_ijk_FT) sans lui...
  Statistiques_dns_ijk_FT(IJK_FT_double& ijk_ft);
  using Statistiques_dns_ijk::initialize;
  void initialize(const IJK_FT_double& ijk_ft,const IJK_Grid_Geometry&);
  int initialize(const IJK_FT_double& ijk_ft,const IJK_Splitting& splitting,
                 const int check_stats);
  Sortie& completer_print(Sortie& os) const override;
  void completer_read(Param& param) override;
  int lire_motcle_non_standard(const Motcle& mot, Entree& is) override;

  void update_stat(IJK_FT_double& cas, const double dt);

  void postraiter(Sortie&, int flag_valeur_instantanee = 0) const;
  void postraiter_thermique(const double t) const;
  double compute_desequil_alpha(const IJK_Grid_Geometry& geom_NS,
                                const double portee_wall_repulsion) const;

  const FixedVector<IJK_Field_double, 3>& get_IJK_vector_field(const Nom& nom) const;

  // Pour les deriv de U, V et W :
  FixedVector<IJK_Field_double, 3> gradU_;
  FixedVector<IJK_Field_double, 3> gradV_;
  FixedVector<IJK_Field_double, 3> gradW_;

  // Pour les deriv secondes de P, U, V et W :
  FixedVector<IJK_Field_double, 3> grad2Pi_; // Partie diagonale de la jacobienne
  FixedVector<IJK_Field_double, 3> grad2Pc_; // contient les deriv croisees
  FixedVector<IJK_Field_double, 3> grad2Ui_; // Partie diagonale de la jacobienne
  FixedVector<IJK_Field_double, 3> grad2Uc_; // contient les deriv croisees
  FixedVector<IJK_Field_double, 3> grad2Vi_; // Partie diagonale de la jacobienne
  FixedVector<IJK_Field_double, 3> grad2Vc_; // contient les deriv croisees
  FixedVector<IJK_Field_double, 3> grad2Wi_; // Partie diagonale de la jacobienne
  FixedVector<IJK_Field_double, 3> grad2Wc_; // contient les deriv croisees
protected:
  int check_stats_;
  int nb_thermal_fields_; // Number of objects thermique_ in the list
  int nvalt_;             // Number of variables post-processed per field
  REF(IJK_FT_double) ref_ijk_ft_;
  // Last instantaneous value of the space average (only on processor 0)
  DoubleTab moyenne_spatiale_instantanee_temperature_; // (i,j,k) ->  (nvalt_,nb_elem_k_tot,nb_thermal_fields_)
  // Temporal integral of statistics variables
  DoubleTab integrale_temporelle_temperature_; // (i,j,k) ->  (nvalt_,nb_elem_k_tot,nb_thermal_fields_)
  VECT(Nom) noms_moyennes_temperature_;
};
#endif
