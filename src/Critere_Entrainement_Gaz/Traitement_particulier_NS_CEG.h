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
// File:        Traitement_particulier_NS_CEG.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Critere_Entrainement_Gaz/src
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Traitement_particulier_NS_CEG_included
#define Traitement_particulier_NS_CEG_included

#include <Traitement_particulier_NS_base.h>
#include <TRUST_Ref.h>

class Front_VF;

class Traitement_particulier_NS_CEG : public Traitement_particulier_NS_base
{
  Declare_base_sans_constructeur_ni_destructeur(Traitement_particulier_NS_CEG);

public :
  inline Traitement_particulier_NS_CEG() {};
  inline ~Traitement_particulier_NS_CEG() override {};

  Entree& lire(Entree& ) override;
  void preparer_calcul_particulier(void) override ;
  void post_traitement_particulier(void) override ;
  inline void en_cours_de_resolution(int , DoubleTab&, DoubleTab& ,double) override {};
  inline void sauver_stat(void) const override {};
  inline void reprendre_stat(void) override {};

protected :
  // Donnees
  Nom la_surface_libre_nom_; 			// Nom de la surface libre
  REF(Front_VF) la_surface_libre_;		// Pointe vers la frontiere
  double haspi_;				// Hauteur d'aspiration
  double C_;					// Constante AREVA
  int calculer_critere_areva_,calculer_critere_cea_jaea_; // Flags pour savoir si on calcule ou pas les criteres
  int critere_cea_jaea_normalise_;		// Normalise les criteres
  int nb_mailles_mini_;				// Nombre de mailles minimal dans le vortex
  double t_deb_;				// Debut de la recherche
  double t_fin_;				// Fin de la recherche
  double dernier_temps_;
  double dt_post_;				// Periode
  int debug_;					// Impression debuggage
  double min_critere_Q_sur_max_critere_Q_;	// Pour autoriser des critere_Q negatifs dans le cercle...
  // Methodes
  void critere_areva();
  void critere_cea_jaea();
  void imprimer(const double , const Nom& , const ArrOfDouble& , const double );
  int lpost(double temps_courant, double dt_post) const;
};

#endif
