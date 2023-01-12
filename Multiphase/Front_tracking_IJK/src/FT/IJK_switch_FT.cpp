//TRUST_NO_INDENT
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

#include <IJK_switch_FT.h>
#include <init_forcage_THI.h>
#include <corrections_qdm.h>
#include <Force_sp.h>


Implemente_instanciable_sans_constructeur(Switch_FT_double, "Switch_FT_double", Switch_double);

Switch_FT_double::Switch_FT_double() :
    old_ijk_splitting_ft_extension_(0)
{}

Sortie & Switch_FT_double::printOn(Sortie&s) const
{
  return s;
}
Entree & Switch_FT_double::readOn(Entree&s)
{
  return s;
}

void Switch_FT_double::set_param(Param& param)
{
  Switch_double::set_param(param);
  // Parametres pour le FT:
  param.ajouter("interfaces", &interfaces_);
  param.ajouter("old_ijk_splitting_ft_extension", &old_ijk_splitting_ft_extension_);
  // Parmetres pour la thermique:
  param.ajouter("thermique", &thermique_);
  // GAB : gabriel.ramirez@cea.fr
  // Parametres pour le forcage spectral 
  /* Voir reprendre probleme dans IJK_FT.cpp */
  param.ajouter("forcage", &old_forcage_);
  // Parametres pour la correction de qdm 
  param.ajouter("corrections_qdm", &old_qdm_corrections_);
  param.ajouter("reprise_vap_velocity_tmoy", &vap_velocity_tmoy_);
  param.ajouter("reprise_liq_velocity_tmoy", &liq_velocity_tmoy_);
  }

void Switch_FT_double::initialise()
{
  Cout << que_suis_je() <<"::initialise() Pb of type IJK_FT detected." << finl;

  // Probleme of type FT :
  // old_mesh_ and new_mesh_ are acutally splittings...
  const IJK_Grid_Geometry& geom = old_mesh_.get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const IJK_Grid_Geometry& geom_new = new_mesh_.get_grid_geometry();
  const double new_dx = geom_new.get_constant_delta(DIRECTION_I);
  const double new_dy = geom_new.get_constant_delta(DIRECTION_J);
  const Nom& vdf_name = geom_new.le_nom()+Nom("_VDF");

  const double old_to_new_ratio = std::min(dx/new_dx, dy/new_dy);
  const int new_ijk_splitting_ft_extension = (int)std::lrint(std::ceil(old_ijk_splitting_ft_extension_* old_to_new_ratio));
  Cerr << "Extended splitting dimensions. old = " << old_ijk_splitting_ft_extension_
      << " new = " <<new_ijk_splitting_ft_extension << finl;
  Cerr << "Construction de la zone VDF..." << finl;
  IJK_Splitting splitting_ft;
  build_extended_splitting(new_mesh_, splitting_ft, new_ijk_splitting_ft_extension);
  // Le probleme ft disc qui porte le maillage vdf pour les algorithmes front-tracking
  Probleme_base& refprobleme_ft_disc = creer_zone_vdf(splitting_ft, vdf_name);
  const Zone_dis& zone_dis = refprobleme_ft_disc.domaine_dis().zone_dis(0);
  interfaces_.initialize(splitting_ft /* splitting_FT */,
      new_mesh_ /* splitting_Ns */,
      zone_dis, 
      false);
  interfaces_.set_reprise(1);
  interfaces_.lire_maillage_ft_dans_lata();
  
  Switch_double::initialise();

  // GAB
  const int nproc_tot = Process::nproc();
  old_forcage_.compute_initial_chouippe(nproc_tot,geom,old_ni_,old_nj_,old_nk_,splitting_ft,nom_sauvegarde_);
  new_forcage_.compute_initial_chouippe(nproc_tot,geom_new,new_ni_,new_nj_,new_nk_,splitting_ft,nom_sauvegarde_);
  Cout << "new_forcage_.get_semi_gen()" << new_forcage_.get_semi_gen() << finl;
  Cout << "old_forcage_.get_semi_gen()" << old_forcage_.get_semi_gen() << finl;
  Cout << "new_forcage_.get_b_flt()" << new_forcage_.get_b_flt() << finl;
}

/////////////////////////////////////:
// MODIFICATIONS GAB : gabriel.ramirez@cea.fr
/*  Cout << "forcage_.get_type_forcage() : " <<forcage_.get_type_forcage() << finl;
  if (forcage_.get_type_forcage() > 0)
    {
      const IJK_Splitting& gbz_splitting = velocity_[0].get_splitting();
      const IJK_Grid_Geometry& my_geom = velocity_[0].get_splitting().get_grid_geometry();

      const int my_ni = velocity_[0].ni();
      const int my_nj = velocity_[0].nj();
      const int my_nk = velocity_[0].nk();
      const int nproc_tot = Process::nproc();
      Cout << "BF compute_initial_chouippe" << finl;
      Cout << "ni : "<<my_ni<<" ,nj : "<<my_nj<<" ,nk : "<<my_nk << finl;
      std::cout << "in initialise i_offset : " << gbz_splitting.get_offset_local(DIRECTION_I) << std::endl;
      std::cout << "Process::me()" << Process::me() << std::endl;
      forcage_.compute_initial_chouippe(nproc_tot,geom_new,new_ni_,new_nj_,new_nk_,splitting_ft,nom_sauvegarde_);
      statistiques().begin_count(m2);
      Cout << "AF compute_initial_chouippe" << finl;
    }
    */
/////////////////////////////////////: gabriel.ramirez@cea.fr


int Switch_FT_double::init_thermique()
// The Thermique.initialize does both the allocate and the initialize:
{
  int nb_allocated_arrays = 0;
  
  int idx =0;
    for (auto& itr : thermique_)
    {
      Cout << "Reading the old temperature field from " << Nom(curseur->get_fichier_sauvegarde())
               << " to fill the thermique_ field."<< finl;
      nb_allocated_arrays += itr.initialize(old_mesh_, idx);
      idx++;
    }
  return nb_allocated_arrays;
}

void Switch_FT_double::prepare_thermique(const Nom lata_name)
{
  for (auto& itr : thermique_)
      itr.set_fichier_sauvegarde(lata_name);
}


// flag and_lata to know if we also create the associated lata
void Switch_FT_double::ecrire_fichier_reprise(const char *fichier_sauvegarde, const bool and_lata)
{
  Nom lata_name(fichier_sauvegarde);
  lata_name += ".lata";

  if (Process::je_suis_maitre())
    {
      Cerr << "T= " << current_time_ << " Checkpointing dans le fichier " << fichier_sauvegarde << finl;
      SFichier fichier(fichier_sauvegarde);
      fichier.precision(17);
      fichier.setf(std::ios_base::scientific);
      fichier << "{\n"
          << " tinit " << current_time_ << "\n"
          << " terme_acceleration_init " << terme_source_acceleration_ << "\n"
          << " reprise_vap_velocity_tmoy " << vap_velocity_tmoy_ << "\n"
          << " reprise_liq_velocity_tmoy " << liq_velocity_tmoy_ << "\n"
          << " fichier_reprise_vitesse " << lata_name << "\n"
          << " timestep_reprise_vitesse 1\n";

      interfaces_.set_fichier_sauvegarde(lata_name);
      Cerr << "  saving interfaces... " << finl;
      fichier << " interfaces " << interfaces_  << "\n";
      fichier <<  " forcage " << new_forcage_;
      fichier <<  " corrections_qdm " << new_qdm_corrections_;
      Cerr << "  potentially saving temperature fields... " << finl;
      int flag_list_not_empty = 0;
      if ((int) thermique_.size() > 0)
        {
          fichier << " thermique {\n" ;
          flag_list_not_empty = 1;
        }
      int idx =0;
      for (auto itr = thermique_.begin(); itr != thermique_.end(); )
        {
          //curseur->sauvegarder_temperature(lata_name, idx);
          (*itr).set_fichier_sauvegarde(lata_name);
          fichier << *itr ;
          ++itr;
          if (itr)
            fichier << ", \n" ;
          else
            fichier << "\n" ;
          Cerr << "  end of temperature field #" << idx << "... " << finl;
          ++idx;
        }
      if (flag_list_not_empty)
        fichier << " } \n" ;

      fichier << "}\n";
    }

  if (and_lata)
    {
      ecrire_header_lata(lata_name);
      write_velocity(lata_name);
    }
}

void Switch_FT_double::ecrire_header_lata(const Nom lata_name) // const
{
  Switch_double::ecrire_header_lata(lata_name);

  Cout << "Adding interfaces into " << lata_name << finl;
  interfaces_.sauvegarder_interfaces(lata_name);
  // Le bloc suivant aurait pour consequence d'ecrire la temperature initiale (sur le maillage relu)
  // (car c'est pour lui qu'on a fait un allocate dans l'objet thermique pour le relire)
  // dans le fichier de sauvegarde des champs exrits pour la reprise. Ce n'est pas ce que l'on veut.
  // On avait pu se permettre cela pour l'interface car on ecrit directement celui qu'on a relu, meme
  // s'il n'a pas le bon support (le vieux geom).

}

void Switch_FT_double::set_param_reprise(Param& param)
{
  Switch_double::set_param_reprise(param);
  param.ajouter("interfaces", & interfaces_);
  param.ajouter("thermique", & thermique_);

  // GAB : gabriel.rmairez@cea.fr
  /* Voir reprendre probleme dans IJK_FT.cpp */
  param.ajouter("forcage", &new_forcage_);
  param.ajouter("corrections_qdm", &new_qdm_corrections_);
  param.ajouter("reprise_vap_velocity_tmoy", &vap_velocity_tmoy_);
  param.ajouter("reprise_liq_velocity_tmoy", &liq_velocity_tmoy_);
  // fin GAB : gabriel.ramirez@cea.fr

}

void Switch_FT_double::compute_and_write_extra_fields(const Nom& lata_name,
                                                    DoubleTab& coeff_i, IntTab Indice_i, 
                                                    DoubleTab& coeff_j, IntTab Indice_j,
                                                    DoubleTab& coeff_k, IntTab Indice_k)
{
	
	// GAB : gabriel.ramirez@cea.fr
  	/*Cout << "forcage2_.get_type_forcage() : " <<forcage2_.get_type_forcage() << finl;
  	if (forcage2_.get_type_forcage() > 0)
    {
      const IJK_Splitting& gbz_splitting = velocity_[0].get_splitting();
      const IJK_Grid_Geometry& my_geom = velocity_[0].get_splitting().get_grid_geometry();

      const int my_ni = velocity_[0].ni();
      const int my_nj = velocity_[0].nj();
      const int my_nk = velocity_[0].nk();
      const int nproc_tot = Process::nproc();
      forcage2_.compute_initial_chouippe(nproc_tot,my_geom,my_ni,my_nj,my_nk,gbz_splitting,nom_sauvegarde_);
      statistiques().begin_count(m2);
      Cout << "AF compute_initial_chouippe" << finl;
    }*/
    // GAB : gabriel.ramirez@cea.fr
  IJK_Field_double new_temperature;

  if (thermique_.size() > 0)
    {
      calculer_coords_elem();
      calculer_coeff(coeff_i,Indice_i,
                     coeff_j,Indice_j,
                     coeff_k,Indice_k);
      new_temperature.allocate(new_mesh_ /* it is in fact a splitting */, IJK_Splitting::ELEM, 0);
    }

  int idx = 0;
  for (auto& itr : thermique_)
    {
      const IJK_Field_double& old_temperature = itr.get_temperature();
      //IJK_Field_double& new_temperature = curseur->set_temperature();
      switch_scalar_field(old_temperature, new_temperature,
                          coeff_i, Indice_i,
                          coeff_j ,Indice_j,
                          coeff_k ,Indice_k);

      Cout << "Writing " <<Nom("TEMPERATURE_")+Nom(idx) << " into " << lata_name << finl;
      // Process::exit();
      //const int latastep = 0;
      //      std::ostringstream oss;
      //              oss << "TEMPERATURE_" << idx;
      //            Nom nom_temp(oss.str().c_str());
      dumplata_scalar(lata_name,Nom("TEMPERATURE_")+Nom(idx) , new_temperature, 0 /*we store a 0 */);
      //              oss.str("");
      ++idx;
    }
}

void Switch_FT_double::compute_and_write_extra_fields_direct(SFichier& file,
                                                           DoubleTab& coeff_i, IntTab Indice_i, 
                                                           DoubleTab& coeff_j, IntTab Indice_j,
                                                           DoubleTab& coeff_k, IntTab Indice_k)
{
}



