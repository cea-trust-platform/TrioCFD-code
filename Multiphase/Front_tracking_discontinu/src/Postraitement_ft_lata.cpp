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
// File:        Postraitement_ft_lata.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/23
//
//////////////////////////////////////////////////////////////////////////////
#include <Postraitement_ft_lata.h>
#include <Probleme_FT_Disc_gen.h>
#include <Transport_Marqueur_FT.h>
#include <Motcle.h>
#include <Domaine.h>
#include <EcrFicPartage.h>
#include <EcrFicPartageBin.h>
#include <Fichier_Lata.h>
#include <TRUSTTab.h>
#include <communications.h>
#include <SFichier.h>
#include <Param.h>

Implemente_instanciable_sans_constructeur(Postraitement_ft_lata,"Postraitement_ft_lata",Postraitement_lata);

// Description: Constructeur par defaut
Postraitement_ft_lata::Postraitement_ft_lata()
{
}

// Description: Appel a Postraitement_lata::readOn
// Precondition: mon_probleme.valeur() doit etre de type :
// Probleme_FT_Disc_gen
// Pb_Thermohydraulique_Especes_QC
// Pb_Thermohydraulique_Especes_Turbulent_QC
Entree& Postraitement_ft_lata::readOn(Entree& is)
{
  // Verification du type du probleme
  const Nom& type_pb =  mon_probleme.valeur().que_suis_je();
  if ( (!sub_type(Probleme_FT_Disc_gen, mon_probleme.valeur()))
       && (type_pb!="Pb_Thermohydraulique_Especes_QC")
       && (type_pb!="Pb_Thermohydraulique_Especes_Turbulent_QC") )
    {

      Cerr << " Reading Postraitement_ft_lata\n";
      Cerr << " postraitement_ft_lata is not accepted for a problem of type "<<type_pb << finl;
      Cerr << " The recognized problems are :" << finl;
      Cerr << " Probleme_FT_Disc_gen " << finl;
      Cerr << " Pb_Thermohydraulique_Especes_QC " << finl;
      Cerr << " Pb_Thermohydraulique_Especes_Turbulent_QC " << finl;
      exit();
    }

  Postraitement_lata::readOn(is);
  return is;
}

Sortie& Postraitement_ft_lata::printOn(Sortie& os) const
{
  return os;
}

void Postraitement_ft_lata::set_param(Param& param)
{
  Postraitement_lata::set_param(param);
  param.ajouter_non_std("interfaces",(this));
}

int Postraitement_ft_lata::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  Motcle motlu;
  if (mot=="interfaces")
    {
      is >> motlu;
      if (refequation_interfaces.non_nul())
        {
          Cerr<<" Only one transport interface equation name can be specified "<<finl;
          Cerr<<" for a Postraitement_ft_lata post-process."<<finl;
          Cerr<<" The "<<Motcle(refequation_interfaces->le_nom())<<" has already been readen for the post-process "<<(*this).le_nom()<<finl;
          Cerr<<" Please, create a new Postraitement_ft_lata post-process"<<finl;
          Cerr<<" for the "<<motlu<<" transport interface equation."<<finl;
          exit();
        }

      if (Process::je_suis_maitre())
        Cerr << "Post-processing for the interface of transport equation : " << motlu << finl;
      if (sub_type(Probleme_FT_Disc_gen, mon_probleme.valeur()))
        {
          const Probleme_FT_Disc_gen& pb =
            ref_cast(Probleme_FT_Disc_gen, mon_probleme.valeur());
          refequation_interfaces = pb.equation_interfaces(motlu);
        }
      else
        for (int i=0; i<mon_probleme->nombre_d_equations(); i++)
          {
            const Nom& nom_eq =  mon_probleme->equation(i).le_nom();
            if (sub_type(Transport_Interfaces_FT_Disc,mon_probleme->equation(i))
                && (motlu==Motcle(nom_eq)))
              {
                refequation_interfaces=ref_cast(Transport_Interfaces_FT_Disc,mon_probleme->equation(i));
              }
          }

      if (!refequation_interfaces.non_nul())
        {
          Cerr<<" No interface equation name  "<<motlu<<" has been found. "<<finl;
          exit();
        }
      is >> motlu;
      if (motlu != "{")
        {
          Cerr << " Postraitement_ft_lata::lire_champ\n";
          Cerr << " { was expected after the keyword interfaces\n";
          Cerr << " It has been found " << motlu << finl;
          exit();
        }

      while (1)
        {
          is >> motlu;
          if (motlu == "champs")
            lire_champ_interface(is);
          else if (motlu == "}")
            break;
          else
            {
              Cerr << "Postraitement_ft_lata::lire_champ\n";
              Cerr << "The keyword champs was expected.";
              Cerr << "It has been found " << motlu << finl;
              exit();
            }
        }
      return 1;
    }
  else
    return Postraitement_lata::lire_motcle_non_standard(mot,is);
  return 1;
}

// Description: lecture de la liste de champs aux interfaces a postraiter
void Postraitement_ft_lata::lire_champ_interface(Entree& is)
{
  Motcle motcle;
  is >> motcle;
  Localisation loc = SOMMETS;
  if (motcle == "sommets")
    loc = SOMMETS;
  else if (motcle == "elements")
    loc = ELEMENTS;
  else
    {
      Cerr << "Error for Postraitement_ft_lata::lire_champ_interface :\n";
      Cerr << motcle <<" has been readen. "<< finl;
      Cerr << " Kewords sommets or elements were expected after the keyword champs " << finl;
      exit();
    }
  Motcles& liste =
    (loc == SOMMETS) ? liste_champs_i_aux_sommets : liste_champs_i_aux_elements;

  is >> motcle;
  if (motcle != "{")
    {

      Cerr << "Error for Postraitement_ft_lata::lire_champ_interface  :\n";
      Cerr << motcle << " has been readen. "<<finl;
      Cerr << " { was expected after sommets/elements." << finl;

      exit();
    }
  const Transport_Interfaces_FT_Disc& eq_interfaces = refequation_interfaces.valeur();

  while (1)
    {
      is >> motcle;
      if (motcle == "}")
        break;
      if (! eq_interfaces.get_champ_post_FT(motcle, loc, (FloatTab*) 0))
        {
          if (! eq_interfaces.get_champ_post_FT(motcle, loc, (IntTab*) 0))
            {
              Cerr << "Error for Postraitement_ft_lata::lire_champ_interface :\n";
              Cerr << " The field " << motcle;
              Cerr << " is not understood for the " << (eq_interfaces.que_suis_je()=="Transport_Marqueur_FT"?"particules":"interfaces") << " or not authorized at ";
              Nom tmp = ((loc == SOMMETS) ? "sommets" : "elements");
              Cerr << tmp <<finl;
              eq_interfaces.get_champ_post_FT(demande_description, loc, (FloatTab*) 0);
              eq_interfaces.get_champ_post_FT(demande_description, loc, (IntTab*) 0);
              exit();
            }
        }
      if (!liste.contient_(motcle))
        liste.add(motcle);
    }
}

// Description:
//  Appel a Postraitement_lata (on ecrit l'entete la premiere fois et
//   on postraite les champs euleriens),
//  puis ecriture des interfaces et des champs aux interfaces.
//  Voir aussi Postraitement_base
void Postraitement_ft_lata::postraiter(int forcer)
{
  const double temps_dernier_post_before = temps_dernier_post_;
  Postraitement_lata::postraiter(forcer);

  Format_Post_Lata::Options_Para option;
  if (fichiers_multiples_ == 0)
    option = Format_Post_Lata::SINGLE_FILE;
  else
    option = Format_Post_Lata::MULTIPLE_FILES;

  Format_Post_Lata::Format fmt;
  if (format_ == ASCII)
    fmt = Format_Post_Lata::ASCII;
  else
    fmt = Format_Post_Lata::BINAIRE;

  if (! forcer)
    if (! dt_post_ecoule(dt_post_, temps_dernier_post_before))
      // L'intervalle de temps entre postraitements n'est pas ecoule
      return;
  if (Process::je_suis_maitre())
    Cerr << "Postraitement_ft_lata::postraiter time =" << temps_ << finl;

  // Pour eviter un plantage si l'utilisateur specifie un Postraitement_ft_lata
  // sans interfaces specifiees
  if (!refequation_interfaces.non_nul())
    return;

  const Maillage_FT_Disc& maillage =
    refequation_interfaces.valeur().maillage_interface_pour_post();

  // Determination du type de maillage en fonction du type de l equation
  Nom type_maillage;
  if (sub_type(Transport_Marqueur_FT,refequation_interfaces.valeur()))
    type_maillage = "PARTICULES";
  else if (sub_type(Transport_Interfaces_FT_Disc,refequation_interfaces.valeur()))
    type_maillage = "INTERFACES";
  else
    {
      Cerr<<"Type "<<refequation_interfaces.valeur().que_suis_je()<<" not recognized"<<finl;
      exit();
    }

  // Ouverture du fichier lata
  Nom filename_interfaces;
  {
    char str_temps[100] = "0.0";
    if (temps_ >= 0.)
      sprintf(str_temps, "%.10f", temps_);

    Nom basename_interf(lata_basename_);
    Nom extension_interf(".lata.");
    extension_interf+=type_maillage+".";
    const Probleme_base& pb = mon_probleme.valeur();
    extension_interf += pb.domaine().le_nom();
    extension_interf += ".";
    extension_interf += pb.le_nom();
    extension_interf += ".";
    extension_interf += str_temps;

    Fichier_Lata fichier_interf(basename_interf, extension_interf,
                                Fichier_Lata::ERASE, fmt, option);

    filename_interfaces = fichier_interf.get_filename();

    ecrire_maillage(fichier_interf, maillage, fichiers_multiples_);
  }

  Fichier_Lata_maitre fichier(lata_basename_, extension_lata(),
                              Fichier_Lata::APPEND, option);

  if (fichier.is_master())
    {
      fichier.get_SFichier() << "Champ " << type_maillage << " "
                             << remove_path(filename_interfaces)
                             << finl;
    }
  fichier.syncfile();
}

// Format lata pour le maillage :
// Si fichiers_multiples_==1, on sort un fichier par processeur, avec la numerotation locale,
// sinon:
//  nb_sommets dimension
//                (nombre total de sommets sur l'ensemble des procs)
//  x y (z)
//  x y (z)       tous les sommets du proc 1, puis du proc 2, ...
//  ...
//  nb_faces dimension         nombre total de faces
//  sommet1 sommet2 (sommet3)
//  sommet1 sommet2 (sommet3)  toutes les faces...
//  ...
//     (les sommets sont numerotes a partir de zero,
//      numerotation globale dans l'ordre d'apparition des sommets)
//  som 1 drapeau
//  drapeau
//  drapeau         les drapeaux des sommets
//  ...
//  som 1 num_pe
//  pe
//  pe              les numeros des pe de chaque sommet
//  ...              (0,0,...,0,1,1,  ...  ,nproc-1)

void Postraitement_ft_lata::ecrire_maillage(Fichier_Lata& file,
                                            const Maillage_FT_Disc& mesh,
                                            int fichiers_multiples) const
{
  SFichier& sfichier = file.get_SFichier();

  const DoubleTab& sommets = mesh.sommets();
  const IntTab& facettes = mesh.facettes();

  // Quelques calculs impliquant des communications
  const int nb_sommets_local = sommets.dimension(0);
  const int nb_sommets_total =
    fichiers_multiples ? nb_sommets_local : mp_sum(nb_sommets_local);
  const int nb_facettes_local = facettes.dimension(0);
  const int nb_facettes_total =
    fichiers_multiples ? nb_facettes_local : mp_sum(nb_facettes_local);
  const int offset_sommets =
    fichiers_multiples ? 0 : mppartial_sum(nb_sommets_local);

  // Ecriture des sommets sur tous les processeurs
  if (file.is_master())
    sfichier << sommets.dimension(1) << space << nb_sommets_total << finl;
  {
    float data[3];
    const int n = sommets.dimension(0);
    const int m = sommets.dimension(1);
    for (int i = 0; i < n; i++)
      {
        for (int j = 0; j < m; j++)
          data[j] = (float)sommets(i,j); // downcast to float ...
        sfichier.put(data, m, m+1);
      }
  }
  file.syncfile();

  // Ecriture des facettes (on ajoute un offset aux numeros de sommets
  // pour obtenir la numerotation globale: offset = somme des nombres
  // de noeuds sur les processeurs precedents)
  if (file.is_master())
    sfichier << facettes.dimension(1) << space << nb_facettes_total << finl;

  {
    IntTab facettes_renum(facettes);
    facettes_renum += offset_sommets;
    sfichier.put(facettes_renum.addr(),
                 facettes_renum.size_array(), facettes_renum.dimension(1));
  }
  file.syncfile();

  // Ecriture des champs
  const Transport_Interfaces_FT_Disc& eq_interfaces =
    refequation_interfaces.valeur();

  for (int num_loc = 0; num_loc < 2; num_loc++)
    {
      Localisation loc = SOMMETS;
      switch(num_loc)
        {
        case 0:
          loc=SOMMETS;
          break;
        case 1:
          loc=ELEMENTS;
          break;
        }

      const Motcles& liste = (loc == SOMMETS
                              ? liste_champs_i_aux_sommets
                              : liste_champs_i_aux_elements);

      const char * som_elem = (loc == SOMMETS ? "SOM" : "ELEM");

      const int nb_champs = liste.size();
      if (nb_champs == 0)
        return;

      FloatTab ftab;
      IntTab itab;

      int i;
      for (i = 0; i < nb_champs; i++)
        {
          const Motcle& nom_du_champ = liste[i];

          if (eq_interfaces.get_champ_post_FT(nom_du_champ, loc, &ftab))
            {
              // ok, le champ est dans ftab
            }
          else if (eq_interfaces.get_champ_post_FT(nom_du_champ, loc, &itab))
            {
              // copie du champ dans ftab
              if (itab.nb_dim() == 1)
                {
                  const int n = itab.dimension(0);
                  ftab.resize(n,1);
                  for (int j = 0; j < n; j++)
                    ftab(j,0) = (float)itab(j);
                }
              else
                {
                  const int n = itab.dimension(0);
                  const int m = itab.dimension(1);
                  ftab.resize(n, m);
                  for (int j = 0; j < n; j++)
                    for (int k = 0; k < m; k++)
                      ftab(j,k) = (float)itab(j,k);
                }
            }
          else
            {
              // Normalement on aurait deja du detecter cette erreur a la lecture
              // du postraitement :

              Cerr << "Error for Postraitement_ft_lata::ecrire_maillage" << finl;
              Cerr << "Unknown field : " << nom_du_champ << finl;
              exit();
            }

          const int nb_noeuds = ftab.dimension(0);
          const int nb_compo = (ftab.nb_dim() == 1 ? 1 : ftab.dimension(1));
          const int nb_noeuds_attendus =
            (loc == SOMMETS ? nb_sommets_local : nb_facettes_local);
          if (nb_noeuds != nb_noeuds_attendus)
            {
              Cerr << "Error for Postraitement_ft_lata::ecrire_maillage" << finl;
              Cerr << " nb_noeuds = " << nb_noeuds << finl;
              Cerr << " nb_noeuds expected = " << nb_noeuds_attendus << finl;
              exit();
            }

          if (file.is_master())
            {
              sfichier << som_elem << finl;
              sfichier << nb_compo << finl;
              sfichier << nom_du_champ << finl;
            }
          sfichier.put(ftab.addr(), nb_noeuds * nb_compo, nb_compo);
          file.syncfile();

        }
    }
}

