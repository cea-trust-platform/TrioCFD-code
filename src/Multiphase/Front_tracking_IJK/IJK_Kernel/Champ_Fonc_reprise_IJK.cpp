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
// File:        Champ_Fonc_reprise_IJK.cpp
// Directory:   $TRUST_ROOT/src/Kernel/VF/Champs
// Version:     /main/29
//
//////////////////////////////////////////////////////////////////////////////

#include <Champ_Fonc_reprise_IJK.h>
#include <Probleme_base.h>
#include <EcritureLectureSpecial.h>
#include <Champ_Generique.h>
#include <Entree_complete.h>
#include <LecFicDistribueBin.h>
#include <EFichierBin.h>
#include <Domaine_VF.h>
#include <MD_Vector.h>
#include <MD_Vector_std.h>
#include <MD_Vector_tools.h>
#include <Octree_Double.h>
#include <stat_counters.h>
#include <SFichier.h>

extern void convert_to(const char *s, double& ob);
Implemente_instanciable(Champ_Fonc_reprise_IJK,"Champ_Fonc_reprise_IJK",Champ_Fonc_base);


//     printOn()
/////

Sortie& Champ_Fonc_reprise_IJK::printOn(Sortie& s) const
{
  return s << que_suis_je() << " " << le_nom();
}

//// readOn
//
Entree& Champ_Fonc_reprise_IJK::readOn(Entree& s)
{
  Cerr<<"Usage : Champ_Fonc_reprise_IJK fichier.xyz nom_pb nom_inco"<<finl;
  Nom nom_fic,nom_pb;
  Nom nom_champ_inc;

  // Lecture

  s>>nom_fic;
  s>>nom_pb>>nom_champ_inc;

  // On recupere le pb, puis ensuite on cherche le champ; on recupere le domaine_dis
  const Probleme_base& pb =ref_cast(Probleme_base,Interprete::objet(nom_pb));
  REF(Champ_base) ref_ch;

  ref_ch = pb.get_champ(Motcle(nom_champ_inc));
  if (sub_type(Champ_Inc_base,ref_ch.valeur()) && nom_champ_inc == Nom("vitesse"))
    Cerr << nom_champ_inc << " is an unknown of problem " << nom_pb << " and is 'vitesse'!" << finl;
  else
    {
      Cerr << nom_champ_inc << " is not 'vitesse' (this is the only unknown supported for now)!! " << nom_pb << finl;
      exit();
    }

  // Ouverture du fichier
  statistiques().begin_count(temporary_counter_);
  EFichier fic_rep;
  fic_rep.set_bin(1);
  fic_rep.ouvrir(nom_fic);
  if(fic_rep.fail())
    {
      Cerr<<"Error while opening the file of resumption : " <<nom_fic<<finl;
      exit();
    }

  EFichier& fich = fic_rep;

  if (fich.eof())
    {
      Cerr << "Error in Champ_Fonc_reprise_IJK::reprendre" << finl;
      Cerr << "The resumption file does not exist" << finl;
      Cerr << "or could not be opened correctly." << finl;
      exit();
    }

//  Nom ident_lu;
//  fich >> ident_lu;
//  if (ident_lu != "xyz_ijk:")
//    {
//      Cerr << "Wrong file format. Expected XYZ_IJK."<<finl;
//      exit();
//    }

  associer_domaine_dis_base(pb.domaine_dis());
  // on cree un champ comme le ch_ref;
  vrai_champ_.typer(ref_ch.valeur().que_suis_je());
  const Champ_Inc_base& ch_inc=ref_cast(Champ_Inc_base,ref_ch.valeur());
  Champ_Inc_base& v_champ=vrai_champ_.valeur();
  le_champ().associer_domaine_dis_base(pb.domaine_dis());

  v_champ.fixer_nb_valeurs_temporelles(2);
  v_champ.nommer(ch_inc.le_nom());
  v_champ.fixer_nb_comp(ch_inc.nb_comp());
  v_champ.fixer_nb_valeurs_nodales(ch_inc.nb_valeurs_nodales());

  nb_compo_ = ch_inc.nb_comp();

  // Reading file data:
  reprendre_IJK(fich, le_champ());

  statistiques().end_count(temporary_counter_);
  Cerr << "End of resuming the file " << nom_fic << " after " << statistiques().last_time(temporary_counter_) << " s" << finl;

  return s ;
}

/*! @brief Reciproque de la methode ecrit(.
 *
 * ..), lit uniquement les items sequentiels (donc pas les items communs recus d'un autre processeur)
 *    On verifie a la fin qu'on a bien lu exactement le nombre d'items attendus, s'il en manque
 *    c'est que le epsilon n'est pas bon (ou qu'on a change le maillage...)
 *  Valeur de retour: nombre total d'items sequentiels lus (sur tous les procs)
 *
 */
static int lire_special(Entree& fich, const DoubleTab& coords, DoubleTab& val, const double epsilon)
{
  const int dim = coords.dimension(1);
  const int nb_dim = val.nb_dim();
  const int nb_comp = (nb_dim == 1) ? 1 : val.dimension(1);

//  SFichier toto;
//  toto.set_bin(0);
//  toto.ouvrir("/tmp/toto.dbg");
//  coords.ecrit(toto);

  const MD_Vector& md_vect = val.get_md_vector();
  // Dans un premier temps, 1 si l'item est a lire, 0 s'il est lu par un autre processeur.
  // Une fois que l'item est lu, on met le flag a 2.
  ArrOfInt items_to_read;
  const int n_to_read = MD_Vector_tools::get_sequential_items_flags(md_vect, items_to_read);
  Octree_Double octree;
  // Build an octree with "thick" nodes (epsilon size)
  octree.build_nodes(coords, 0 /* do not include virtual elements */, epsilon);
  const ArrOfInt& floor_elements = octree.floor_elements();

  // Le fichier contient ce nombre de lignes pour cette partie du tableau (nombre total d'items sequentiels)
  const int ntot = Process::mp_sum(n_to_read);

  // On lit dans le fichier par blocs de buflines_max parce qu'il y a
  //  un broadcast reseau a chaque comm:
  const int buflines_max = 2048; // pas trop, histoire d'avoir plusieurs blocs dans les cas tests
  DoubleTab buffer(buflines_max, dim + nb_comp);
  int bufptr = buflines_max;
  ArrOfInt items;
  items.set_smart_resize(1);

  double max_epsilon_needed = epsilon;
  // Combien de fois on a trouve plusieurs candidats a moins de epsilon ?
  int error_too_many_matches = 0;
  // Combien de fois on est tombe plusieurs fois sur le meme sommet a lire ?
  int error_duplicate_read = 0;
  // Combien d'items a-t-on lu ?
  int count_items_read = 0;

  // Boucle sur les items sequentiels du fichier:
  int pourcent=0;
  for (int i = 0; i < ntot; i++)
    {
      int tmp=(i*10)/(ntot-1);
      if (tmp>pourcent || i==0)
        {
          pourcent=tmp;
          Cerr<<"\r"<<pourcent*10<<"% of data has been found."<<flush;
        }
      if (bufptr == buflines_max)
        {
          bufptr = 0;
          int n = std::min(buflines_max, ntot - i) * (dim + nb_comp);
          assert(n <= buffer.size_array());
          fich.get(buffer.addr(), n);
        }
      const double x = buffer(bufptr, 0);
      const double y = buffer(bufptr, 1);
      const double z = (dim == 3) ? buffer(bufptr, 2) : 0.;
      double eps = 1.0e-04;
      if ((std::fabs(x - (-0.2 + 0.00625*0.5)) < eps) && (std::fabs(y - (-0.2 + 0.00625*0.5)) < eps) && (std::fabs(z) < eps))
        Cerr << "tiens tiens ..." << finl;
      // Recherche des items correspondant potentiellement au point (x,y,z)
      int index = -1;
      int nb_items_proches = octree.search_elements(x, y, z, index);
      if (nb_items_proches == 0)
        {
          int a = 0;
          a++;
          Cerr << "pas ditems proches" << finl;
        }

      if (nb_items_proches > 0)
        {
          items.resize_array(nb_items_proches, ArrOfInt::NOCOPY_NOINIT);
          // Voir doc de Octree_Double::search_elements: on copie les indices des items proches dans items:
          for (int j = 0; j < nb_items_proches; j++)
            items[j] = floor_elements[index++];
          // On reduit la liste pour avoir uniquement les items a moins de epsilon
          const int item_le_plus_proche = octree.search_nodes_close_to(x, y, z, coords, items, epsilon);
          nb_items_proches = items.size_array();
          if (nb_items_proches == 1)
            {
              const int flag = items_to_read[item_le_plus_proche];
              if (flag == 1)
                {
                  // Ok, il faut lire cette valeur
                  items_to_read[item_le_plus_proche] = 2;
                  count_items_read++;
                  if (nb_dim == 1)
                    {
                      val(item_le_plus_proche) = buffer(bufptr, dim);
                    }
                  else
                    {
                      for (int j = 0; j < nb_comp; j++)
                        val(item_le_plus_proche, j) = buffer(bufptr, dim + j);
                    }
                }
              else if (flag == 0)
                {
                  // Cet item n'est pas a moi, ne pas le lire
                  int a =0;
                  a++;
                }
              else
                {
                  // Erreur, on a deja lu cet item !!! epsilon est trop grand (ou erreur a la sauvegarde ???)
                  error_duplicate_read++;
                }
            }
          else if (nb_items_proches == 0)
            {
              // ok, le sommet est sur un autre processeur (ou epsilon trop petit ??)
              int a =0;
              a++;
              Cerr << "pas ditems proches 2" << finl;
            }
          else
            {
              // Erreur: epsilon est trop grand, on a plusieurs candidats a moins de epsilon
              // Calcul de la distance avec le deuxieme plus proche pour afficher un message d'erreur a la fin:
              for (int ii = 0; ii < nb_items_proches; ii++)
                {
                  const int i_coord = items[ii];
                  if (i_coord == item_le_plus_proche)
                    continue; // celui-la est sans doute le bon, il faut un epsilon superieur a cette valeur la...
                  double xx = 0;
                  for (int j = 0; j < dim; j++)
                    {
                      double yy = coords(i_coord, j) - buffer(bufptr, j);
                      xx += yy * yy;
                    }
                  // On propose de mettre un epsilon au maximum egal a 1/10 de la distance avec le deuxieme point le plus proche:
                  xx = sqrt(xx) * 0.1;
                  if (max_epsilon_needed > xx)
                    max_epsilon_needed = xx;
                }
              error_too_many_matches++;
            }
        }
      bufptr++;
    }
  Cerr << finl;
  // Erreurs ?
  int err = (count_items_read != n_to_read) || (error_too_many_matches > 0) || (error_duplicate_read > 0);
  err = Process::mp_sum(err);
  if (err)
    {
      error_too_many_matches = Process::mp_sum(error_too_many_matches);
      error_duplicate_read = Process::mp_sum(error_duplicate_read);
      max_epsilon_needed = Process::mp_min(max_epsilon_needed);
      if (Process::je_suis_maitre())
        {
          if (error_too_many_matches)
            {
              Cerr << "Error in EcritureLectureSpecial: error_too_many_matches = " << error_too_many_matches
                   << ", epsilon is too large. Suggested value: " << max_epsilon_needed << finl;
              if (max_epsilon_needed==0)
                {
                  Cerr << "It could be because your mesh has two boundaries which match exactly." << finl;
                  Cerr << "It is possible to do calculation with this property but xyz restart process" << finl;
                  Cerr << "is impossible because it can't detect the differences between faces of the two boundaries..." << finl;
                  Cerr << "Try to do a classic restart with .sauv files." << finl;
                }
            }
          else if (error_duplicate_read)
            {
              Cerr << "Error in EcritureLectureSpecial: error_duplicate_read = " << error_duplicate_read
                   << ", probably epsilon too large. " << finl;
            }
          else
            {
              Cerr << "Error in EcritureLectureSpecial: Some items were not found: epsilon too small (or changed the mesh ?)" << finl;
            }
        }
      Process::barrier();
      Process::exit();
    }
  return ntot;
}


void Champ_Fonc_reprise_IJK::reprendre_IJK(Entree& fich, Champ_base& ch)
{
  int nb_val_nodales_old = nb_valeurs_nodales();

  const Domaine_VF& zvf=ref_cast(Domaine_VF,ch.domaine_dis_base());
  DoubleTab& val = ch.valeurs();

  /* [ABN] ce qui suit est une version allegee de
         static int lecture_special_part2(const Domaine_VF& zvf, Entree& fich, DoubleTab& val)
     en depilant tous les appels sous-jacents.
   */

  const MD_Vector& md = val.get_md_vector();
  int ntot = 0;

  if (!md.non_nul())
    {
      Cerr << "Champ_Fonc_reprise_IJK::reprendre_IJK: error, cannot read intno an array with no metadata" << finl;
      Process::exit();
    }
  const int nb_items_seq = md.valeur().nb_items_seq_tot();
  if (nb_items_seq == 0)
    return;

  if (sub_type(MD_Vector_std, md.valeur()))
    {

      if (md != zvf.face_sommets().get_md_vector())
        {
          Cerr << "Champ_Fonc_reprise_IJK::reprendre_IJK: Problem! Expected face descriptor."<< finl;
          exit();
        }
      const DoubleTab& coords = zvf.xv(); // Descripteur des face
      const double epsilon = zvf.domaine().epsilon();
      ntot += lire_special(fich, coords, val, epsilon);
    }
  else
    {
      Cerr << "Champ_Fonc_reprise_IJK::reprendre_IJK: Problem in the resumption. Not a MD_vector_std." << finl;
      exit();
    }

  if (ntot != nb_items_seq)
    {
      Cerr << "Champ_Fonc_reprise_IJK::reprendre_IJK: Internal error" << finl;
      exit();
    }


  if (nb_val_nodales_old != nb_valeurs_nodales())
    {
      Cerr << "Champ_Fonc_reprise_IJK::reprendre_IJK: Problem in the resumption "<< finl;
      Cerr << "The field wich is read, does not have same number of nodal values" << finl;
      Cerr << "that the field created by the discretization " << finl;
      exit();
    }
  Cerr << "Resume of the field " << nom_ << " performed." << finl;
}

void Champ_Fonc_reprise_IJK::mettre_a_jour(double t)
{
  Champ_Fonc_base::mettre_a_jour(t);
}

