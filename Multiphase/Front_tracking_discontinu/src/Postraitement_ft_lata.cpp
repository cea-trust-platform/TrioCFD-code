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

#include <Postraitement_ft_lata.h>
#include <Probleme_FT_Disc_gen.h>
#include <Transport_Marqueur_FT.h>
#include <Motcle.h>
#include <Zone.h>
#include <EcrFicPartage.h>
#include <EcrFicPartageBin.h>
#include <Fichier_Lata.h>
#include <TRUSTTab.h>
#include <communications.h>
#include <SFichier.h>
#include <Param.h>
#include <Format_Post_Lata.h>

Implemente_instanciable(Postraitement_ft_lata,"Postraitement_ft_lata",Postraitement);

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

  Postraitement::readOn(is);
  if (!champs_demande_)
    {
      Cerr << "*********************************************************************" << finl;
      Cerr << "Warning: in Postraitement_ft_lata block, you specified interfaces to post-process" << finl;
      Cerr << "without specifying fields. Interfaces will not be post-processed unless you post-process a field also." << finl;
      Cerr << "Contact TRUST/TrioCFD support team or look for examples in TrioCFD databases" << finl;
      Cerr << "*********************************************************************" << finl;
    }

  if (!sub_type(Format_Post_Lata, format_post.valeur()))
    Process::exit("ERROR: In Postraitement_ft_lata, only the LATA (V2) format is supported! Use directive 'format lata'.");

  return is;
}

Sortie& Postraitement_ft_lata::printOn(Sortie& os) const
{
  return os;
}

void Postraitement_ft_lata::set_param(Param& param)
{
  Postraitement::set_param(param);
  param.ajouter_non_std("interfaces",(this));
}

int Postraitement_ft_lata::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  int ret_val = Postraitement::lire_motcle_non_standard(mot, is);
  if (ret_val != -1) // all good we've hit standard postprocessing keywords
    return ret_val;

  Motcle motlu;
  if (mot=="interfaces")
    {
      is >> motlu;
      if (refequation_interfaces.non_nul())
        {
          Cerr<<" Only one transport interface equation name can be specified "<<finl;
          Cerr<<" for a Postraitement_ft_lata post-process."<<finl;
          Cerr<<" The "<<Motcle(refequation_interfaces->le_nom())<<" has already been read for the post-process "<<(*this).le_nom()<<finl;
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
      if (motlu == "no_virtuals")
        {
          if (Process::nproc() > 1)
            no_virtuals_ = true;
          is >> motlu;
        }
      if (motlu != "{")
        {
          Cerr << " Postraitement_ft_lata::lire_champ\n";
          Cerr << " { was expected after the keyword interfaces\n";
          Cerr << " It has been found " << motlu << finl;
          exit();
        }

      lire_champ_interface(is);
      return 1;
    }
  else
    return -11;
}

/*! @brief lecture de la liste de champs aux interfaces a postraiter
 */
void Postraitement_ft_lata::lire_champ_interface(Entree& is)
{
  Motcle nom_champ, loc_lu;
  const Transport_Interfaces_FT_Disc& eq_interfaces = refequation_interfaces.valeur();

  while (1)
    {
      is >> nom_champ;
      if (nom_champ == "}")  break;

      is >> loc_lu;
      Localisation loc = SOMMETS;
      if (loc_lu == "som")
        loc = SOMMETS;
      else if (loc_lu == "elem")
        loc = ELEMENTS;
      else
        {
          Cerr << "Error for Postraitement_ft_lata::lire_champ_interface :\n";
          Cerr << loc_lu <<" has been readen. "<< finl;
          Cerr << " Keywords 'som' or 'elem' were expected after the field name '" << nom_champ << "'" << finl;
          exit();
        }

      if (!eq_interfaces.get_champ_post_FT(nom_champ, loc, (DoubleTab*) 0) && !eq_interfaces.get_champ_post_FT(nom_champ, loc, (IntTab*) 0))
        {
          Cerr << "Error for Postraitement_ft_lata::lire_champ_interface :\n";
          Cerr << " The field '" << nom_champ << "' is not understood for the " << (eq_interfaces.que_suis_je()=="Transport_Marqueur_FT"?"particules":"interfaces") << " or not authorized at localisation '";
          Nom tmp = ((loc == SOMMETS) ? "sommets" : "elements");
          Cerr << tmp << "'" << finl;
          eq_interfaces.get_champ_post_FT(demande_description, loc, (DoubleTab*) 0);
          eq_interfaces.get_champ_post_FT(demande_description, loc, (IntTab*) 0);
          exit();
        }
      Motcles& liste = loc == SOMMETS ? liste_champs_i_aux_sommets : liste_champs_i_aux_elements;
      if (!liste.contient_(nom_champ))
        liste.add(nom_champ);
    }
}

/*! @brief Build a reduced version of the facettes array, excluding virtual ones
 * Also update the internal renumbering array for later usage when writing out field values.
 */
int Postraitement_ft_lata::filter_out_virtual_fa7(IntTab& new_fa7)
{
  const Maillage_FT_Disc& mesh = refequation_interfaces.valeur().maillage_interface_pour_post();
  const IntTab& fa7 = mesh.facettes();
  int nl=fa7.dimension(0), nc=fa7.dimension(1);
  new_fa7.resize(0, nc);
  int sz = 0;

  // Reset renumbering
  renum_.clear();

  for(int i = 0; i < nl; i++)
    {
      if (!mesh.facette_virtuelle(i))
        {
          sz++;
          new_fa7.resize(sz, nc);
          renum_.push_back(i);
          for(int j = 0; j < nc; j++)
            new_fa7(sz-1, j) = fa7(i, j);
        }
    }
  return sz;
}

void Postraitement_ft_lata::filter_out_array(const DoubleTab& dtab, DoubleTab& new_dtab) const
{
  int nl = (int)renum_.size(), nc = dtab.dimension(1);
  new_dtab.resize(nl, nc);

  for(int i = 0; i < nl; i++)
    for(int j = 0; j < nc; j++)
      new_dtab(i, j) = dtab(renum_[i], j);
}

/*! @brief Write the Maillage_FT_Disc object into a LATA file in V2 format.
 *
 */
int Postraitement_ft_lata::ecrire_maillage_ft_disc()
{
  // Determining mesh type according to equation type:
  if (sub_type(Transport_Marqueur_FT,refequation_interfaces.valeur()))
    id_domaine_ = "PARTICULES";
  else if (sub_type(Transport_Interfaces_FT_Disc,refequation_interfaces.valeur()))
    id_domaine_ = "INTERFACES";
  else
    {
      Cerr<<"Type "<<refequation_interfaces.valeur().que_suis_je()<<" not recognized"<<finl;
      exit();
    }
  const Maillage_FT_Disc& mesh = refequation_interfaces.valeur().maillage_interface_pour_post();
  const DoubleTab& sommets = mesh.sommets();
  const IntTab& fa7 = mesh.facettes();
  int dim = mesh.sommets().dimension(1);
  Motcle type_elem = dim == 2 ? "Segment" : "Triangle";

  Format_Post_Lata& fpl = ref_cast(Format_Post_Lata, format_post.valeur());  // this is check above
  if (no_virtuals_)
    {
      IntTab real_fa7;
      real_fa7.set_smart_resize(1);
      filter_out_virtual_fa7(real_fa7);
      return fpl.ecrire_domaine_low_level(id_domaine_, sommets, real_fa7, type_elem);
    }
  else
    return  fpl.ecrire_domaine_low_level(id_domaine_, sommets, fa7, type_elem);
}

/*! @brief Override. Add the interfaces to the LATA output
 *
 */
int Postraitement_ft_lata::write_extra_mesh()
{
  if (refequation_interfaces.non_nul())
    {
      ecrire_maillage_ft_disc();
      return 1;
    }
  return -1;
}

void Postraitement_ft_lata::postprocess_field_values()
{
  // The standard fields
  Postraitement::postprocess_field_values();

  // Now specific FT fields:
  if(!refequation_interfaces.non_nul()) return;

  double temps_courant = mon_probleme->schema_temps().temps_courant();

  const Transport_Interfaces_FT_Disc& eq_interfaces = refequation_interfaces.valeur();
  const Maillage_FT_Disc& mesh = eq_interfaces.maillage_interface_pour_post();
  const int nb_sommets_local = mesh.sommets().dimension(0);
  const int nb_facettes_local = mesh.facettes().dimension(0);

  for (int num_loc = 0; num_loc < 2; num_loc++)
    {
      Localisation loc = SOMMETS;
      if (num_loc == 1) loc = ELEMENTS;
      const Motcles& liste = loc == SOMMETS ? liste_champs_i_aux_sommets : liste_champs_i_aux_elements;

      const int nb_champs = liste.size();
      if (nb_champs == 0) continue;

      const char * som_elem =  loc == SOMMETS ? "SOM" : "ELEM";

      DoubleTab dtab;
      IntTab itab;

      for (int i = 0; i < nb_champs; i++)
        {
          const Motcle& nom_du_champ = liste[i];
          if (eq_interfaces.get_champ_post_FT(nom_du_champ, loc, &dtab))
            {
              // ok, le champ est dans ftab
            }
          else if (eq_interfaces.get_champ_post_FT(nom_du_champ, loc, &itab))
            {
              const int n = itab.dimension(0);
              const int m = itab.nb_dim() == 1 ? 1 : itab.dimension(1);
              dtab.resize(n, m);
              // copie du champ dans dtab
              if (itab.nb_dim() == 1)
                for (int j = 0; j < n; j++)
                  dtab(j,0) = (double)itab(j);
              else
                for (int j = 0; j < n; j++)
                  for (int k = 0; k < m; k++)
                    dtab(j,k) = (double)itab(j,k);
            }
          else
            {
              // Normalement on aurait deja du detecter cette erreur a la lecture du postraitement :
              Cerr << "Error for Postraitement_ft_lata::ecrire_maillage" << finl;
              Cerr << "Unknown field : " << nom_du_champ << finl;
              exit();
            }

          const int nb_items = dtab.dimension(0);
          const int nb_items_attendus = loc == SOMMETS ? nb_sommets_local : nb_facettes_local;
          const char * msg = loc == SOMMETS ? " nb_noeuds " : " nb_facettes ";
          const int nb_compo = (dtab.nb_dim() == 1 ? 0 : dtab.dimension(1));
          if (nb_items != nb_items_attendus)
            {
              Cerr << "Error for Postraitement_ft_lata::ecrire_maillage" << finl;
              Cerr << msg << "= " << nb_items << finl;
              Cerr << msg << "expected = " << nb_items_attendus << finl;
              exit();
            }
          // [ABN] should I really bother with this?:
          Noms unites(nb_compo), noms_compo;
          std::string cp[3] = {"X", "Y", "Z"};
          for(int nc=0; nc<nb_compo; nc++)
            {
              unites.add("");
              noms_compo.add(nc < 3 ? cp[nc] : Nom(nc));
            }
          Zone dummy_dom;
          dummy_dom.nommer(id_domaine_);
          Nom nature = nb_compo == 1 ? "scalar" : "vectorial";
          const int component_to_process = -1; // meaning that we always want all components
          // For ELEMENT processing, we need to filter out virtual facettes (if requested):
          if (loc == ELEMENTS && no_virtuals_)
            {
              DoubleTab new_dtab;
              filter_out_array(dtab, new_dtab);
              postraiter_tableau(dummy_dom, unites, noms_compo, component_to_process, temps_courant, nom_du_champ, som_elem, nature,
                                 new_dtab);
            }
          else
            postraiter_tableau(dummy_dom, unites, noms_compo, component_to_process, temps_courant, nom_du_champ, som_elem, nature,
                               dtab);
        }
    }
}


