/****************************************************************************
* Copyright (c) 2022, CEA
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
// File      : DebogFT.cpp
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////
#include <DebogFT.h>
// #include <IJK_Field.h>
#include <Param.h>
#include <TRUSTTabs.h>
#include <communications.h>
#include <Maillage_FT_Disc.h>

Implemente_instanciable(DebogFT,"DebogFT",Interprete) ;
int DebogFT::debog_mode_ = DebogFT::DISABLED;
Nom DebogFT::filename_;
EFichier DebogFT::infile_;
SFichier DebogFT::outfile_;
double DebogFT::seuil_absolu_ = 1e-4;
double DebogFT::seuil_relatif_ = 1e-8;
double DebogFT::seuil_minimum_relatif_ = 1e-4;
Sortie& DebogFT::printOn(Sortie& os) const
{
  Objet_U::printOn(os);
  return os;
}
Entree& DebogFT::readOn(Entree& is)
{
  return is;
}
Entree& DebogFT::interpreter(Entree& is)
{
  filename_ = "DEBOG.txt";
  Param param(que_suis_je());
  param.ajouter("mode", &debog_mode_);
  param.dictionnaire("disabled", (int)DISABLED);
  param.dictionnaire("write_pass", (int)WRITE_PASS);
  param.dictionnaire("check_pass", (int)CHECK_PASS);
  param.ajouter("filename", &filename_);
  param.ajouter("seuil_absolu", &seuil_absolu_);
  param.ajouter("seuil_relatif", &seuil_relatif_);
  param.ajouter("seuil_minimum_relatif", &seuil_minimum_relatif_);
  param.lire_avec_accolades(is);
  if (je_suis_maitre())
    {
      if (debog_mode_ == WRITE_PASS)
        outfile_.ouvrir(filename_);
      if (debog_mode_ == CHECK_PASS)
        infile_.ouvrir(filename_);
    }
  return is;
}

// Calcule une signature numerique du champ. On veut une signature qui varie
// lineairement avec les valeurs du champ et qui detecte de petites variations.
// On va calculer quelques sommes ponderees des valeurs du champ avec des
// ponderations pseudo-aleatoires.
void DebogFT::compute_signature_sommets(const Maillage_FT_Disc& mesh, const ArrOfDouble& data, ArrOfDouble& signature)
{
  const int sig_size = 5;
  signature.resize_array(sig_size);
  ArrOfDouble facteurs(sig_size);

  // Generation de facteurs, nombres irrationnels et non multiples entre eux
  facteurs[0] = 1.35914091422952;
  for (int i = 1; i < sig_size; i++)
    facteurs[i] = facteurs[i-1] * facteurs[0];
  for (int i = 0; i < sig_size; i++)
    signature[i] = 0.;
  const int nbsommets = mesh.nb_sommets();
  const DoubleTab& sommets = mesh.sommets();
  for (int i = 0; i < nbsommets; i++)
    {
      if (mesh.sommet_virtuel(i))
        continue;
      const double d = data[i];
      const double coord_x = sommets(i,0);
      const double coord_y = sommets(i,1);
      const double coord_z = sommets(i,2);
      for (int l = 0; l < sig_size; l++)
        {
          const double f = facteurs[l];
          signature[l] += d * cos(f * coord_x) * cos(f * coord_y) * cos(f * coord_z);
        }
    }
  mp_sum_for_each_item(signature);
}

void DebogFT::verifier_maillage_ft(const char *msg, const Maillage_FT_Disc& mesh)
{
  if (debog_mode_ == DISABLED)
    return;
  ArrOfDouble sig;
  ArrOfDouble data(mesh.nb_sommets());
  data = 1.;
  compute_signature_sommets(mesh, data, sig);

  DebogFT::verifier_(msg, sig);
}
void DebogFT::verifier_tableau_sommets(const char *msg, const Maillage_FT_Disc& mesh, const ArrOfDouble& data)
{
  if (debog_mode_ == DISABLED)
    return;
  ArrOfDouble sig;
  compute_signature_sommets(mesh, data, sig);

  DebogFT::verifier_(msg, sig);
}
void DebogFT::verifier_tableau_sommets(const char *msg, const Maillage_FT_Disc& mesh, const DoubleTab& tab)
{
  if (debog_mode_ == DISABLED)
    return;
  ArrOfDouble sig0, sig1;
  ArrOfDouble data(mesh.nb_sommets());
  for (int compo = 0; compo < tab.dimension(1); compo++)
    {
      const int n = mesh.nb_sommets();
      for (int i = 0; i < n; i++)
        data[i] = tab(i,compo);
      compute_signature_sommets(mesh, data, sig1);
      if (compo == 0)
        sig0 = sig1;
      else
        sig0 += sig1;
    }

  DebogFT::verifier_(msg, sig0);
}

int debogft_ijk_compteur_stop_if = -1; // mettre ici le numero du passage ou on veut s'arreter avec gdb

void DebogFT::verifier_(const char *msg, const ArrOfDouble& sig)
{
  static int compteur = 0;
  compteur++;
  // (faire "set compteur_stop_if=valeur" dans gdb)
  // et mettre un point d'arret sur le Cerr:
  if (debogft_ijk_compteur_stop_if == compteur)
    {
      Cerr << "DebogFT::verifier_, compteur = " << compteur << "On y est !" << finl;
    }

  if (Process::je_suis_maitre())
    {
      Nom s("");
      char ss[1000];
      for (int i = 0; i < sig.size_array(); i++)
        {
          sprintf(ss, "%24.16e ", sig[i]);
          s += ss;
        }
      s += msg;
      Journal() << "DEBOG1:" << s << finl;
      if (debog_mode_ == WRITE_PASS)
        outfile_ << s << finl;
      else
        {
          ArrOfDouble sig2(sig.size_array());
          Nom s2("");
          for (int i = 0; i < sig.size_array(); i++)
            {
              infile_ >> sig2[i];
              sprintf(ss, "%24.16e ", sig2[i]);
              s2 += ss;
            }
          std::string ligne;
          std::getline(infile_.get_ifstream(), ligne);
          Journal() << "DEBOG2:" << s2 << ligne.c_str() << finl;
          bool erreur = false;
          for (int i = 0; i < sig.size_array(); i++)
            {
              double m = std::max(fabs(sig[i]),fabs(sig2[i]));
              m = std::max(seuil_minimum_relatif_, m);
              if (fabs(sig[i] - sig2[i]) > seuil_absolu_
                  || fabs(sig[i] - sig2[i]) / m > seuil_relatif_)
                erreur = true;
            }
          if (erreur)
            {
              Cerr << "DEBOG: erreur" << finl << "THIS:" << s << finl << "REF: " << s2 << finl;
            }
        }
    }
}
