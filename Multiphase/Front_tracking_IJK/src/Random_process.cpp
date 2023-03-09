/****************************************************************************
* Copyright (c) 2021, CEA
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
// File      : Random_process.cpp
// Directory : $IJK_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Random_process.h>
#include <random>
#include <math.h>
#include <fstream>

#include <IJK_Splitting.h> // Pour que la classe de splitting soit connue
#include <communications.h>

#include <Param.h>
#include <Interprete_bloc.h>
#include <EFichier.h>
#include <SFichier.h>
#include <IJK_Lata_writer.h>
#include <IJK_Navier_Stokes_tools.h>
#include <LecFicDiffuse_JDD.h>
#include <MaillerParallel.h>
#include <EChaine.h>
#include <SChaine.h>
#include <Probleme_base.h>
#include <Domaine_VF.h>
#include <Sonde_IJK.h>
#include <Ouvrir_fichier.h>
#include <EcritureLectureSpecial.h>


Implemente_instanciable_sans_constructeur_ni_destructeur( Random_process, "Random_process", Objet_U ) ;

Sortie& Random_process::printOn( Sortie& os ) const
{

  // gr-21.07.21 - 23.07.21 :
  // Je garde une trace de toutes les idees tant que de vraies belles simus ne sont pas passees sur Occigen.
  // Idee 1 : os << gen << "\n";
  //   -> l'operateur << ne fonctionne pas entre un Sortie et un minstd_rand ni entre un sortie et un minstd_rand
  // Idee 2 : std::fstream my_os; my_os << gen; os << my_os;
  //   -> l'operateur << ne fonctionne pas entre un Sortie et un std:fstream...
  // Idee 3 : Stocker les 625 valeur de l'etat de gen dans un DoubleTab (pour un mt19937),
  //   -> je ne sais pas acceder un a un aux 625 elements de gen.
  // Idee 4 : Ecrire dans un fichier la donnée, la stocker dans un DoubleTab depuis la lecture de ce fichier,
  //          donner ce DoubleTab a la Sortie os.
  // Idee AK: Passer par un ostringstream puis convertir la string en long avec std::stol
  //
  // REMARQUES :
  //   - avec mt19937, gen contient 625 entiers codes sur 32 bits.
  //     l'ecriture dans un fichier ecrit correctement le premier chiffre. ecrit un deuxieme chiffre qui
  //     n'est pas present dans le fichier puis complete avec des zeros. MAUVAIS. (IntTab de long int ?)
  //   - avec minstd_rand alors gen contient un integer code sur 32 bit. A stocker dans un long int
  // PROBLEMES :
  //   - Passer les variables pour lecture/ecriture en attribut de classe
  //   - Travailler sur l'attribut de classe semi_gen_et_modulo_reprise_ des le printOn.

  //  *gen_write_->open(nomFichierReprise_.c_str());
  //  *gen_write_ << gen;
  //  gen_write_->close();

  std::ofstream gen_write(nom_sauvegarde_ + ".gen");
//  SFichier gen_write;
//  gen_write.ouvrir(nom_sauvegarde_);
  gen_write << gen;
  gen_write.close();

//  IntTab semi_gen_et_modulo_reprise;
  ArrOfInt semi_gen_et_modulo_reprise(2);
  long long_gen;
//  gen_read_->open(nomFichierReprise_.c_str(), std::eam::out);
//  *gen_read_ >> long_gen;
//  gen_read_->close();

//  std::ifstream gen_read(nom_du_cas()+".sauv.gen");
  EFichier gen_read;
  gen_read.ouvrir(nom_sauvegarde_ + ".gen");
  gen_read >> long_gen;
  gen_read.close();

  semi_gen_et_modulo_reprise(0) = (int)(long_gen/2); // 10 (en binaire) : "10 >> 1 = 01"
  semi_gen_et_modulo_reprise(1) = (int)(long_gen%2);
  //  IntTab like_gen;
  //  like_gen.resize(625);
  //  for (int i=0; i<625; i++)
  //    {
  //      gen_read >> like_gen(i);
  //    }
  os << "{\n";
  os << "     semi_gen_et_modulo_reprise_ " << semi_gen_et_modulo_reprise << "\n";
  os << "     process_b " << process_flt << "\n";
  os << "   }\n" ;
  return os;
}

Entree& Random_process::readOn( Entree& is )
{
  // param.ajouter n'existe pas pour un long int.
  Param param(que_suis_je());
  param.ajouter("process_b",&process_flt);

//  gen_write_->open(nomFichierReprise_.c_str());
//  *gen_write_ << like_gen;
//  gen_write_->close();
//  gen_read_->open(nomFichierReprise_.c_str());
//  *gen_read_ >> long_gen;
//  gen_read_->close();

//  std::ifstream gen_read("reprise_gen.txt");
//  gen_read >> long_gen;
//  gen_read >> gen;
//  gen_read.close();

  param.ajouter("semi_gen_et_modulo_reprise_",&semi_gen_et_modulo_reprise_,Param::REQUIRED);
  //Cout << "hello" << finl;
  //param.ajouter("distribution",&distribution);
  param.lire_avec_accolades(is);
  return is;
}


Random_process::Random_process() : nl(0),nm(0),nn(0),n_lmn(nl*nn*nm),kmin(0),kmax(0),semi_gen_et_modulo_reprise_(2)
{
  semi_gen_et_modulo_reprise_[0]=0;
  semi_gen_et_modulo_reprise_[1]=0;
}


Random_process::~Random_process()
{

}


void Random_process::initialise(double a_eps_etoile, double a_tL,
                                int a_nl, int a_nm, int a_nn, std::string nom_fichier,//, int a_random_fixed_,
                                Nom nom_sauvegarde
                               )
{
// int ind,l,m,n;
  nl = a_nl;
  nm = a_nm;
  nn = a_nn;
  n_lmn = (2*nl+1)*(2*nm+1)*(2*nn+1);

  eps_etoile = a_eps_etoile;
  tL = a_tL;
  process = set_dimensions(process,2,3,n_lmn);
// process_flt.resize(2*3*n_lmn,0.0);
  process_flt.resize(2*3*n_lmn);
  semi_gen_et_modulo_reprise_.resize(2);
  moke_gen_ = semi_gen_et_modulo_reprise_[0] ;

  // Initialisation de la serie aleatoire qui sera tiree.
  if (moke_gen_ < 0 )
    {
      // La valeur de l'etat initial est prise au hasard
      // (dans le jdd on met semi_gen_et_modulo_reprise_ et process_flt_ a 0, ou on ne les met pas du tout)
      // comme moke_gen_ est issu d'un flag, il ne doit jamais pouvoir etre <0.
      std::random_device rd {};
      std::minstd_rand gen_support(rd());
      gen = gen_support;
    }
  else if (moke_gen_ >= 0)
    {
      // L'état initial de la serie est connu mais pas le champ b. Solution simple pour determinisme
      // (dans le jeu de donnees on renseigne normalement semi_gen_et_modulo_reprise_ et on ne met pas process_flt_)
      long long_gen(2*semi_gen_et_modulo_reprise_[0] + semi_gen_et_modulo_reprise_[1]);
      std::minstd_rand gen_support(long_gen);
      gen = gen_support;
    }
  else
    {
//	  // (dans le jdd, semi_gen_et_modulo_reprise_ et process_flt_ sont connus et renseignes)
//      // l'etat initial est connu a partir de semi_gen_et_modulo_reprise_ et le champ b est connu a partir de process_flt
//      long long_gen(2*semi_gen_et_modulo_reprise_[0] + semi_gen_et_modulo_reprise_[1]);
//      std::minstd_rand gen_support(long_gen);
//      gen = gen_support;
//
      Cerr << "moke_gen a une valeur anormale" << finl;
    }

  std::normal_distribution < > dist_support(0.,1.);
  distribution = dist_support;

//  std::ofstream detail_gen(nom_fichier);
  nom_sauvegarde_ = nom_sauvegarde;
  nom_fichier_ = nom_fichier;
  std::ofstream Detail_gen(nom_fichier_.c_str());
  if(Detail_gen)
    {
      Detail_gen << "-------- RANDOM_PROCESS; detail de gen --------" << std::endl;
//	  Detail_gen << std::endl << "l,m,n \t : b_x, \tb_y, \tb_z,\t";
    }

//  nomFichierReprise_ = "reprise_gen.sauv";
//  *gen_write_ = nomFichierReprise_.c_str();
//  *gen_read_ = nomFichierReprise_.c_str();

}


void Random_process::next_step(double dt, int it)
{
  // Ajouter le lien d'ou j'ai pompe
  // --> Plus bien voir comment FIXER l'etat aleatoire que l'on va parcourir
  // --> Plus donner la possibilite dans le jdd de jouer sur la serie aleatoire tiree
  // std::random_device rd;
  // std::minstd_rand gen(rd());
  // std::minstd_rand gen(it);
  // std::normal_distribution<double> distribution(0.0,1.0);
  double Gaussian[2][3];
  std::vector< std::vector< std:: vector < double > > > old_process;
  int cpx, dir, l, m, n, ind;

  old_process = process;
  for (n=0; n<2*nn+1; n++)
    for (m=0; m<2*nm+1; m++)
      // TODO : Distinguer cas pair de cas impair !!! ATTN
      for (l=nl; l<2*nl+1; l++) // La force sp est a symmetrie Hermitienne, donc le rp aussi. On ne parcourt donc que la moitiee du tout
        {
          ind = (n*(2*nm+1) + m) * (2*nl+1) +l;
          // std::cout << "ind : " << ind << "/" << n_lmn << std::endl;
          // std::cout << "otr : " << n_lmn-ind << "/" << n_lmn << std::endl;
          for (dir=0; dir<3; dir++)
            {
              for (cpx=0; cpx<2; cpx++)
                {
                  Gaussian[cpx][dir] = distribution(gen);
                  process[cpx][dir][ind] = old_process[cpx][dir][ind]*(1-dt/tL);
                  process[cpx][dir][ind] += Gaussian[cpx][dir]*sqrt(2*eps_etoile*dt*pow(tL,2));
                }
              // Symetrie Hermitienne. Sous nos conventions, -k[ind] = k[n_lmn - ind]
              // Remarque : on ne se sert meme pas de l'autre moitie.
              //            on la calcule quand meme pour faciliter les ecritures
              process[0][dir][n_lmn-ind] =   process[0][dir][ind];
              process[1][dir][n_lmn-ind] = - process[1][dir][ind];
            }
        }
}



void Random_process::next_step2(double dt, int it)
{
  if (Process::je_suis_maitre())
    {
      // Si sorties supplementaires est mis a 1, on ecrit les etats valeurs de seed dans un fichier.
      //  --> Utile que pour debuguer/ attester de façon extremement précise que les reprises sont "justes"
      int sorties_supplementaires(0);
      std::ofstream Detail_gen(nom_fichier_.c_str(),std::ios_base::app);
      int n_dir(3);
      // source : https://en.cppreference.com/w/cpp/numeric/random/normal_distribution
      // Attention, initialiser de nouvelles graines a chaque fois c'est pas la bonne utilisation
      double Gaussian[2][3];

      // std:: vector < double > old_process;
      ArrOfDouble old_process(process_flt);
      // old_process = process_flt;
      // process_flt.resize(2*3*n_lmn,0.0);
      process_flt = 0.;

      for (int n=0; n<2*nn+1; n++)
        for (int m=0; m<2*nm+1; m++)
          // TODO : Distinguer cas pair de cas impair !!! ATTN
          for (int l=0; l<2*nl+1; l++) // La force sp est a symmetrie Hermitienne, donc le rp aussi. On ne parcourt donc que la moitiee du tou
            {
              int ind_lmn = (n*(2*nm+1) + m) * (2*nl+1) +l;
              for (int dir=0; dir<3; dir++)
                {
                  // double a(distribution(gen));
                  // Gaussian[0][dir] = cos(a);
                  // Gaussian[1][dir] = sin(a);
                  for (int cpx=0; cpx<2; cpx++)
                    {
                      int ind_CDIlmn((cpx*n_dir+dir)*n_lmn+ind_lmn);
                      if (sorties_supplementaires)
                        Detail_gen << gen << " ";
                      Gaussian[cpx][dir] = distribution(gen);  // Pour avoir une variance à 1, pas à 2.
                      if (sorties_supplementaires)
                        Detail_gen << gen << std::endl;
                      process_flt[ind_CDIlmn] = old_process[ind_CDIlmn]*(1-dt/tL);
                      process_flt[ind_CDIlmn] += Gaussian[cpx][dir]*sqrt(2.0*eps_etoile*(0.0+dt)/pow(tL,2));
                    }
                  // Symetrie Hermitienne. Sous nos conventions, -k[ind] = k[n_lmn - ind]
                  // Remarque : on ne se sert meme pas de l'autre moitie.
                  //            on la calcule quand meme pour faciliter les ecritures
                  // process_flt[(0*n_dir+dir)*n_lmn+(ind_lmn_moins)] =   process_flt[(0*n_dir+dir)*n_lmn+(ind_lmn)];
                  // process_flt[(1*n_dir+dir)*n_lmn+(ind_lmn_moins)] = - process_flt[(1*n_dir+dir)*n_lmn+(ind_lmn)];
                }
            }
      if (sorties_supplementaires)
        Detail_gen.close();
    }
  envoyer_broadcast(process_flt,0);
}

void Random_process::next_step3(ArrOfDouble& advection_velocity, double dt, int it)
{
  // On ajoute l'advection du forcage uniquement sur processus aleatoire
  // On ajoute aussi la condition quon mene une action uniquement si |k| > k_min
  if (Process::je_suis_maitre())
    {
      std::ofstream Detail_gen(nom_fichier_.c_str(),std::ios_base::app);
      int n_dir(3);
      // source : https://en.cppreference.com/w/cpp/numeric/random/normal_distribution
      // Attention, initialiser de nouvelles graines a chaque fois c'est pas la bonne utilisation
      double Gaussian[2][3];

      // std:: vector < double > old_process;
      ArrOfDouble old_process(process_flt);
      // old_process = process_flt;
      // process_flt.resize(2*3*n_lmn,0.0);
      process_flt = 0.;
      double kappa[3];
      double coefficient_translation[2*3];

      for (int n=0; n<2*nn+1; n++)
        for (int m=0; m<2*nm+1; m++)
          // TODO : Distinguer cas pair de cas impair !!! ATTN
          for (int l=0; l<2*nl+1; l++) // La force sp est a symmetrie Hermitienne, donc le rp aussi. On ne parcourt donc que la moitiee du tout
            {
              int ind_lmn = (n*(2*nm+1) + m) * (2*nl+1) +l;
              int ind_lmn_moins = (n*(2*nm+1) + m) * (2*nl+1) + (2*nl-l);

              // Position dans l'espace spectral. k_x va de -kmin a +kmin,
              // d'ou le l-nl
              kappa[0] = - kmax + (l)*(2*kmax)/(2*nl);
              kappa[1] = - kmax + (m)*(2*kmax)/(2*nm);
              kappa[2] = - kmax + (n)*(2*kmax)/(2*nn);


              if (!(std::fabs(kappa[0])<kmin && std::fabs(kappa[1])<kmin && std::fabs(kappa[2])<kmin))
                {
                  for (int dir=0; dir<3; dir++)
                    {
                      // double a(distribution(gen));
                      // Gaussian[0][dir] = cos(a);
                      // Gaussian[1][dir] = sin(a);
                      int ind_RD(0*n_dir+dir);
                      int ind_CD(1*n_dir+dir);
                      int ind_CDIlmn((1*n_dir+dir)*n_lmn+ind_lmn);
                      int ind_RDIlmn((0*n_dir+dir)*n_lmn+ind_lmn);

                      // ATTENTION ! On fait ici une translation de u*dt or il faut faire une translation de u*t
                      coefficient_translation[ind_RD] = cos(-1.*kappa[dir]*dt*advection_velocity[dir]);
                      coefficient_translation[ind_CD] = sin(-1.*kappa[dir]*dt*advection_velocity[dir]);

                      // Translation de l'etat precedent
                      process_flt[ind_RDIlmn] = old_process[ind_RDIlmn]*coefficient_translation[ind_RD] - old_process[ind_CDIlmn]*coefficient_translation[ind_CD];
                      process_flt[ind_CDIlmn] = old_process[ind_RDIlmn]*coefficient_translation[ind_CD] + old_process[ind_CDIlmn]*coefficient_translation[ind_RD];

                      for (int cpx=0; cpx<2; cpx++)
                        {
                          int ind_RCDIlmn((cpx*n_dir+dir)*n_lmn+ind_lmn);
                          Detail_gen << gen << " ";
                          Gaussian[cpx][dir] = distribution(gen);
                          Detail_gen << gen << std::endl;
                          process_flt[ind_RCDIlmn] = process_flt[ind_RCDIlmn]*(1-dt/tL);
                          process_flt[ind_RCDIlmn] += Gaussian[cpx][dir]*sqrt(2*eps_etoile*(0.0+dt)/pow(tL,2));
                        }
                      // Symetrie Hermitienne. Sous nos conventions, -k[ind] = k[n_lmn - ind]
                      // Remarque : on ne se sert meme pas de l'autre moitie.
                      //            on la calcule quand meme pour faciliter les ecritures
                      process_flt[(0*n_dir+dir)*n_lmn+(ind_lmn_moins)] =   process_flt[(0*n_dir+dir)*n_lmn+(ind_lmn)];
                      process_flt[(1*n_dir+dir)*n_lmn+(ind_lmn_moins)] = - process_flt[(1*n_dir+dir)*n_lmn+(ind_lmn)];
                    }
                }
            }
      Detail_gen.close();
    }
  envoyer_broadcast(process_flt,0);
}


void Random_process::write(std::string nom_fichier_sortie, double t)
{
  std::ofstream Random_flux(nom_fichier_sortie.c_str(), std::ios::app);
  if (Random_flux)
    {
      // int cpx, dir;
      int l, m, n, ind;

      Random_flux << std::endl << "time : " << t << std::endl << std::endl;
      for (n=0; n<2*nn+1; n++)
        for (m=0; m<2*nm+1; m++)
          for (l=0; l<2*nl+1; l++)
            {
              ind = (n*(2*nm+1) + m) * (2*nl+1) +l;
              Random_flux << l<<","<<m<<","<<n<<"\t : ";
              Random_flux << process[0][0][ind] << " + i" << process[1][0][ind]<<", \t";
              Random_flux << process[0][1][ind] << " + i" << process[1][1][ind]<<", \t";
              Random_flux << process[0][2][ind] << " + i" << process[1][2][ind]<<", \t";
              Random_flux << std::endl;
            }
    }
}

void Random_process::write_separate(std::string nom_fichier_sortie, double t)
{
  std::ofstream Random_flux(nom_fichier_sortie.c_str());
  if (Random_flux)
    {
      // int cpx, dir;
      int l, m, n, ind;
      Random_flux << std::endl << "l,m,n, \t  rb_x, \t cb_x, \t\t rb_y, \t cb_y, \t\t rb_z, \t cb_z \t";
      Random_flux << std::endl;

      // Random_flux << std::endl << "time : " << t << std::endl << std::endl;
      for (n=0; n<2*nn+1; n++)
        for (m=0; m<2*nm+1; m++)
          for (l=0; l<2*nl+1; l++)
            {
              ind = (n*(2*nm+1) + m) * (2*nl+1) +l;
              Random_flux << l<<","<<m<<","<<n<<"\t,";
              Random_flux << process[0][0][ind] << ",\t" << process[1][0][ind]<<", \t\t";
              Random_flux << process[0][1][ind] << ",\t" << process[1][1][ind]<<", \t\t";
              Random_flux << process[0][2][ind] << ",\t" << process[1][2][ind]<<"\t";
              Random_flux << std::endl;
            }
    }
}

std::vector< std::vector< std:: vector < double >>> Random_process::get_b()
{
  return process;
}

// std:: vector < double > Random_process::get_b_flt()
ArrOfDouble& Random_process::get_b_flt()
{
  return process_flt;
}

int Random_process::get_semi_gen()
{
  return semi_gen_et_modulo_reprise_[0];
}
