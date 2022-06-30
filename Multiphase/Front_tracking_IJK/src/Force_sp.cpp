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
// File      : Force_sp.cpp
// Directory : $IJK_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Force_sp.h>
#include <fstream>
#include <math.h>
#include <IJK_Splitting.h>
#include <communications.h>
#include <MaillerParallel.h>
#include <EChaine.h>
#include <SChaine.h>
#include <Probleme_base.h>
#include <Zone_VF.h>
#include <IJK_Navier_Stokes_tools.h>

Implemente_instanciable_sans_constructeur_ni_destructeur( Force_sp, "Force_sp", Objet_U ) ;

Sortie& Force_sp::printOn( Sortie& os ) const
{
  Objet_U::printOn( os );
  return os;
}

Entree& Force_sp::readOn( Entree& is )
{
  Objet_U::readOn( is );
  return is;
}


Force_sp::Force_sp() : nl(0),nm(0),nn(0),n_lmn(nl*nn*nm),kmin(0),kmax(0), energie(0)
{

}

Force_sp::~Force_sp()
{

}

void Force_sp::initialise(int a_nl,int a_nn,int a_nm, int a_momin, int a_momax, double a_kmin, double a_kmax, std::string nom_fichier
                         )
{
  nl = a_nl;
  nm = a_nm;
  nn = a_nn;
  n_lmn = (2*nl+1)*(2*nm+1)*(2*nn+1);
  momin = a_momin;
  momax = a_momax;
  kmin = a_kmin;
  kmax = a_kmax;
  amplitude = 0.;


// force[0] -> Partie reelle;
// force[1] -> partie imaginaire;
  force.resize(2);
// force[_][i] -> Composant i = x, y ou z
  force[0].resize(3);
  force[1].resize(3);
  force[0][0].resize(n_lmn,0.0);
  force[1][0].resize(n_lmn,0.0);
  force[0][1].resize(n_lmn,0.0);
  force[1][1].resize(n_lmn,0.0);
  force[0][2].resize(n_lmn,0.0);
  force[1][2].resize(n_lmn,0.0);

  force_flt.resize_array(2*3*n_lmn);

  std::ofstream Spectral_flux(nom_fichier.c_str());
  if(Spectral_flux)
    {
      Spectral_flux << "-------- SECTRAL_FORCE --------" << std::endl;
      Spectral_flux << std::endl << "l,m,n \t : f_x, \tf_y, \tf_z,\t";
    }

  energie = 0.;

}


void Force_sp::initialise(int a_nl,int a_nn,int a_nm, int a_momin, int a_momax, double a_kmin, double a_kmax, double a_amplitude, std::string nom_fichier)
{

  nl = a_nl;
  nm = a_nm;
  nn = a_nn;
  n_lmn = (2*nl+1)*(2*nm+1)*(2*nn+1);
  momin = a_momin;
  momax = a_momax;
  kmin = a_kmin;
  kmax = a_kmax;
  amplitude = a_amplitude;

// force[0] -> Partie reelle;
// force[1] -> partie imaginaire;
  force.resize(2);
// force[_][i] -> Composant i = x, y ou z
  force[0].resize(3);
  force[1].resize(3);

  force[0][0].resize(n_lmn,0.0);
  force[1][0].resize(n_lmn,0.0);
  force[0][1].resize(n_lmn,0.0);
  force[1][1].resize(n_lmn,0.0);
  force[0][2].resize(n_lmn,0.0);
  force[1][2].resize(n_lmn,0.0);

  force_flt.resize_array(2*3*n_lmn);

  std::ofstream Spectral_flux(nom_fichier.c_str());
  if(Spectral_flux)
    {
      Spectral_flux << "-------- SECTRAL_FORCE --------" << std::endl;
      Spectral_flux << std::endl << "l,m,n \t : f_x, \tf_y, \tf_z,\t";
    }

}



void Force_sp::compute_step2(ArrOfDouble& process_flt)
{
  if (Process::je_suis_maitre())
    {
      /*
      Version avancée du calcul du processus aléatoire :
        - Les champs sont stockés dans des tableaux à une dimension,
          l'indice est reconstruit.
      */
      int n_dir(3);
      double kappa[3];
      double terme[2];
      for (int n=0; n<2*nn+1; n++)
        for (int m=0; m<2*nm+1; m++)
          for (int l=nl; l<2*nl+1; l++) // La TF-1 de cette force est réelle. Donc cette force est à symétrie Hermitienne
            {
              int ind_lmn = (n*(2*nm+1) + m) * (2*nl+1) +l;
              int ind_lmn_moins = (n*(2*nm+1) + m) * (2*nl+1) + (2*nl-l);
              // Position dans l'espace spectral. k_x va de -kmin a +kmin,
              // d'ou le l-nl
              kappa[0] = - kmax + (l)*(2*kmax)/(2*nl);
              kappa[1] = - kmax + (m)*(2*kmax)/(2*nm);
              kappa[2] = - kmax + (n)*(2*kmax)/(2*nn);

              if (abs(kappa[0])<kmin && abs(kappa[1])<kmin && abs(kappa[2])<kmin)
                {
                  /* Le nombre d'onde est trop petit pour etre dans le domaine force*/
                  for (int cpx=0; cpx<2; cpx++)
                    {
                      for (int dir=0; dir<3; dir++)
                        {
                          int ind_CDI((cpx*n_dir+dir)*n_lmn+ind_lmn);
                          force_flt[ind_CDI] = 0;
                          force_flt[(0*n_dir+dir)*n_lmn+(ind_lmn_moins)] =   force_flt[(0*n_dir+dir)*n_lmn+(ind_lmn)];
                          force_flt[(1*n_dir+dir)*n_lmn+(ind_lmn_moins)] = - force_flt[(1*n_dir+dir)*n_lmn+(ind_lmn)];
                        }
                    }
                }
              else
                {
                  /* On est dans le domaine force */
                  double norme_kappa = sqrt(kappa[0]*kappa[0] + kappa[1]*kappa[1] + kappa[2]*kappa[2]);
                  for (int cpx=0; cpx<2; cpx++)
                    {
                      terme[cpx] =  kappa[0]*process_flt[(cpx*n_dir+0)*n_lmn+ind_lmn];
                      terme[cpx] += kappa[1]*process_flt[(cpx*n_dir+1)*n_lmn+ind_lmn];
                      terme[cpx] += kappa[2]*process_flt[(cpx*n_dir+2)*n_lmn+ind_lmn];
                      terme[cpx] /= pow(norme_kappa,2);
                      for (int dir=0; dir<3; dir++)
                        {
                          int ind_CDI((cpx*n_dir+dir)*n_lmn+ind_lmn);
                          force_flt[ind_CDI] = process_flt[ind_CDI] - kappa[dir]*terme[cpx];
                          // Symétrie Hermitienne : Sous nos conventions, -k[ind] = k[n_lmn - ind]
                          force_flt[(0*n_dir+dir)*n_lmn+(ind_lmn_moins)] =   force_flt[(0*n_dir+dir)*n_lmn+(ind_lmn)];
                          force_flt[(1*n_dir+dir)*n_lmn+(ind_lmn_moins)] = - force_flt[(1*n_dir+dir)*n_lmn+(ind_lmn)];
                        }
                    }
                }
            }
    }
  envoyer_broadcast(force_flt,0);
}

void Force_sp::field_advection(const ArrOfDouble& advection_length, const double dt, const int it)
{
  if (it>0)
    {
      if (Process::je_suis_maitre())
        {
          /*
          Version avancée du calcul du processus aléatoire :
            - Les champs sont stockés dans des tableaux à une dimension,
              l'indice est reconstruit.
          */
          int n_dir(3);
          double kappa[3];
          double coefficient_translation[2];
          for (int n=0; n<2*nn+1; n++)
            for (int m=0; m<2*nm+1; m++)
              for (int l=0; l<2*nl+1; l++) //for (int l=nl; l<2*nl+1; l++) // La TF-1 de cette force est réelle. Donc cette force est à symétrie Hermitienne
                {
                  int ind_lmn = (n*(2*nm+1) + m) * (2*nl+1) +l;
                  // Position dans l'espace spectral. k_x va de -kmin a +kmin,
                  // d'ou le l-nl
                  kappa[0] = - kmax + (l)*(2*kmax)/(2*nl);
                  kappa[1] = - kmax + (m)*(2*kmax)/(2*nm);
                  kappa[2] = - kmax + (n)*(2*kmax)/(2*nn);

                  if (!(abs(kappa[0])<kmin && abs(kappa[1])<kmin && abs(kappa[2])<kmin))
                    {
                      /* Le nombre d'onde est dans le domaine force*/


                      double argument_advection(0);
                      for (int dor=0; dor<n_dir; dor++)
                        {
                          /* Produit scalaire : kappa . advection_length */
                          Cout << dor << "argument_advection" << argument_advection << finl;
                          argument_advection += kappa[dor]*advection_length[dor];
                        }
                      Cout << "argument_advection" << argument_advection << finl;
                      coefficient_translation[0] = cos(-1.*argument_advection);
                      coefficient_translation[1] = sin(-1.*argument_advection);

                      for (int dir=0; dir<n_dir; dir++)
                        {
                          // Translation dans le domaine physique <=> multiplication par expo cpx en spectral
                          // TF(f(x-u Dt)) = exp(-i kappa u Dt) * TF(f(x))

                          {
                            int ind_CDIlmn((1*n_dir+dir)*n_lmn+ind_lmn);
                            int ind_RDIlmn((0*n_dir+dir)*n_lmn+ind_lmn);

                            force_flt[ind_RDIlmn] = force_flt[ind_RDIlmn]*coefficient_translation[0] - force_flt[ind_CDIlmn]*coefficient_translation[1];
                            force_flt[ind_CDIlmn] = force_flt[ind_CDIlmn]*coefficient_translation[0] + force_flt[ind_RDIlmn]*coefficient_translation[1];

                            // Symétrie Hermitienne : Sous nos conventions, -k[ind] = k[n_lmn - ind]
//                            force_flt[(0*n_dir+dir)*n_lmn+(ind_lmn_moins)] =   force_flt[(0*n_dir+dir)*n_lmn+(ind_lmn)];
//                            force_flt[(1*n_dir+dir)*n_lmn+(ind_lmn_moins)] = - force_flt[(1*n_dir+dir)*n_lmn+(ind_lmn)];
                          }

                        }
                    }
                }
        }
      envoyer_broadcast(force_flt,0);
    }
}


void Force_sp::set_zero()
{
  /*Mise à zero de tout le champ de force spectral*/
  n_lmn = (2*nl+1)*(2*nm+1)*(2*nn+1);

  force_flt.resize(2*3*n_lmn);//,0.0);
}


void Force_sp::compute_force_kappa()
{
  // Seul le proc 0 fait le calcul, comme c'est le mm clc pour tous les procs
  if (Process::je_suis_maitre())
    {
      /*
      Fonction de test, non opérationnelle pour le moment. Peut facilement
      ettre opérationnelle.
      Construit une force spectrale : f_sp(k) = k
      */
      int cpx, dir, l, m, n, ind;
      double kappa[3];
      double norme_kappa;
      for (n=0; n<2*nn+1; n++)
        {
          for (m=0; m<2*nm+1; m++)
            for (l=0; l<2*nl+1; l++) // La TF-1 de cette force est réelle. Donc cette force est à symétrie Hermitienne
              {
                ind = (n*(2*nm+1) + m) * (2*nl+1) +l;

                kappa[0] = - kmax + (l)*(2*kmax)/(2*nl);
                kappa[1] = - kmax + (m)*(2*kmax)/(2*nm);
                kappa[2] = - kmax + (n)*(2*kmax)/(2*nn);
                norme_kappa = sqrt(kappa[0]*kappa[0] + kappa[1]*kappa[1] + kappa[2]*kappa[2]);
                if (norme_kappa >= kmin)
                  {
                    for (dir=0; dir<3; dir++)
                      for (cpx=0; cpx<2; cpx++)
                        {
                          force[cpx][dir][ind] = kappa[0];
                        }
                  }
              }
        }
    }
  envoyer_broadcast(force_flt,0);
}

void Force_sp::compute_dirac_board()
{
  // Seul le proc 0 fait le calcul, comme c'est le mm clc pour tous les procs
  if (Process::je_suis_maitre())
    {
      /*
      Fonction de test. Non opérationnelle à ce jour
      Construit une force spectral qui vaut 1 sur un plan et 0 ailleurs
      */
      int ind;
      double kappa[3];
      double norme_kappa;
      for (int n=0; n<2*nn+1; n++)
        {

          for (int m=0; m<2*nm+1; m++)
            for (int l=0; l<2*nl+1; l++) // La TF-1 de cette force est réelle. Donc cette force est à symétrie Hermitienne
              {
                ind = (n*(2*nm+1) + m) * (2*nl+1) +l;
                kappa[0] = - kmax + (l)*(2*kmax)/(2*nl);
                kappa[1] = - kmax + (m)*(2*kmax)/(2*nm);
                kappa[2] = - kmax + (n)*(2*kmax)/(2*nn);
                norme_kappa = sqrt(kappa[0]*kappa[0] + kappa[1]*kappa[1] + kappa[2]*kappa[2]);

                if (norme_kappa >= kmin)
                  {
                    for (int dir=0; dir<3; dir++)
                      {
                        // ICI ON FORCE SUR LES 2 PLAQUES QUI SONT  X=Xo. Or ce n'est pas ce qu'on veut !
                        if (l==nl+nl/2 || l==nl/2)
                          {
                            force[0][dir][ind] = amplitude;
                            force[1][dir][ind] = 0;
                          }
                      }
                  }
              }
        }
    }
  envoyer_broadcast(force_flt,0);
}




void Force_sp::compute_dirac_point()
{
  // Seul le proc 0 fait le calcul, comme c'est le mm clc pour tous les procs
  if (Process::je_suis_maitre())
    {
      /*
      Fonciton de test. Operationnelle
      Construit la force spectrale qui vaut 1 en deux point et 0 ailleurs
      de sorte à ce que la fonction de forçage physique soit : cos(x+y) ex
      */
      int dir(0);

      int roc(momin+1);

      int l1(nl-roc), l2(nl+roc);
      int n1(nn-roc), n2(nn+roc);
      int fsp_m1(nm), fsp_m2(nm);

      int ind1((n1*(2*nm+1) + fsp_m1) * (2*nl+1) +l1), ind2((n2*(2*nm+1) + fsp_m2) * (2*nl+1) +l2);
      force[0][dir][ind1]=amplitude;
      force[0][dir][ind2]=amplitude;
    }
  envoyer_broadcast(force_flt,0);
}

void Force_sp::compute_dirac_point_div_nulle()
{
  // Seul le proc 0 fait le calcul, comme c'est le mm clc pour tous les procs
  if (Process::je_suis_maitre())
    {
      /*
      Fonction de test. Operationnelle
      Construit la force spectrale qui vaut 1 en deux point et 0 ailleurs
      de sorte à ce que la fonction de forçage physique soit : cos(x-y) (ex + ey)
      */

      int n_dir(3);
      int roc(momin+0);
      int n_ind_lmn(n_lmn);

      int l_moins(nl-roc), l_plus(nl+roc);
      int n_moins(nn-roc), n_plus(nn+roc);
      int m_moins(nm);
      int m_plus(nm);

      int ind_moins((n_moins*(2*nm+1) + m_moins) * (2*nl+1) +l_plus);
      int ind_plus((n_plus*(2*nm+1) + m_plus) * (2*nl+1) +l_moins);

      for (int dir=0; dir<2; dir++)
        {
          dir = 2*dir;
          double coeff(0);
          if (dir==0)
            coeff=1;
          else
            coeff=-1;

          int ind_RDImoins((0*n_dir+dir)*n_ind_lmn+ind_moins);
          int ind_RDIplus((0*n_dir+dir)*n_ind_lmn+ind_plus);
          force_flt[ind_RDImoins]=amplitude*(1./sqrt(2))*coeff;
          force_flt[ind_RDIplus]=amplitude*(1./sqrt(2))*coeff;

          ind_RDImoins = (0*n_dir+dir+2)*n_ind_lmn+ind_moins;
          ind_RDIplus = (0*n_dir+dir+2)*n_ind_lmn+ind_plus;
          force_flt[ind_RDImoins]=amplitude*(1./sqrt(2))*coeff;
          force_flt[ind_RDIplus]=amplitude*(1./sqrt(2))*coeff;
        }
    }
  envoyer_broadcast(force_flt,0);
}

void Force_sp::compute_dirac_point_uniZ()
{

  // Seul le proc 0 fait le calcul, comme c'est le mm clc pour tous les procs
  if (Process::je_suis_maitre())
    {
      // donne cos(z) (ex)
      int dir(2);
      int n_dir(3);

      int roc(momin+0);
      int n_ind_lmn(n_lmn);

      // int l_moins(nl-roc);
      int l_zero(nl);
      // int l_plus(nl+roc);
      // int m_moins(nm-roc);
      int m_zero(nm);
      // ints m_plus(nm+roc);
      int n_moins(nn-roc);
      // int n_zero(nn);
      int n_plus(nn+roc);

      int ind_moins((n_moins*(2*nm+1) + m_zero) * (2*nl+1) +l_zero);
      int ind_plus((n_plus*(2*nm+1) + m_zero) * (2*nl+1) +l_zero);


      int ind_RDImoins((0*n_dir+dir)*n_ind_lmn+ind_moins);
      int ind_RDIplus((0*n_dir+dir)*n_ind_lmn+ind_plus);
      force_flt[ind_RDImoins]=amplitude;
      force_flt[ind_RDIplus]=amplitude;

      force[0][dir][ind_moins]=amplitude;
      force[0][dir][ind_plus]=amplitude;
    }
  envoyer_broadcast(force_flt,0);
}


void Force_sp::compute_dirac_point_uniY()
{
  // Seul le proc 0 fait le calcul, comme c'est le mm clc pour tous les procs
  if (Process::je_suis_maitre())
    {
      // donne cos(y) (ex)
      int dir(1);
      int n_dir(3);
      int roc(momin+0);
      int n_ind_lmn(n_lmn);

      // int l_moins(nl-roc);
      int l_zero(nl);
      // int l_plus(nl+roc);
      int m_moins(nm-roc);
      // int m_zero(nm);
      int m_plus(nm+roc);
      // int n_moins(nn-roc);
      int n_zero(nn);
      // int n_plus(nn+roc);

      int ind_moins((n_zero*(2*nm+1) + m_moins) * (2*nl+1) +l_zero);
      int ind_plus((n_zero*(2*nm+1) + m_plus) * (2*nl+1) +l_zero);

      int ind_RDImoins((0*n_dir+dir)*n_ind_lmn+ind_moins);
      int ind_RDIplus((0*n_dir+dir)*n_ind_lmn+ind_plus);
      force_flt[ind_RDImoins]=amplitude;
      force_flt[ind_RDIplus]=amplitude;


      force[0][dir][ind_moins]=amplitude;
      force[0][dir][ind_plus]=amplitude;
    }
  envoyer_broadcast(force_flt,0);
}


void Force_sp::compute_dirac_point_uniX_alongX()
{
  // Seul le proc 0 fait le calcul, comme c'est le mm clc pour tous les procs
  if (Process::je_suis_maitre())
    {
      // donne cos(x) (ex)
      int dir(0);
      int n_dir(3);
      int n_ind_lmn(n_lmn);
      int roc(momin+0);

      int l_moins(nl-roc);
      // int l_zero(nl);
      int l_plus(nl+roc);
      // int n_moins(nn-roc);
      int n_zero(nn);
      // int n_plus(nn+roc);
      // int m_moins(nm-roc);
      int m_zero(nm);
      // ints m_plus(nm+roc);

      int ind_moins((n_zero*(2*nm+1) + m_zero) * (2*nl+1) +l_moins);
      int ind_plus((n_zero*(2*nm+1) + m_zero) * (2*nl+1) +l_plus);
      std::cout << "ind_moins : " << ind_moins << std::endl;
      std::cout << "ind_plus : " << ind_plus << std::endl;
      std::cout << "in force_sp : " << Process::me() << std::endl;
      int ind_RDImoins((0*n_dir+dir)*n_ind_lmn+ind_moins);
      int ind_RDIplus((0*n_dir+dir)*n_ind_lmn+ind_plus);
      force_flt[ind_RDImoins]=amplitude;
      force_flt[ind_RDIplus]=amplitude;

      force[0][dir][ind_moins]=amplitude;
      force[0][dir][ind_plus]=amplitude;
    }
  envoyer_broadcast(force_flt,0);
}


void Force_sp::compute_dirac_point_uniX_alongY()
{
  // Seul le proc 0 fait le calcul, comme c'est le mm clc pour tous les procs
  if (Process::je_suis_maitre())
    {
      // donne cos(x) (ex)
      int dir(1);
      int n_dir(3);
      int n_ind_lmn(n_lmn);
      int roc(momin+0);

      int l_moins(nl-roc);
      // int l_zero(nl);
      int l_plus(nl+roc);
      // int n_moins(nn-roc);
      int n_zero(nn);
      // int n_plus(nn+roc);
      // int m_moins(nm-roc);
      int m_zero(nm);
      // ints m_plus(nm+roc);

      int ind_moins((n_zero*(2*nm+1) + m_zero) * (2*nl+1) +l_moins);
      int ind_plus((n_zero*(2*nm+1) + m_zero) * (2*nl+1) +l_plus);
      std::cout << "ind_moins : " << ind_moins << std::endl;
      std::cout << "ind_plus : " << ind_plus << std::endl;
      std::cout << "in force_sp : " << Process::me() << std::endl;
      int ind_RDImoins((0*n_dir+dir)*n_ind_lmn+ind_moins);
      int ind_RDIplus((0*n_dir+dir)*n_ind_lmn+ind_plus);
      force_flt[ind_RDImoins]=amplitude;
      force_flt[ind_RDIplus]=amplitude;

      force[0][dir][ind_moins]=amplitude;
      force[0][dir][ind_plus]=amplitude;
    }
  envoyer_broadcast(force_flt,0);
}




void Force_sp::compute_door_rope()
{
  // Seul le proc 0 fait le calcul, comme c'est le mm clc pour tous les procs
  if (Process::je_suis_maitre())
    {
      /*
      Fonction de test. Met le champ de force spectral à 1 sur un segment défini.
      Met le champ de force spectral à 0 ailleurs
      */
      int dir, l, m, n, ind;
      int roc(momin+0);

      // int l1(nl-roc), l2(nl+roc);
      int fsp_m1(nm-roc), fsp_m2(nm+roc);
      int n1(nn-roc), n2(nn+roc);

      int l1(nl), l2(nl);
      // int fsp_m1(nm), fsp_m2(nm);
      // int n1(nn), n2(nn);

      for (n=n1; n<n2+1; n++)
        for (m=fsp_m1; m<fsp_m2+1; m++)
          for (l=l1; l<l2+1; l++) // La TF-1 de cette force est réelle. Donc cette force est à symétrie Hermitienne
            {
              ind = (n*(2*nm+1) + m) * (2*nl+1) +l;
              for (dir=0; dir<3; dir++)
                {
                  force[0][dir][ind] = amplitude;
                  force[1][dir][ind] = 0.0;
                }
            }
    }
  envoyer_broadcast(force_flt,0);
}



void Force_sp::compute_door_cube()
{
  // Seul le proc 0 fait le calcul, comme c'est le mm clc pour tous les procs
  if (Process::je_suis_maitre())
    {
      int cpx, dir, l, m, n, ind;

      int roc(momin+0);

      int l1(nl-roc), l2(nl+roc);
      // int fsp_m1(nm-roc), fsp_m2(nm+roc);
      int n1(nn-roc), n2(nn+roc);

      // int l1(nl), l2(nl);
      int fsp_m1(nm), fsp_m2(nm);
      // int n1(nn), n2(nn);

      for (n=n1; n<n2+1; n++)
        for (m=fsp_m1; m<fsp_m2+1; m++)
          for (l=l1; l<l2+1; l++) // La TF-1 de cette force est réelle. Donc cette force est à symétrie Hermitienne
            {
              ind = (n*(2*nm+1) + m) * (2*nl+1) +l;
              for (dir=0; dir<3; dir++)
                {
                  for (cpx=0; cpx<2; cpx++)
                    {
                      force[0][dir][ind] = amplitude;
                      force[1][dir][ind] = 0;
                    }
                }
            }
    }
  envoyer_broadcast(force_flt,0);
}



void Force_sp::write(std::string nom_fichier_sortie, double t)
{
  // TODO :  il faut certainement mettre du envoyer_broadcast ou quoi
  std::ofstream Spectral_flux(nom_fichier_sortie.c_str(), std::ios::app);
  if (Spectral_flux)
    {
      int l, m, n, ind;
      double norme_kappa;
      double kappa[3];
      Spectral_flux << std::endl << "time : " << t << std::endl << std::endl;
      for (n=0; n<2*nn+1; n++)
        for (m=0; m<2*nm+1; m++)
          for (l=0; l<2*nl+1; l++)
            {
              ind = (n*(2*nm+1) + m) * (2*nl+1) +l;
              kappa[0] = - kmax + (l)*(2*kmax)/(2*nl);
              kappa[1] = - kmax + (m)*(2*kmax)/(2*nm);
              kappa[2] = - kmax + (n)*(2*kmax)/(2*nn);
              norme_kappa = sqrt(kappa[0]*kappa[0] + kappa[1]*kappa[1] + kappa[2]*kappa[2]);
              if (norme_kappa >= kmin)
                {
                  Spectral_flux << l<<","<<m<<","<<n<<"\t : ";
                  Spectral_flux << force[0][0][ind] << " + i" << force[1][0][ind]<<", \t";
                  Spectral_flux << force[0][1][ind] << " + i" << force[1][1][ind]<<", \t";
                  Spectral_flux << force[0][2][ind] << " + i" << force[1][2][ind]<<", \t";
                  Spectral_flux << std::endl;
                }
            }
    }
}


void Force_sp::write_separate(std::string nom_fichier_sortie, double t)
{
  std::ofstream Spectral_flux(nom_fichier_sortie.c_str());
  if (Spectral_flux)
    {
      double kappa[3];
      Spectral_flux << std::endl << "l,m,n, \t kx,ky,kz, \t  rf_x, \t cf_x, \t\t rf_y, \t cf_y, \t\t rf_z, \t cf_z \t";
      Spectral_flux << std::endl;

      for (int n=0; n<2*nn+1; n++)
        for (int m=0; m<2*nm+1; m++)
          {
            for (int l=0; l<2*nl+1; l++)
              {
                int ind = (n*(2*nm+1) + m) * (2*nl+1) +l;

                kappa[0] = - kmax + (l)*(2*kmax)/(2*nl);
                kappa[1] = - kmax + (m)*(2*kmax)/(2*nm);
                kappa[2] = - kmax + (n)*(2*kmax)/(2*nn);
                double norme_kappa = sqrt(kappa[0]*kappa[0] + kappa[1]*kappa[1] + kappa[2]*kappa[2]);
                if (norme_kappa >= kmin)
                  {
                    Spectral_flux << l<<","<<m<<","<<n<<"\t,";
                    Spectral_flux << kappa[0] << ",\t" << kappa[1] << ",\t" << kappa[2] <<", \t\t";
                    Spectral_flux << force[0][0][ind] << ",\t" << force[1][0][ind]<<", \t\t";
                    Spectral_flux << force[0][1][ind] << ",\t" << force[1][1][ind]<<", \t\t";
                    Spectral_flux << force[0][2][ind] << ",\t" << force[1][2][ind]<<"\t";
                    Spectral_flux << std::endl;
                  }
              }

          }
    }
}

void Force_sp::compute_energie()
{
  // Pas vérifié pour parallèle
  int l,m,n,dir,ind;
  energie = 0;
  for (n=0; n<2*nn+1; n++)
    for (m=0; m<2*nm+1; m++)
      for (l=0; l<2*nl+1; l++)
        {
          ind = (n*(2*nm+1) + m) * (2*nl+1) +l;
          for (dir=0; dir<3; dir++)
            energie += force[0][dir][ind]*force[0][dir][ind] + force[1][dir][ind]*force[1][dir][ind];
        }
}

double Force_sp::get_energie()
{
  return energie;
}

double Force_sp::get_force(int cpx, int dir, int ind)
{
  return force[cpx][dir][ind];
}

std::vector< std::vector< std:: vector <double >>> Force_sp::get_coeff()
{
  return force;
}


ArrOfDouble& Force_sp::get_coeff_flt()
{
  return force_flt;
}
