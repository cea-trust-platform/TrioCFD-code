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
// File:        Algorithmes_Transport_FT_Disc.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/6
//
//////////////////////////////////////////////////////////////////////////////

#include <Algorithmes_Transport_FT_Disc.h>
#include <Parser.h>
#include <TRUSTArray.h>
#include <TRUSTTab.h>
#include <Zone_VF.h>

Implemente_instanciable(Algorithmes_Transport_FT_Disc,"Algorithmes_Transport_FT_Disc",Objet_U);


/*! @brief readOn : appel interdit
 *
 */
Entree& Algorithmes_Transport_FT_Disc::readOn(Entree& is)
{
  assert(0);
  return is;
}

/*! @brief printOn : appel interdit
 *
 */
Sortie& Algorithmes_Transport_FT_Disc::printOn(Sortie& os) const
{
  assert(0);
  return os;
}

/*! @brief Calcul des grandeurs suivantes en fonction de l'indicatrice de phase.
 *
 * On calcule la somme sur tous les processeurs.
 *
 * @param (indicatrice) le tableau des valeurs d'un champ aux elements de zone_vf (entree) Le champ est suppose constant par element.
 * @param (volume_total) (sortie) le volume total de la zone = INTEGRALE(d omega)
 * @param (volume_phase_1) (sortie) integrale de volume de l'indicatrice de phase, egale au volume de la phase 1 = INTEGRALE(indicatrice * d omega)
 * @param (barycentre) En entree: on attend un tableau de meme dimension que les coordonnees En sortie: iso-barycentre de la zone = INTEGRALE(X * d omega) / volume_phase_1
 * @param (barycentre_phase_1) En entree: idem que barycentre En sortie: barycentre de la phase_1 calcule comme INTEGRALE(indicatrice * X * d omega) / volume_phase_1
 */
void Algorithmes_Transport_FT_Disc::calculer_moments_indicatrice(
  const Zone_VF& zone_vf,
  const DoubleTab& indicatrice,
  double& volume_total,
  double& volume_phase_1,
  ArrOfDouble& barycentre,
  ArrOfDouble& barycentre_phase_1)
{
  const DoubleTab& xp = zone_vf.xp();
  const int dim = xp.line_size();
  assert(barycentre.size_array() == dim);
  assert(barycentre_phase_1.size_array() == dim);
  // Volumes des elements
  const DoubleVect& volumes = zone_vf.volumes();
  const int nb_elem = indicatrice.size();
  assert(nb_elem == zone_vf.nb_elem());
  int i;
  double somme_v = 0.; // Somme des volumes
  double somme_vI = 0.; // Somme des volumes de phase 1
  double somme_v_x[3] = { 0., 0., 0. }; // Somme des xp[i] * volume
  double somme_vI_x[3] = { 0., 0., 0. }; // Somme des xp[i] * volume * indicatrice
  for (i = 0; i < nb_elem; i++)
    {
      const double I = indicatrice[i];
      const double v = volumes[i];
      const double vI = v * I;
      somme_v += v;
      somme_vI += vI;
      int j;
      for (j = 0; j < dim; j++)
        {
          const double x = xp(i, j);
          somme_v_x[j] += v * x;
          somme_vI_x[j] += vI * x;
        }
    }
  volume_total = mp_sum(somme_v);
  volume_phase_1 = mp_sum(somme_vI);
  for (i = 0; i < dim; i++)
    {
      const double sv = mp_sum(somme_v_x[i]);
      const double svI = mp_sum(somme_vI_x[i]);
      double b = 0.;
      double bI = 0.;
      if (volume_total > 0.)
        b = sv / volume_total;
      if (volume_phase_1 > 0.)
        bI = svI / volume_phase_1;
      barycentre[i] = b;
      barycentre_phase_1[i] = bI;
    }
}

/*! @brief evalue les fonctions px, py et pz en (x,y,z,t) (essentiellement utilise dans integrer_vitesse_imposee)
 *
 * @param (px, py, pz) Trois fonctions qui comprennent au moins 4 parametres.
 * @param (vx, vy, vz) References aux variables ou on stocke le resultat vx=px(x,y,z,t) vy=py(x,y,z,t) vz=pz(x,y,z,t)
 */
void eval_vitesse(double x, double y, double z, double t,
                  Parser& px, Parser& py, Parser& pz,
                  double& vx, double& vy, double& vz)
{
  int i0=0, i1=1, i2=2, i3=3;
  px.setVar(i0,x);
  px.setVar(i1,y);
  px.setVar(i2,z);
  px.setVar(i3,t);

  vx = px.eval();
  py.setVar(i0,x);
  py.setVar(i1,y);
  py.setVar(i2,z);
  py.setVar(i3,t);

  vy = py.eval();
  pz.setVar(i0,x);
  pz.setVar(i1,y);
  pz.setVar(i2,z);
  pz.setVar(i3,t);
  vz = pz.eval();
}

/*! @brief Integre le systeme differentiel d/dt(x)=vx(x,y,z,t)
 *
 *    d/dt(y)=vy(x,y,z,t)
 *    d/dt(z)=vy(x,y,z,t)
 *   entre "t=temps" et "t=temps+dt" par un unique pas de Runge Kutta ordre 4
 *  Parametres: x,y,z
 *  Signification: Position initiale (en t=temps) en entree, position finale
 *   (en t=temps+dt) en sortie.
 *
 */
void integrer_vitesse_imposee(
  Parser& parser_vx, Parser& parser_vy, Parser& parser_vz,
  double temps, double dt, double& x, double& y, double& z)
{
  // Runge Kutta ordre 3:
  double vx, vy, vz;
  eval_vitesse(x,
               y,
               z,
               temps,
               parser_vx,parser_vy,parser_vz,vx,vy,vz);
  double vx1,vy1,vz1;
  eval_vitesse(x + vx * dt * 0.5,
               y + vy * dt * 0.5,
               z + vz * dt * 0.5,
               temps + dt * 0.5,
               parser_vx,parser_vy,parser_vz,vx1,vy1,vz1);
  double vx2,vy2,vz2;
  eval_vitesse(x + vx1 * dt * 0.5,
               y + vy1 * dt * 0.5,
               z + vz1 * dt * 0.5,
               temps + dt * 0.5,
               parser_vx,parser_vy,parser_vz,vx2,vy2,vz2);
  double vx3,vy3,vz3;
  eval_vitesse(x + vx2 * dt,
               y + vy2 * dt,
               z + vz2 * dt,
               temps + dt,
               parser_vx,parser_vy,parser_vz,vx3,vy3,vz3);

  x += (vx + 2.*vx1 + 2.*vx2 + vx3) / 6. * dt;
  y += (vy + 2.*vy1 + 2.*vy2 + vy3) / 6. * dt;
  z += (vz + 2.*vz1 + 2.*vz2 + vz3) / 6. * dt;
}

