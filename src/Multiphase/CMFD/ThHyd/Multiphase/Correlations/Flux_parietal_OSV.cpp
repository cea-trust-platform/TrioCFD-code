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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Flux_parietal_OSV.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Flux_parietal_OSV.h>
#include <Flux_parietal_adaptatif.h>
#include <Loi_paroi_adaptative.h>
#include <Correlation.h>
#include <Pb_Multiphase.h>
#include <Domaine_dis.h>
#include <TRUSTTrav.h>
#include <Milieu_composite.h>
#include <Saturation_base.h>

#include <math.h>

Implemente_instanciable(Flux_parietal_OSV, "Flux_parietal_OSV", Flux_parietal_base);

Sortie& Flux_parietal_OSV::printOn(Sortie& os) const { return Flux_parietal_base::printOn(os); }

Entree& Flux_parietal_OSV::readOn(Entree& is)
{
  const Pb_Multiphase& pbm = ref_cast(Pb_Multiphase, pb_.valeur());
  correlation_monophasique_.typer_lire(pbm, "Flux_parietal", is);
  Cout << que_suis_je() << " : single-phase wall heat flux is " << correlation_monophasique_->que_suis_je() << finl;
  correlation_loi_paroi_ = ref_cast(Pb_Multiphase, pb_.valeur()).get_correlation("Loi_paroi");

  Param param(que_suis_je());
  //param.ajouter("Dh", &Dh_, Param::REQUIRED);
  param.ajouter("Water", &Water_);
  param.ajouter("R12", &R12_);
  param.ajouter("R134a", &R134a_);
  param.ajouter("R22", &R22_);
  param.lire_avec_accolades(is);

  if (Water_) fac_pression_ = 1, fac_puissance_ = 1.;
  else if (R12_) fac_pression_ = 5.9, fac_puissance_ = 0.06792763348649916; // Ahmad
  else if (R22_) fac_pression_ = 4.8, fac_puissance_ = 0.09411288750689908; // Ahmad
  else if (R134a_) fac_pression_=5.8, fac_puissance_ = 0.07353463409726364; // Ahmad
  else Process::exit(que_suis_je() + " : you must define a fluid for the simulation to run !");

  if ( !sub_type(Milieu_composite, pb_->milieu())) Process::exit(que_suis_je() + " : the medium must be composite !");
  if (!pbm.nom_phase(0).debute_par("liquide")) Process::exit(que_suis_je() + " : the first phase must be liquid !");

  const Milieu_composite& milc = ref_cast(Milieu_composite, pb_->milieu());
  int n_g = -1;
  for (int k = 1 ; k<pbm.nb_phases() ; k++)
    if (milc.has_saturation(0, k)) n_g += 1;
  if (n_g > 0) Process::exit(que_suis_je() + " : there can only be one evaporating phase for the carrying liquid for now ! Please feel free to update the code if you need.");

  return is;
}

void Flux_parietal_OSV::completer()
{
  correlation_monophasique_->completer();
}

void Flux_parietal_OSV::qp(const input_t& in, output_t& out) const
{
  // On met tout a 0 a tout hasard
  if (out.qpk)     (*out.qpk)    = 0.;
  if (out.da_qpk)  (*out.da_qpk) = 0.;
  if (out.dp_qpk)  (*out.dp_qpk) = 0.;
  if (out.dv_qpk)  (*out.dv_qpk) = 0.;
  if (out.dTf_qpk) (*out.dTf_qpk)= 0.;
  if (out.dTp_qpk) (*out.dTp_qpk)= 0.;
  if (out.qpi)     (*out.qpi)    = 0.;
  if (out.da_qpi)  (*out.da_qpi) = 0.;
  if (out.dp_qpi)  (*out.dp_qpi) = 0.;
  if (out.dv_qpi)  (*out.dv_qpi) = 0.;
  if (out.dTf_qpi) (*out.dTf_qpi)= 0.;
  if (out.dTp_qpi) (*out.dTp_qpi)= 0.;

  // On remplit le monophasique ; pas besoin du flux interfacial normalement
  ref_cast(Flux_parietal_base, correlation_monophasique_.valeur()).qp(in, out);

  // Ici la phase liquide est forcement la phase 0 car la correlation monophasique ne remplit que la phase 0
  int n_l = 0 ;
  const Milieu_composite& milc = ref_cast(Milieu_composite, pb_->milieu());

  if (out.nonlinear) (*out.nonlinear) = 1; // we turn on nonlinear in case it goes to two-phase during the newton

  for (int k = 0 ; k < in.N ; k++)
    if (n_l != k)
      if (milc.has_saturation(n_l, k))
        {
          // on part de ind_sat = (k*(in.N-1)-(k-1)*(k)/2) + (l-k-1) // avec k<l

          int ind_sat = k<n_l ? ( k *(in.N-1)-( k -1)*( k )/2) + (n_l- k -1) :
                        (n_l*(in.N-1)-(n_l-1)*(n_l)/2) + ( k -n_l-1);

          if (in.Tp - in.Tsat[ind_sat] > 0.) // Else : no wall superheat => no nucleation => single phase heat transfer only, we leave coefficients filled by single phase
            {
              // Flux diphasique ; on linearise car Tp commence a Tl qui peut etre proche de Tsat...
              double q_Jens_Lottes     = in.Tp-in.Tsat[ind_sat] > 1. ?
                                         fac_puissance_ * 1.e6 * std::pow( (in.Tp-in.Tsat[ind_sat])/25.*std::exp(in.p/62.e5 * fac_pression_), 4. ) :
                                         fac_puissance_ * 1.e6 * (in.Tp-in.Tsat[ind_sat]) * std::pow( 1./25.*std::exp(in.p/62.e5 * fac_pression_), 4. ) ;
              double dTp_q_Jens_Lottes = in.Tp-in.Tsat[ind_sat] > 1. ?
                                         fac_puissance_ * 1.e6 * 4.* std::pow(in.Tp-in.Tsat[ind_sat], 3) * std::pow(1./25.*std::exp(in.p/62.e5* fac_pression_), 4. ) :
                                         fac_puissance_ * 1.e6 * std::pow(1./25.*std::exp(in.p/62.e5* fac_pression_), 4. ) ;

              if ((*out.qpk)(n_l) < q_Jens_Lottes ) // other case : single phase greater than JL => single phase heat tranfer only
                {
                  //double h_mono = (*out.dTp_qpk)(n_l);

                  // On met tout a 0 avant remplissage de la partie diphasique
                  if (out.qpk)     (*out.qpk)    = 0.;
                  if (out.da_qpk)  (*out.da_qpk) = 0.;
                  if (out.dp_qpk)  (*out.dp_qpk) = 0.;
                  if (out.dv_qpk)  (*out.dv_qpk) = 0.;
                  if (out.dTf_qpk) (*out.dTf_qpk)= 0.;
                  if (out.dTp_qpk) (*out.dTp_qpk)= 0.;
                  if (out.qpi)     (*out.qpi)    = 0.;
                  if (out.da_qpi)  (*out.da_qpi) = 0.;
                  if (out.dp_qpi)  (*out.dp_qpi) = 0.;
                  if (out.dv_qpi)  (*out.dv_qpi) = 0.;
                  if (out.dTf_qpi) (*out.dTf_qpi)= 0.;
                  if (out.dTp_qpi) (*out.dTp_qpi)= 0.;

                  if (in.Tsat[ind_sat] <= in.T[n_l]) // Everything in the evaporative term
                    {
                      // Evaporation only
                      if (out.qpi)     (*out.qpi)(n_l, k)     = q_Jens_Lottes;
                      if (out.dTp_qpi) (*out.dTp_qpi)(n_l, k) = dTp_q_Jens_Lottes;
                    }

                  else // Local Saha Zuber thing
                    {

                      // Flux monophasique maximal
                      const Loi_paroi_base& corr_loi_paroi = ref_cast(Loi_paroi_base, correlation_loi_paroi_.valeur());
                      const double u_tau = corr_loi_paroi.get_utau(in.f);
                      double theta_plus = std::max( 1., beta_ + 2.12 * std::log( corr_loi_paroi.get_y_plus(in.f) ) );
                      double h_mono_max = in.rho[0] * in.Cp[0] * u_tau / theta_plus ; // Cf flux parietal adaptatif

                      double q_tot     = q_Jens_Lottes ;
                      double dTp_q_tot = dTp_q_Jens_Lottes ;
                      double dTl_q_tot = 0. ;

                      double q_evap     = ( q_tot - h_mono_max * (in.Tp-in.T[n_l]) > 0. ) ? q_tot - h_mono_max * (in.Tp-in.T[n_l]) : 0. ;
                      double dTp_q_evap = ( q_tot - h_mono_max * (in.Tp-in.T[n_l]) > 0. ) ? dTp_q_tot - h_mono_max : 0. ;
                      double dTl_q_evap = ( q_tot - h_mono_max * (in.Tp-in.T[n_l]) > 0. ) ? dTl_q_tot + h_mono_max : 0. ;

                      double q_mono = q_tot-q_evap ;
                      double dTp_q_mono = dTp_q_tot - dTp_q_evap ;
                      double dTl_q_mono = dTl_q_tot - dTl_q_evap ;

                      // Fill in the outputs : heat flux towards liquid
                      if (out.qpk)     (*out.qpk)(n_l)        = q_mono;
                      if (out.dTp_qpk) (*out.dTp_qpk)(n_l)    = dTp_q_mono;
                      if (out.dTf_qpk) (*out.dTf_qpk)(n_l,n_l)= dTl_q_mono;

                      // Evaporation
                      if (out.qpi)     (*out.qpi)(n_l, k)          = q_evap;
                      if (out.dTp_qpi) (*out.dTp_qpi)(n_l, k)      = dTp_q_evap;
                      if (out.dTf_qpi) (*out.dTf_qpi)(n_l, k, n_l) = dTl_q_evap;
                    }
                }
            }
        }
}

