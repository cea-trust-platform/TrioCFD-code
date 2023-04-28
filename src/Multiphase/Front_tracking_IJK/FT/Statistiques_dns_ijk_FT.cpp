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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Statistiques_dns_ijk_FT.cpp
// Directory : $IJK_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////
#include <Statistiques_dns_ijk_FT.h>
#include <IJK_Grid_Geometry.h>
#include <TRUSTTab.h>
#include <communications.h>
#include <Param.h>
#include <IJK_Navier_Stokes_tools.h>
#include <IJK_FT.h>

#define CREATE_STR(sum_id) #sum_id    // turn the enum value into a string.

// #define STAT_VERBOSE
// #define STAT_EXIT
#define PRECISION_DERIVEES (100.1) // (1e-10)
#define PRECISION_DIVU (1e-10)
#define PRECISION_DDIVU (1e-10)
#define PRECISION_LAPLP (1e-10)

Implemente_instanciable_sans_constructeur(Statistiques_dns_ijk_FT,"Statistiques_dns_ijk_FT",Statistiques_dns_ijk);
#include <Statistiques_defines.h>
#include <Statistiques_defines_thermique.h>
// Pour generer la liste XXX_MOY, envoyer la liste des #define dans
//  sed -e '/\/\//d' -e 's/#define /"/g;s/_MOY.*/",/g' /tmp/titi >/tmp/resultat
// et copier ci-dessous:
// Ou inversement :
// sed -e "s/\"//g" -e "s/,/_MOY/g" list  | awk '{a=NR+476;print "#define " $0 " " a}' > list_define

Statistiques_dns_ijk_FT::Statistiques_dns_ijk_FT():
  check_stats_(0),
  nb_thermal_fields_(0),
  nvalt_(0)
{
  Cerr << "The ref to the correct IJK_FT is filled later" << finl;
}

Statistiques_dns_ijk_FT::Statistiques_dns_ijk_FT(IJK_FT_double& ijk_ft):
  check_stats_(0),
  nb_thermal_fields_(0),
  nvalt_(0),
  ref_ijk_ft_(ijk_ft)
{
// I : Indicatrice.
// Iv : Indicatrice vapeur (1. - chi)
  const char *noms_moyennes_prov[] =
  {
    "I",
    "UI",
    "VI",
    "WI",
    "PI",
    "UIv",
    "VIv",
    "WIv",
    "PIv",
    "UUI",
    "VVI",
    "WWI",
    "UVI",
    "VWI",
    "UWI",
    "UUUI",
    "UUVI",
    "UUWI",
    "UVVI",
    "UVWI",
    "UWWI",
    "VVVI",
    "VVWI",
    "VWWI",
    "WWWI",
    "UPI",
    "VPI",
    "WPI",
    "UdIdx",
    "VdIdx",
    "WdIdx",
    "UdIdy",
    "VdIdy",
    "WdIdy",
    "UdIdz",
    "VdIdz",
    "WdIdz",
    "UPdIdx",
    "VPdIdx",
    "WPdIdx",
    "UPdIdy",
    "VPdIdy",
    "WPdIdy",
    "UPdIdz",
    "VPdIdz",
    "WPdIdz",
    "UUdIdx",
    "VVdIdx",
    "WWdIdx",
    "UUdIdy",
    "VVdIdy",
    "WWdIdy",
    "UUdIdz",
    "VVdIdz",
    "WWdIdz",
    "UVdIdx",
    "UWdIdx",
    "VWdIdx",
    "UVdIdy",
    "UWdIdy",
    "VWdIdy",
    "UVdIdz",
    "UWdIdz",
    "VWdIdz",
    "PdIdx",
    "PdIdy",
    "PdIdz",
    "IdIdx",
    "IdIdy",
    "IdIdz",
    "IdUdx",
    "IdVdx",
    "IdWdx",
    "IdUdy",
    "IdVdy",
    "IdWdy",
    "IdUdz",
    "IdVdz",
    "IdWdz",
    "IPdUdx",
    "IPdVdx",
    "IPdWdx",
    "IPdUdy",
    "IPdVdy",
    "IPdWdy",
    "IPdUdz",
    "IPdVdz",
    "IPdWdz",
    "dUdxdIdx",
    "dVdxdIdx",
    "dWdxdIdx",
    "dUdydIdy",
    "dVdydIdy",
    "dWdydIdy",
    "dUdzdIdz",
    "dVdzdIdz",
    "dWdzdIdz",
    "UdUdxdIdx",
    "UdVdxdIdx",
    "UdWdxdIdx",
    "UdUdydIdy",
    "UdVdydIdy",
    "UdWdydIdy",
    "UdUdzdIdz",
    "UdVdzdIdz",
    "UdWdzdIdz",
    "VdUdxdIdx",
    "VdVdxdIdx",
    "VdWdxdIdx",
    "VdUdydIdy",
    "VdVdydIdy",
    "VdWdydIdy",
    "VdUdzdIdz",
    "VdVdzdIdz",
    "VdWdzdIdz",
    "WdUdxdIdx",
    "WdVdxdIdx",
    "WdWdxdIdx",
    "WdUdydIdy",
    "WdVdydIdy",
    "WdWdydIdy",
    "WdUdzdIdz",
    "WdVdzdIdz",
    "WdWdzdIdz",
    "IdUdxdUdx",
    "IdUdxdUdy",
    "IdUdxdUdz",
    "IdUdxdVdx",
    "IdUdxdVdy",
    "IdUdxdVdz",
    "IdUdxdWdx",
    "IdUdxdWdy",
    "IdUdxdWdz",
    "IdUdydUdy",
    "IdUdydUdz",
    "IdUdydVdx",
    "IdUdydVdy",
    "IdUdydVdz",
    "IdUdydWdx",
    "IdUdydWdy",
    "IdUdydWdz",
    "IdUdzdUdz",
    "IdUdzdVdx",
    "IdUdzdVdy",
    "IdUdzdVdz",
    "IdUdzdWdx",
    "IdUdzdWdy",
    "IdUdzdWdz",
    "IdVdxdVdx",
    "IdVdxdVdy",
    "IdVdxdVdz",
    "IdVdxdWdx",
    "IdVdxdWdy",
    "IdVdxdWdz",
    "IdVdydVdy",
    "IdVdydVdz",
    "IdVdydWdx",
    "IdVdydWdy",
    "IdVdydWdz",
    "IdVdzdVdz",
    "IdVdzdWdx",
    "IdVdzdWdy",
    "IdVdzdWdz",
    "IdWdxdWdx",
    "IdWdxdWdy",
    "IdWdxdWdz",
    "IdWdydWdy",
    "IdWdydWdz",
    "IdWdzdWdz",
    "IdPdx",
    "IdPdy",
    "IdPdz",
    "IUdPdx",
    "IUdPdy",
    "IUdPdz",
    "IVdPdx",
    "IVdPdy",
    "IVdPdz",
    "IWdPdx",
    "IWdPdy",
    "IWdPdz",
    "IUdUdx",
    "IUdVdx",
    "IUdWdx",
    "IUdUdy",
    "IUdVdy",
    "IUdWdy",
    "IUdUdz",
    "IUdVdz",
    "IUdWdz",
    "IVdUdx",
    "IVdVdx",
    "IVdWdx",
    "IVdUdy",
    "IVdVdy",
    "IVdWdy",
    "IVdUdz",
    "IVdVdz",
    "IVdWdz",
    "IWdUdx",
    "IWdVdx",
    "IWdWdx",
    "IWdUdy",
    "IWdVdy",
    "IWdWdy",
    "IWdUdz",
    "IWdVdz",
    "IWdWdz",
    "IddUdxx",
    "IddVdxx",
    "IddWdxx",
    "IddUdyy",
    "IddVdyy",
    "IddWdyy",
    "IddUdzz",
    "IddVdzz",
    "IddWdzz",
    "Fx",
    "Fy",
    "Fz",
    "Frx",
    "Fry",
    "Frz",
    "Frax",
    "Fray",
    "Fraz",
    "DISSIP",
    "UUIv",
    "VVIv",
    "WWIv",
    "UVIv",
    "VWIv",
    "UWIv",
    "UUIbv",
    "VVIbv",
    "WWIbv",
    "UUIb",
    "VVIb",
    "WWIb",
    "UUUIb",
    "UUVIb",
    "UUWIb",
    "UVVIb",
    "UWWIb",
    "VVVIb",
    "VVWIb",
    "VWWIb",
    "WWWIb",
    "AI",
    "NK",
    "IdUdxdUdxdUdx",
    "IdUdydUdydUdx",
    "IdUdzdUdzdUdx",
    "IdUdxdVdxdUdy",
    "IdUdydVdydUdy",
    "IdUdzdVdzdUdy",
    "IdUdxdWdxdUdz",
    "IdUdydWdydUdz",
    "IdUdzdWdzdUdz",
    "IdVdxdUdxdVdx",
    "IdVdydUdydVdx",
    "IdVdzdUdzdVdx",
    "IdVdxdVdxdVdy",
    "IdVdydVdydVdy",
    "IdVdzdVdzdVdy",
    "IdVdxdWdxdVdz",
    "IdVdydWdydVdz",
    "IdVdzdWdzdVdz",
    "IdWdxdUdxdWdx",
    "IdWdydUdydWdx",
    "IdWdzdUdzdWdx",
    "IdWdxdVdxdWdy",
    "IdWdydVdydWdy",
    "IdWdzdVdzdWdy",
    "IdWdxdWdxdWdz",
    "IdWdydWdydWdz",
    "IdWdzdWdzdWdz",
    "IdUdxdUdxW",
    "IdUdydUdyW",
    "IdUdzdUdzW",
    "IdVdxdVdxW",
    "IdVdydVdyW",
    "IdVdzdVdzW",
    "IdWdxdWdxW",
    "IdWdydWdyW",
    "IdWdzdWdzW",
    "IdUdxddPdxdx",
    "IdUdyddPdxdy",
    "IdUdzddPdxdz",
    "IdVdxddPdydx",
    "IdVdyddPdydy",
    "IdVdzddPdydz",
    "IdWdxddPdzdx",
    "IdWdyddPdzdy",
    "IdWdzddPdzdz",
    "IddUdxdxddUdxdx",
    "IddUdxdyddUdxdy",
    "IddUdxdzddUdxdz",
    "IddUdydxddUdydx",
    "IddUdydyddUdydy",
    "IddUdydzddUdydz",
    "IddUdzdxddUdzdx",
    "IddUdzdyddUdzdy",
    "IddUdzdzddUdzdz",
    "IddVdxdxddVdxdx",
    "IddVdxdyddVdxdy",
    "IddVdxdzddVdxdz",
    "IddVdydxddVdydx",
    "IddVdydyddVdydy",
    "IddVdydzddVdydz",
    "IddVdzdxddVdzdx",
    "IddVdzdyddVdzdy",
    "IddVdzdzddVdzdz",
    "IddWdxdxddWdxdx",
    "IddWdxdyddWdxdy",
    "IddWdxdzddWdxdz",
    "IddWdydxddWdydx",
    "IddWdydyddWdydy",
    "IddWdydzddWdydz",
    "IddWdzdxddWdzdx",
    "IddWdzdyddWdzdy",
    "IddWdzdzddWdzdz",
    "dIdxddUdxdxdUdx",
    "dIdxddUdxdydUdy",
    "dIdxddUdxdzdUdz",
    "dIdyddUdydxdUdx",
    "dIdyddUdydydUdy",
    "dIdyddUdydzdUdz",
    "dIdzddUdzdxdUdx",
    "dIdzddUdzdydUdy",
    "dIdzddUdzdzdUdz",
    "dIdxddVdxdxdVdx",
    "dIdxddVdxdydVdy",
    "dIdxddVdxdzdVdz",
    "dIdyddVdydxdVdx",
    "dIdyddVdydydVdy",
    "dIdyddVdydzdVdz",
    "dIdzddVdzdxdVdx",
    "dIdzddVdzdydVdy",
    "dIdzddVdzdzdVdz",
    "dIdxddWdxdxdWdx",
    "dIdxddWdxdydWdy",
    "dIdxddWdxdzdWdz",
    "dIdyddWdydxdWdx",
    "dIdyddWdydydWdy",
    "dIdyddWdydzdWdz",
    "dIdzddWdzdxdWdx",
    "dIdzddWdzdydWdy",
    "dIdzddWdzdzdWdz",
    "dIdxddUdxdz",
    "dIdyddUdydz",
    "dIdzddUdzdz",
    "dIdzdUdxdUdx",
    "dIdzdUdydUdy",
    "dIdzdUdzdUdz",
    "dIdzdVdxdVdx",
    "dIdzdVdydVdy",
    "dIdzdVdzdVdz",
    "dIdzdWdxdWdx",
    "dIdzdWdydWdy",
    "dIdzdWdzdWdz",
    "Nx",
    "Ny",
    "Nz",
    "UaiNx",
    "VaiNx",
    "WaiNx",
    "UaiNy",
    "VaiNy",
    "WaiNy",
    "UaiNz",
    "VaiNz",
    "WaiNz",
    "UPaiNx",
    "VPaiNx",
    "WPaiNx",
    "UPaiNy",
    "VPaiNy",
    "WPaiNy",
    "UPaiNz",
    "VPaiNz",
    "WPaiNz",
    "UUaiNx",
    "VVaiNx",
    "WWaiNx",
    "UUaiNy",
    "VVaiNy",
    "WWaiNy",
    "UUaiNz",
    "VVaiNz",
    "WWaiNz",
    "UVaiNx",
    "UWaiNx",
    "VWaiNx",
    "UVaiNy",
    "UWaiNy",
    "VWaiNy",
    "UVaiNz",
    "UWaiNz",
    "VWaiNz",
    "PaiNx",
    "PaiNy",
    "PaiNz",
    "IaiNx",
    "IaiNy",
    "IaiNz",
    "dUdxaiNx",
    "dVdxaiNx",
    "dWdxaiNx",
    "dUdyaiNy",
    "dVdyaiNy",
    "dWdyaiNy",
    "dUdzaiNz",
    "dVdzaiNz",
    "dWdzaiNz",
    "UdUdxaiNx",
    "UdVdxaiNx",
    "UdWdxaiNx",
    "UdUdyaiNy",
    "UdVdyaiNy",
    "UdWdyaiNy",
    "UdUdzaiNz",
    "UdVdzaiNz",
    "UdWdzaiNz",
    "VdUdxaiNx",
    "VdVdxaiNx",
    "VdWdxaiNx",
    "VdUdyaiNy",
    "VdVdyaiNy",
    "VdWdyaiNy",
    "VdUdzaiNz",
    "VdVdzaiNz",
    "VdWdzaiNz",
    "WdUdxaiNx",
    "WdVdxaiNx",
    "WdWdxaiNx",
    "WdUdyaiNy",
    "WdVdyaiNy",
    "WdWdyaiNy",
    "WdUdzaiNz",
    "WdVdzaiNz",
    "WdWdzaiNz",
    "aiNxddUdxdxdUdx",
    "aiNxddUdxdydUdy",
    "aiNxddUdxdzdUdz",
    "aiNyddUdydydUdy",
    "aiNyddUdydzdUdz",
    "aiNzddUdzdzdUdz",
    "aiNxddVdxdxdVdx",
    "aiNxddVdxdydVdy",
    "aiNxddVdxdzdVdz",
    "aiNyddVdydydVdy",
    "aiNyddVdydzdVdz",
    "aiNzddVdzdzdVdz",
    "aiNxddWdxdxdWdx",
    "aiNxddWdxdydWdy",
    "aiNxddWdxdzdWdz",
    "aiNyddWdydydWdy",
    "aiNyddWdydzdWdz",
    "aiNzddWdzdzdWdz",
    "aiNxddUdxdz",
    "aiNyddUdydz",
    "aiNzddUdzdz",
    "aiNzdUdxdUdx",
    "aiNzdUdydUdy",
    "aiNzdUdzdUdz",
    "aiNzdVdxdVdx",
    "aiNzdVdydVdy",
    "aiNzdVdzdVdz",
    "aiNzdWdxdWdx",
    "aiNzdWdydWdy",
    "aiNzdWdzdWdz",
    "aiNx",
    "aiNy",
    "aiNz",
    "kaiNx",
    "kaiNy",
    "kaiNz",
    "aiNxx",
    "aiNxy",
    "aiNxz",
    "aiNyy",
    "aiNyz",
    "aiNzz",
    "aiNNxxxx",
    "aiNNxxxy",
    "aiNNxxxz",
    "aiNNxxyy",
    "aiNNxxyz",
    "aiNNxxzz",
    "aiNNxyyy",
    "aiNNxyyz",
    "aiNNxyzz",
    "aiNNxzyy",
    "aiNNxzyz",
    "aiNNxzzz",
    "aiNNyyyy",
    "aiNNyyyz",
    "aiNNyyzz",
    "aiNNyzzz",
    "aiNNzzzz",
    "kai",
// };
//  const char *noms_moyennes_prov_fixe[] =
// {
    "UI_INT",
    "VI_INT",
    "WI_INT",
    "UI_INTUI_INT",
    "UI_INTVI_INT",
    "UI_INTWI_INT",
    "VI_INTVI_INT",
    "VI_INTWI_INT",
    "WI_INTWI_INT",
    "UI_INTUI_INTUI_INT",
    "UI_INTUI_INTVI_INT",
    "UI_INTUI_INTWI_INT",
    "UI_INTVI_INTUI_INT",
    "UI_INTVI_INTVI_INT",
    "UI_INTVI_INTWI_INT",
    "UI_INTWI_INTUI_INT",
    "UI_INTWI_INTVI_INT",
    "UI_INTWI_INTWI_INT",
    "VI_INTUI_INTUI_INT",
    "VI_INTUI_INTVI_INT",
    "VI_INTUI_INTWI_INT",
    "VI_INTVI_INTUI_INT",
    "VI_INTVI_INTVI_INT",
    "VI_INTVI_INTWI_INT",
    "VI_INTWI_INTUI_INT",
    "VI_INTWI_INTVI_INT",
    "VI_INTWI_INTWI_INT",
    "WI_INTUI_INTUI_INT",
    "WI_INTUI_INTVI_INT",
    "WI_INTUI_INTWI_INT",
    "WI_INTVI_INTUI_INT",
    "WI_INTVI_INTVI_INT",
    "WI_INTVI_INTWI_INT",
    "WI_INTWI_INTUI_INT",
    "WI_INTWI_INTVI_INT",
    "WI_INTWI_INTWI_INT",
    "UI_INTUIUI",
    "UI_INTUIVI",
    "UI_INTUIWI",
    "UI_INTVIUI",
    "UI_INTVIVI",
    "UI_INTVIWI",
    "UI_INTWIUI",
    "UI_INTWIVI",
    "UI_INTWIWI",
    "VI_INTUIUI",
    "VI_INTUIVI",
    "VI_INTUIWI",
    "VI_INTVIUI",
    "VI_INTVIVI",
    "VI_INTVIWI",
    "VI_INTWIUI",
    "VI_INTWIVI",
    "VI_INTWIWI",
    "WI_INTUIUI",
    "WI_INTUIVI",
    "WI_INTUIWI",
    "WI_INTVIUI",
    "WI_INTVIVI",
    "WI_INTVIWI",
    "WI_INTWIUI",
    "WI_INTWIVI",
    "WI_INTWIWI",
    "VI_INTP_INT",
    "WI_INTP_INT",
    "dUINTdxdUINTdx",
    "dUINTdxdUINTdy",
    "dUINTdxdUINTdz",
    "dUINTdydUINTdy",
    "dUINTdydUINTdz",
    "dUINTdzdUINTdz",
    "dUINTdxdVINTdx",
    "dUINTdxdVINTdy",
    "dUINTdxdVINTdz",
    "dUINTdydVINTdy",
    "dUINTdydVINTdz",
    "dUINTdzdVINTdz",
    "dUINTdxdWINTdx",
    "dUINTdxdWINTdy",
    "dUINTdxdWINTdz",
    "dUINTdydWINTdy",
    "dUINTdydWINTdz",
    "dUINTdzdWINTdz",
    "dVINTdxdUINTdx",
    "dVINTdxdUINTdy",
    "dVINTdxdUINTdz",
    "dVINTdydUINTdy",
    "dVINTdydUINTdz",
    "dVINTdzdUINTdz",
    "dVINTdxdVINTdx",
    "dVINTdxdVINTdy",
    "dVINTdxdVINTdz",
    "dVINTdydVINTdy",
    "dVINTdydVINTdz",
    "dVINTdzdVINTdz",
    "dVINTdxdWINTdx",
    "dVINTdxdWINTdy",
    "dVINTdxdWINTdz",
    "dVINTdydWINTdy",
    "dVINTdydWINTdz",
    "dVINTdzdWINTdz",
    "dWINTdxdUINTdx",
    "dWINTdxdUINTdy",
    "dWINTdxdUINTdz",
    "dWINTdydUINTdy",
    "dWINTdydUINTdz",
    "dWINTdzdUINTdz",
    "dWINTdxdVINTdx",
    "dWINTdxdVINTdy",
    "dWINTdxdVINTdz",
    "dWINTdydVINTdy",
    "dWINTdydVINTdz",
    "dWINTdzdVINTdz",
    "dWINTdxdWINTdx",
    "dWINTdxdWINTdy",
    "dWINTdxdWINTdz",
    "dWINTdydWINTdy",
    "dWINTdydWINTdz",
    "dWINTdzdWINTdz",
    "P_INTdUINTdx",
    "P_INTdUINTdy",
    "P_INTdUINTdz",
    "P_INTdVINTdx",
    "P_INTdVINTdy",
    "P_INTdVINTdz",
    "P_INTdWINTdx",
    "P_INTdWINTdy",
    "P_INTdWINTdz",
    "dissip_int",
    "P_NOPERTURBE",
    "U_NOPERTURBE",
    "V_NOPERTURBE",
    "W_NOPERTURBE",
    "I_NP",
    "UI_INTP_INT",
    // pression
    "PNP",
    "PIUNP",
    "PIVNP",
    "PIWNP",
    "IPdUdxNP",
    "IPdVdxNP",
    "IPdWdxNP",
    "IPdUdyNP",
    "IPdVdyNP",
    "IPdWdyNP",
    "IPdUdzNP",
    "IPdVdzNP",
    "IPdWdzNP",
    "IdPdxNP",
    "IdPdyNP",
    "IdPdzNP",
    "PIseuil",
    "PIUseuil",
    "PIVseuil",
    "PIWseuil",
    "IPdUdxseuil",
    "IPdVdxseuil",
    "IPdWdxseuil",
    "IPdUdyseuil",
    "IPdVdyseuil",
    "IPdWdyseuil",
    "IPdUdzseuil",
    "IPdVdzseuil",
    "IPdWdzseuil",
    "IdPdxseuil",
    "IdPdyseuil",
    "IdPdzseuil",
    "PNPseuil",
    "PIUNPseuil",
    "PIVNPseuil",
    "PIWNPseuil",
    "IPdUdxNPseuil",
    "IPdVdxNPseuil",
    "IPdWdxNPseuil",
    "IPdUdyNPseuil",
    "IPdVdyNPseuil",
    "IPdWdyNPseuil",
    "IPdUdzNPseuil",
    "IPdVdzNPseuil",
    "IPdWdzNPseuil",
    "IdPdxNPseuil",
    "IdPdyNPseuil",
    "IdPdzNPseuil",
    "Ivrappelx",
    "Ivrappely",
    "Ivrappelz",
    "IP_INT",
    "UIvrappelx",
    "VIvrappelx",
    "WIvrappelx",
    "UIvrappely",
    "VIvrappely",
    "WIvrappely",
    "UIvrappelz",
    "VIvrappelz",
    "WIvrappelz",

    // termini aggiuntivi per il calcolo della parte deviatorica del tensore degli sforzi
    "dUdyaiNx",
    "dUdzaiNx",
    "dVdxaiNy",
    "dVdzaiNy",
    "dWdxaiNz",
    "dWdyaiNz",

    // Viscous stress in the vapor phase
    "dUdxIv",
    "dUdyIv",
    "dUdzIv",
    "dVdxIv",
    "dVdyIv",
    "dVdzIv",
    "dWdxIv",
    "dWdyIv",
    "dWdzIv",

    // Attempt for the molecular diffusion
    "IddUdxdxU",
    "IddVdxdxU",
    "IddWdxdxU",
    "IddUdxdxV",
    "IddVdxdxV",
    "IddWdxdxV",
    "IddUdxdxW",
    "IddVdxdxW",
    "IddWdxdxW",
    "IddUdydyU",
    "IddVdydyU",
    "IddWdydyU",
    "IddUdydyV",
    "IddVdydyV",
    "IddWdydyV",
    "IddUdydyW",
    "IddVdydyW",
    "IddWdydyW",
    "IddUdzdzU",
    "IddVdzdzU",
    "IddWdzdzU",
    "IddUdzdzV",
    "IddVdzdzV",
    "IddWdzdzV",
    "IddUdzdzW",
    "IddVdzdzW",
    "IddWdzdzW",

    // To deal with the numerical pressure ..
    "Ix",
    "xaiNx",
    "xaiNy",
    "xaiNz",
    "UIx",
    "VIx",
    "WIx",
// Redistribution corrective term
    "dUdxIx",
    "dUdyIx",
    "dUdzIx",
    "dVdxIx",
    "dVdyIx",
    "dVdzIx",
    "dWdxIx",
    "dWdyIx",
    "dWdzIx",
// Interfacial terms
    "xUaiNx",
    "xVaiNx",
    "xWaiNx",
    "xUaiNy",
    "xVaiNy",
    "xWaiNy",
    "xUaiNz",
    "xVaiNz",
    "xWaiNz",
    //extended pressure

    "P_VAP_Iv",
    "P_VAP_aiNx",
    "P_VAP_aiNy",
    "P_VAP_aiNz",
    "P_LIQ_I",
    "P_LIQ_aiNx",
    "P_LIQ_aiNy",
    "P_LIQ_aiNz",
    "UP_LIQ_aiNx",
    "VP_LIQ_aiNx",
    "WP_LIQ_aiNx",
    "UP_LIQ_aiNy",
    "VP_LIQ_aiNy",
    "WP_LIQ_aiNy",
    "UP_LIQ_aiNz",
    "VP_LIQ_aiNz",
    "WP_LIQ_aiNz",
    "IP_LIQ_dUdx",
    "IP_LIQ_dUdy",
    "IP_LIQ_dUdz",
    "IP_LIQ_dVdx",
    "IP_LIQ_dVdy",
    "IP_LIQ_dVdz",
    "IP_LIQ_dWdx",
    "IP_LIQ_dWdy",
    "IP_LIQ_dWdz",
    "UP_LIQ_I",
    "VP_LIQ_I",
    "WP_LIQ_I",

    // GAB DISSIP
    "TRUE_DISSIP_MOY",
    "DISSIP_VAP_MOY",
    "TRUE_DISSIP_VAP_MOY",
    //"IdP_LIQ_dx",
    //"IdP_LIQ_dy",
    //"IdP_LIQ_dz",
    //"IUdP_LIQ_dx",
    //"IUdP_LIQ_dy",
    //"IUdP_LIQ_dz",
    //"IVdP_LIQ_dx",
    //"IVdP_LIQ_dy",
    //"IVdP_LIQ_dz",
    //"IWdP_LIQ_dx",
    //"IWdP_LIQ_dy",
    //"IWdP_LIQ_dz"

  };

  // int nb_val_1=502;
  //int nb_val_2=196; // Bulles fixes.
  //nval_=nb_val_1+nb_val_2;

  // FONCTIONNE
//  nval_=502+196+27+6+9+25+29;//+6;
  // A VOIR
  nval_=502+196+27+6+9+25+29+3;

  noms_moyennes_.dimensionner_force(nval_);
  for (int i=0; i<nval_; i++)
    noms_moyennes_[i]=noms_moyennes_prov[i];

  const char *noms_moyennes_provt[] =
  {
    "TI",
    "TTI",
    "ITU",
    "ITV",
    "ITW",
    "ITUU",
    "ITUV",
    "ITUW",
    "ITVV",
    "ITVW",
    "ITWW",
    "ITdPdx",
    "ITdPdy",
    "ITdPdz",
    "ITddUdxdx",
    "ITddUdydy",
    "ITddUdzdz",
    "ITddVdxdx",
    "ITddVdydy",
    "ITddVdzdz",
    "ITddWdxdx",
    "ITddWdydy",
    "ITddWdzdz",
    "IUddTdxdx",
    "IUddTdydy",
    "IUddTdzdz",
    "IVddTdxdx",
    "IVddTdydy",
    "IVddTdzdz",
    "IWddTdxdx",
    "IWddTdydy",
    "IWddTdzdz",
    "ITbulles"
  };
  nvalt_ = 33;
  noms_moyennes_temperature_.dimensionner_force(nvalt_);
  for (int i=0; i<nvalt_; i++)
    noms_moyennes_temperature_[i]=noms_moyennes_provt[i];

}

// Les champs doivent tous etre sur le domaine ns normalement...
// y compris ceux-ci :
//    o  field_ai
//    o  field_kappa_ai
//    o  normale_cell
// On pourrait avoir des champs ft tant  qu'on ne les multiplie pas avec des champs de splitting different.
// Mais ce n'est pas le cas car on veut multiplier par P et U...

double Calculer_valeur_seuil(double p_seuil, double pressionijk, double p_seuil_max, double p_seuil_min)
{
  if (pressionijk>p_seuil_max)
    {
      p_seuil = p_seuil_max;
    }
  else if (pressionijk<p_seuil_min)
    {
      p_seuil = p_seuil_min;
    }
  else
    {
      p_seuil = pressionijk;
    }
  return p_seuil;
}

void Statistiques_dns_ijk_FT::update_stat(IJK_FT_double& cas, const double dt)
{

  FixedVector<IJK_Field_double, 3>& vitesse=cas.velocity_;
  IJK_Field_double& extended_pressure_liq= (cas.post_.extended_pressure_computed_) ? cas.post_.extended_pl_: cas.pressure_;
  IJK_Field_double& extended_pressure_vap= (cas.post_.extended_pressure_computed_) ? cas.post_.extended_pv_: cas.pressure_;
  IJK_Field_double& pression=cas.pressure_;
  const IJK_Field_double& indicatrice=cas.itfce().I();
  FixedVector<IJK_Field_double, 3>& gradP=cas.post_.grad_P_;

  // Nombre total de mailles en K
  const int nktot = pression.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);

// For each slice in the domain, the integral is divided by the length of the domain
// facteur 1./(ni*nj) car sommation de ni*nj valeurs sur des mailles de meme taille
// OU facteur delta_z / taille_totale_en_z  si mailles non uniformes en z
  const int nijtot = pression.get_splitting().get_grid_geometry().get_nb_elem_tot(DIRECTION_I)
                     * pression.get_splitting().get_grid_geometry().get_nb_elem_tot(DIRECTION_J);
  double facteur = 1./(double)(nijtot);

  // Nombre local de mailles en K
  const int imax = pression.ni();
  const int jmax = pression.nj();
  const int kmax = pression.nk();
  const int offset = pression.get_splitting().get_offset_local(DIRECTION_K);
  const int periok = pression.get_splitting().get_grid_geometry().get_periodic_flag(DIRECTION_K);

  const bool dipha = (!cas.disable_diphasique_);
  IJK_Field_local_double zero, zeros;
  zero.allocate(imax, jmax, kmax, 1 /* ghost pour le grad */, 1. /* additional_layers */, 1 /*nb compo*/);
  zeros.allocate(imax, jmax, kmax, 1 /* ghost pour le grad */, 1. /* additional_layers */, 3 /*nb compo*/);
  //allocate_velocity()

  FixedVector<IJK_Field_double, 3>& gradI=cas.post_.grad_I_ns_;
  FixedVector<IJK_Field_double, 3>& sourceI=cas.terme_source_interfaces_ns_;
  FixedVector<IJK_Field_double, 3>& repuls=cas.terme_repulsion_interfaces_ns_;
  FixedVector<IJK_Field_double, 3>& absrepuls=cas.terme_abs_repulsion_interfaces_ns_;
  IJK_Field_double& field_ai=cas.post_.ai_ns_;
  IJK_Field_double& field_kappa_ai=cas.post_.kappa_ai_ns_;
  FixedVector<IJK_Field_double, 3>& normale_cell=cas.post_.normale_cell_ns_;
  IJK_Field_double& field_dudx=cas.post_.dudx_ ;
  IJK_Field_double& field_dvdy=cas.post_.dvdy_ ;
  IJK_Field_double& field_dwdx=cas.post_.dwdx_ ;
  IJK_Field_double& field_dudz=cas.post_.dudz_ ;
  IJK_Field_double& field_dvdz=cas.post_.dvdz_ ;
  IJK_Field_double& field_dwdz=cas.post_.dwdz_ ;
  const Motcles& liste_post_instantanes = ref_ijk_ft_->post_.get_liste_post_instantanes();
  FixedVector<IJK_Field_double, 3>& rot=cas.post_.rot_ ;
  IJK_Field_double& critere_Q=cas.post_.critere_Q_;

  double coef_immobilisation_=cas.coef_immobilisation_;
  IJK_Field_double& indicatrice_np=cas.post_.indicatrice_non_perturbe_;
  FixedVector<IJK_Field_double, 3>& integral_vitesse=cas.post_.integrated_velocity_;
  IJK_Field_double& integral_pression=cas.post_.integrated_pressure_;
  FixedVector<IJK_Field_double, 3> force_rappel=cas.force_rappel_;
  double p_seuil_max = cas.p_seuil_max_;
  double p_seuil_min = cas.p_seuil_min_;
  const IJK_Field_double& integral_vitesse_i = integral_vitesse[0];
  const IJK_Field_double& integral_vitesse_j = integral_vitesse[1];
  const IJK_Field_double& integral_vitesse_k = integral_vitesse[2];

  ArrOfDouble coord_i, coord_j, coord_k;
  build_local_coords(pression, coord_i, coord_j, coord_k);

  if (elem_coord_.size_array() == 0)
    {
      Cerr << "Erreur dans Statistiques_dns_ijk::update_stat: non initialise" << finl;
      Process::exit();
    }

  const IJK_Field_double& vitesse_i = vitesse[0];
  const IJK_Field_double& vitesse_j = vitesse[1];
  const IJK_Field_double& vitesse_k = vitesse[2];

  // Calcul le gradient de U aux cellules a partir de la vitesse aux faces :
  compute_and_store_gradU_cell(vitesse_i, vitesse_j, vitesse_k,
                               /* Et les champs en sortie */
                               field_dudx, field_dvdy, field_dwdx,
                               field_dudz, field_dvdz, field_dwdz);

#ifdef STAT_VERBOSE
  // Pour verifications :
  double erreur, erreur_ddP, divU, laplP, d_divU_dx, d_divU_dy, d_divU_dz;
  double L1_erreur=0., L1_erreur_ddP=0., L1_divU=0., L1_laplP=0., L1_d_divU_dx=0., L1_d_divU_dy=0., L1_d_divU_dz=0.;
  double L2_erreur=0., L2_erreur_ddP=0., L2_divU=0., L2_laplP=0., L2_d_divU_dx=0., L2_d_divU_dy=0., L2_d_divU_dz=0.;
#endif

  DoubleTab tmp(nktot, nval_);

  //Invalid value for the interpolation
  const double  p_invalid = 1.e19;
  //Number of invalid cells on a single slice
  IntTab occ(nktot,2); // 2 cols : 0:liq, 1:vap
  occ=0;
  DoubleTab tmpt(nktot, nvalt_,nb_thermal_fields_);
  //Here starts the loop in space
  for (int k = 0; k < kmax; k++)
    {
      const double dz = tab_dz_[k+offset];
      bool on_the_first_cell = false;
      bool on_the_last_cell = false;
      const int kglob=k + offset;
      if (!periok)
        {
          if (kglob == 0)
            on_the_first_cell = true;
          if (kglob == nktot - 1)
            on_the_last_cell = true;
        }
      // Calcul des moyennes spatiales sur le plan ij

      // On y stocke la somme de toutes les valeurs sur le plan ij, on divisera apres
      ArrOfDouble moy(nval_);
      for (int i = 0; i < nval_; i++)
        {
          moy[i] = 0.;
        }
      DoubleTab moyv(nvalt_,nb_thermal_fields_);
      moyv = 0.; // Table initialization.
      for (int j = 0; j < jmax; j++)
        {
          for (int i = 0; i < imax; i++)
            {
              double x = coord_i[i]; //This function find the x-coordinate of the assigned point of the plane
              // Vitesses au centre de la maille i,j,k
              // On peut ici interpoler plus finement si on veut:
              double u = (vitesse_i(i,j,k) + vitesse_i(i+1, j, k)) * 0.5;
              double v = (vitesse_j(i,j,k) + vitesse_j(i, j+1, k)) * 0.5;
              double w = (vitesse_k(i,j,k) + vitesse_k(i, j, k+1)) * 0.5;
              double p = pression(i,j,k);
              double pl = extended_pressure_liq(i,j,k);
              if  (extended_pressure_liq(i,j,k) > p_invalid)
                {
                  pl = 0;
                  occ(kglob,0) +=1;
                }
              double pv = extended_pressure_vap(i,j,k);
              if  (extended_pressure_vap(i,j,k) > p_invalid)
                {
                  pv = 0;
                  occ(kglob,1) +=1;
                }
              double chi = indicatrice(i,j,k);
              double chiv = 1. - chi;

              // Gradient de pression au centre de la maille i,j,k
              double dPdx = (gradP[0](i,j,k) + gradP[0](i+1, j, k)) * 0.5;
              double dPdy = (gradP[1](i,j,k) + gradP[1](i, j+1, k)) * 0.5;
              double dPdz = (gradP[2](i,j,k) + gradP[2](i, j, k+1)) * 0.5;

              // Produit vitesse par gradient de pression au centre de la maille i,j,k
              double UdPdx = (vitesse_i(i,j,k)*gradP[0](i,j,k) + vitesse_i(i+1, j, k)*gradP[0](i+1, j, k)) * 0.5;
              double VdPdy = (vitesse_j(i,j,k)*gradP[1](i,j,k) + vitesse_j(i, j+1, k)*gradP[1](i, j+1, k)) * 0.5;
              double WdPdz = (vitesse_k(i,j,k)*gradP[2](i,j,k) + vitesse_k(i, j, k+1)*gradP[2](i, j, k+1)) * 0.5;

              // Produits vitesse/vitesse au centre de la maille i,j,k
              // (pour la partie symetrique du tenseur uu, on regarde une autre interpolation) :
              // ((ui+uj)/2)**2  != (ui**2 + uj**2)/2
              double uu = (vitesse_i(i,j,k)*vitesse_i(i,j,k) + vitesse_i(i+1, j, k)*vitesse_i(i+1, j, k)) * 0.5;
              double vv = (vitesse_j(i,j,k)*vitesse_j(i,j,k) + vitesse_j(i, j+1, k)*vitesse_j(i, j+1, k)) * 0.5;
              double ww = (vitesse_k(i,j,k)*vitesse_k(i,j,k) + vitesse_k(i, j, k+1)*vitesse_k(i, j, k+1)) * 0.5;

              double uuu = (vitesse_i(i,j,k)*vitesse_i(i,j,k)*vitesse_i(i,j,k) + vitesse_i(i+1, j, k)*vitesse_i(i+1, j, k)*vitesse_i(i+1, j, k)) * 0.5;
              double vvv = (vitesse_j(i,j,k)*vitesse_j(i,j,k)*vitesse_j(i,j,k) + vitesse_j(i, j+1, k)*vitesse_j(i, j+1, k)*vitesse_j(i, j+1, k)) * 0.5;
              double www = (vitesse_k(i,j,k)*vitesse_k(i,j,k)*vitesse_k(i,j,k) + vitesse_k(i, j, k+1)*vitesse_k(i, j, k+1)*vitesse_k(i, j, k+1)) * 0.5;

              // Gradient de vitesse au centre de la maille i,j,k
              double  dUdx=0., dVdx=0., dWdx=0., dUdy=0., dVdy=0., dWdy=0., dUdz=0., dVdz=0., dWdz=0.;
              double  ddUdxdx=0., ddVdxdx=0., ddWdxdx=0., ddUdydy=0., ddVdydy=0., ddWdydy=0., ddUdzdz=0., ddVdzdz=0., ddWdzdz=0.;
              double pseudo_dissip = calculer_gradients_vitesse(vitesse_i, vitesse_j, vitesse_k,
                                                                i,j,k,
                                                                dz,
                                                                dUdx, dVdx, dWdx,
                                                                dUdy, dVdy, dWdy,
                                                                dUdz, dVdz, dWdz,
                                                                ddUdxdx, ddVdxdx, ddWdxdx,
                                                                ddUdydy, ddVdydy, ddWdydy,
                                                                ddUdzdz, ddVdzdz, ddWdzdz,
                                                                on_the_first_cell, on_the_last_cell);
              double true_dissip = calculer_vraie_dissipation(pseudo_dissip,
                                                              dUdx, dUdy, dUdz,
                                                              dVdx, dVdy, dVdz,
                                                              dWdx, dWdy, dWdz);

              // Pour verifier un peu les stats ou pour le post-traitement, on stocke (systematiquement) les grad de vitesse :
              //if ((check_stats_) || (postraitement des gradients demande?))
              //{
              gradU_[0](i,j,k) = dUdx;
              gradU_[1](i,j,k) = dUdy;
              gradU_[2](i,j,k) = dUdz;
              gradV_[0](i,j,k) = dVdx;
              gradV_[1](i,j,k) = dVdy;
              gradV_[2](i,j,k) = dVdz;
              gradW_[0](i,j,k) = dWdx;
              gradW_[1](i,j,k) = dWdy;
              gradW_[2](i,j,k) = dWdz;
              //}
#ifdef STAT_VERBOSE
              // Debug pour verifier que je m'emmele pas les pinceaux dans les variables :
              erreur = std::fabs(field_dudx(i,j,k)  - dUdx)
                       + std::fabs(field_dvdy(i,j,k)  - dVdy)
                       + std::fabs(field_dwdx(i,j,k)  - dWdx)
                       + std::fabs(field_dudz(i,j,k)  - dUdz)
                       + std::fabs(field_dvdz(i,j,k)  - dVdz)
                       + std::fabs(field_dwdz(i,j,k)  - dWdz);

              if (erreur > PRECISION_DERIVEES)
                {
                  Cerr << "Statistiques_dns_ijk_FT::update_stat -- Erreur (gradU) superieure au seuil : "
                       << i << " " << j << " " << k << " "
                       << erreur
                       << finl;
                  Cerr << "dudx: " << field_dudx(i,j,k) << " " << dUdx << " " << field_dudx(i,j,k)  - dUdx << finl;
                  Cerr << "dvdy: " << field_dvdy(i,j,k) << " " << dVdy << " " << field_dvdy(i,j,k)  - dVdy << finl;
                  Cerr << "dwdx: " << field_dwdx(i,j,k) << " " << dWdx << " " << field_dwdx(i,j,k)  - dWdx << finl;
                  Cerr << "dudz: " << field_dudz(i,j,k) << " " << dUdz << " " << field_dudz(i,j,k)  - dUdz << finl;
                  Cerr << "dvdz: " << field_dvdz(i,j,k) << " " << dVdz << " " << field_dvdz(i,j,k)  - dVdz << finl;
                  Cerr << "dwdz: " << field_dwdz(i,j,k) << " " << dWdz << " " << field_dwdz(i,j,k)  - dWdz << finl;
                  Cerr << "On proc " << Process::me() << " Cell: " << i << " " << j << " " << k << finl;
#ifdef STAT_EXIT
                  Process::exit();
#endif
                }

              divU = dUdx + dVdy + dWdz;
              if (std::fabs(divU) > PRECISION_DIVU)
                {
                  Cerr << "Statistiques_dns_ijk_FT::update_stat -- divU : " << divU << finl;
                  Cerr << "dUdx,dVdy,dWdz : " << dUdx << " " <<  dVdy << " " << dWdz << finl;
                  Cerr << "On proc " << Process::me() << " Cell: " << i << " " << j << " " << k << finl;
#ifdef STAT_EXIT
                  Process::exit();
#endif
                }

              // Compute L1 and L2 norms :
              L1_erreur += erreur;
              L2_erreur += erreur*erreur;
              L1_divU += std::fabs(divU);
              L2_divU += divU*divU;
#endif

              // Second gradients croises :
              double ddUdxdy = 0.;
              double ddUdxdz = 0.;
              //double ddUdydx = ddUdxdy;

              double ddUdydz = 0.;
              //double ddUdzdx = ddUdxdz;
              //double ddUdzdy = ddUdydz;

              double ddVdxdy = 0.;
              double ddVdxdz = 0.;
              //double ddVdydx = ddVdxdy;

              double ddVdydz = 0.;
              //double ddVdzdx = ddVdxdz;
              //double ddVdzdy = ddVdydz;

              double ddWdxdy = 0.;
              double ddWdxdz = 0.;
              //double ddWdydx = ddWdxdy;

              double ddWdydz = 0.;
              //double ddWdzdx = ddWdxdz;
              //double ddWdzdy = ddWdydz;
              cell_to_cell_gradient(i,j,k,
                                    field_dudx, field_dvdy, field_dwdx,
                                    field_dudz, field_dvdz, field_dwdz,
                                    /* Et les outputs en ref aussi!! */
                                    ddUdxdy, ddUdxdz, ddUdydz,
                                    ddVdxdy, ddVdxdz, ddVdydz,
                                    ddWdxdy, ddWdxdz, ddWdydz);

              // Pour verifier un peu les stats, on peut stocker le double grad :
              if (check_stats_)
                {
                  grad2Ui_[0](i,j,k) = ddUdxdx;
                  grad2Ui_[1](i,j,k) = ddUdydy;
                  grad2Ui_[2](i,j,k) = ddUdzdz;
                  // Et les croisees :
                  grad2Uc_[0](i,j,k) = ddUdxdy;
                  grad2Uc_[1](i,j,k) = ddUdxdz;
                  grad2Uc_[2](i,j,k) = ddUdydz;
                  //
                  grad2Vi_[0](i,j,k) = ddVdxdx;
                  grad2Vi_[1](i,j,k) = ddVdydy;
                  grad2Vi_[2](i,j,k) = ddVdzdz;
                  // Et les croisees :
                  grad2Vc_[0](i,j,k) = ddVdxdy;
                  grad2Vc_[1](i,j,k) = ddVdxdz;
                  grad2Vc_[2](i,j,k) = ddVdydz;
                  //
                  grad2Wi_[0](i,j,k) = ddWdxdx;
                  grad2Wi_[1](i,j,k) = ddWdydy;
                  grad2Wi_[2](i,j,k) = ddWdzdz;
                  // Et les croisees :
                  grad2Wc_[0](i,j,k) = ddWdxdy;
                  grad2Wc_[1](i,j,k) = ddWdxdz;
                  grad2Wc_[2](i,j,k) = ddWdydz;
                  //
                }

              // Les croises non calcules sont en fait symetriques :
              double ddUdydx = ddUdxdy;
              double ddUdzdx = ddUdxdz;
              double ddUdzdy = ddUdydz;
              double ddVdydx = ddVdxdy;
              double ddVdzdx = ddVdxdz;
              double ddVdzdy = ddVdydz;
              double ddWdydx = ddWdxdy;
              double ddWdzdx = ddWdxdz;
              double ddWdzdy = ddWdydz;
#ifdef STAT_VERBOSE
              d_divU_dx = ddUdxdx + ddVdydx + ddWdzdx;
              d_divU_dy = ddUdxdy + ddVdydy + ddWdzdy;
              d_divU_dz = ddUdxdz + ddVdydz + ddWdzdz;
              if ((std::fabs(d_divU_dx) > PRECISION_DDIVU)
                  || (std::fabs(d_divU_dy) > PRECISION_DDIVU)
                  || (std::fabs(d_divU_dz) > PRECISION_DDIVU) )
                {
                  Cerr << "Statistiques_dns_ijk_FT::update_stat -- grad(divU) : "
                       << d_divU_dx << " "
                       << d_divU_dy << " "
                       << d_divU_dz << " "
                       << finl;
                  Cerr << "ddUdxdz,ddVdydz,ddWdzdz : " << ddUdxdz << " " << ddVdydz << " " << ddWdzdz << finl;
                  Cerr << "On proc " << Process::me() << " Cell: " << i << " " << j << " " << k << finl;
#ifdef STAT_EXIT
                  Process::exit();
#endif
                }

              // Compute L1 and L2 norms :
              L1_d_divU_dx += std::fabs(d_divU_dx);
              L2_d_divU_dx += d_divU_dx*d_divU_dx;
              L1_d_divU_dy += std::fabs(d_divU_dy);
              L2_d_divU_dy += d_divU_dy*d_divU_dy;
              L1_d_divU_dz += std::fabs(d_divU_dz);
              L2_d_divU_dz += d_divU_dz*d_divU_dz;
#endif
              // Si CURL ou CRITERE_Q sont demandes, on les remplis.
              if (liste_post_instantanes.contient_("CURL"))
                {
                  rot[0](i,j,k) = dWdy - dVdz;
                  rot[1](i,j,k) = dUdz - dWdx;
                  rot[2](i,j,k) = dVdx - dUdy;
                }

              if (liste_post_instantanes.contient_("CRITERE_Q"))
                {
                  // Calcul du critere Q selon (Jeong & Hussain 1995)
                  critere_Q(i,j,k) = -0.5 * (dUdx*dUdx + 2.*dUdy*dVdx + 2.*dUdz*dWdx
                                             + dVdy*dVdy + 2.*dVdz*dWdy
                                             + dWdz*dWdz );
                }

              // Derivee seconde de la pression :
              double ddPdxdx = 0.;
              double ddPdxdy = 0.;
              double ddPdxdz = 0.;
              double ddPdydx = 0.;
              double ddPdydy = 0.;
              double ddPdydz = 0.;
              double ddPdzdx = 0.;
              double ddPdzdy = 0.;
              double ddPdzdz = 0.;
              face_to_cell_gradient(gradP[0],gradP[1],gradP[2],
                                    i,j,k,
                                    dz,
                                    ddPdxdx,ddPdydx,ddPdzdx,
                                    ddPdxdy,ddPdydy,ddPdzdy,
                                    ddPdxdz,ddPdydz,ddPdzdz,
                                    on_the_first_cell, on_the_last_cell,
                                    1 /* bc_type for gradP at the wall :
				       	   	   	   	   assumed equal to faces values */);

              // Pour verifier un peu les stats, on peut stocker le double grad :
              if (check_stats_)
                {
                  grad2Pi_[0](i,j,k) = ddPdxdx;
                  grad2Pi_[1](i,j,k) = ddPdydy;
                  grad2Pi_[2](i,j,k) = ddPdzdz;
                  // Et les croisees :
                  grad2Pc_[0](i,j,k) = ddPdxdy;
                  grad2Pc_[1](i,j,k) = ddPdxdz;// ou ddPdzdx (sont equivalents d'ordre 2)
                  grad2Pc_[2](i,j,k) = ddPdzdy;// ou ddPdzdy (sont equivalents d'ordre 2)
                }

#ifdef STAT_VERBOSE
              // Debug pour verifier que je m'emmele pas les pinceaux dans les variables :
              erreur_ddP = std::fabs(ddPdydx  - ddPdxdy)
                           + std::fabs(ddPdzdx - ddPdxdz)
                           + std::fabs(ddPdzdy - ddPdydz);

              if (erreur_ddP > PRECISION_DERIVEES)
                {
                  Cerr << "Statistiques_dns_ijk_FT::update_stat -- Erreur_ddP superieure au seuil : "
                       << i << " " << j << " " << k << " "
                       << erreur_ddP
                       << finl;
                  Cerr << "We should have symmetric behaviour for crossed derivatives : \n"
                       << "ddPdydx: " << ddPdydx << " " << ddPdxdy << " diff: " << ddPdydx-ddPdxdy << "\n";
                  Cerr << "ddPdzdx: " << ddPdzdx << " " << ddPdxdz << " diff: " << ddPdzdx-ddPdxdz << "\n"
                       << "ddPdydz: " << ddPdydz << " " << ddPdzdy << " diff: " << ddPdydz-ddPdzdy << "\n"
                       << finl;
                  Cerr << "On proc " << Process::me() << " Cell: " << i << " " << j << " " << k << finl;
                  Process::exit();
                }

              laplP = ddPdxdx + ddPdydy + ddPdzdz;
              if (std::fabs(dIdx)+std::fabs(dIdy)+std::fabs(dIdz) < PRECISION_LAPLP )
                {
                  // On a gradI nul sur toutes les faces donc la cellule est "franchement" monophasique.
                  // On ne devrait donc pas avoir de surprise, rho aux faces devrait etre constant
                  // et par consequent, on devrait avoir : laplP = 0
                  if (std::fabs(laplP) > PRECISION_LAPLP)
                    {
                      Cerr << "Statistiques_dns_ijk_FT::update_stat -- Laplacien P : " << laplP << finl;
                      Cerr << "On proc " << Process::me() << " Cell: " << i << " " << j << " " << k << finl;
#ifdef STAT_EXIT
                      Process::exit();
#endif
                    }
                }

              // Compute L1 and L2 norms :
              L1_erreur_ddP += erreur_ddP;
              L2_erreur_ddP += erreur_ddP*erreur_ddP;
              L1_laplP += std::fabs(laplP);
              L2_laplP += laplP*laplP;
#endif

              double ai=0., kai=0., Nx=0.,Ny=0.,Nz=0., aiNx=0.,aiNy=0.,aiNz=0., Nxx=0.,Nxy=0.,Nxz=0.;
              double Frx=0.,Fry=0.,Frz=0., Frax=0.,Fray=0.,Fraz=0., Fx=0.,Fy=0.,Fz=0., Nyy=0.,Nyz=0.,Nzz=0.;
              double dIdx=0., dIdy=0., dIdz=0.;
              if (dipha)
                {
                  ai = field_ai(i,j,k);
                  kai = field_kappa_ai(i,j,k);

                  // Normale a l'interface au centre de l'elem :
                  Nx = normale_cell[0](i,j,k);
                  Ny = normale_cell[1](i,j,k);
                  Nz = normale_cell[2](i,j,k);

                  aiNx = ai*Nx;
                  aiNy = ai*Ny;
                  aiNz = ai*Nz;

                  Nxx = Nx*Nx;
                  Nxy = Nx*Ny;
                  Nxz = Nx*Nz;
                  Nyy = Ny*Ny;
                  Nyz = Ny*Nz;
                  Nzz = Nz*Nz;

                  // Force de repulsion :
                  Fx = (sourceI[0](i,j,k) + sourceI[0](i+1, j, k)) * 0.5;
                  Fy = (sourceI[1](i,j,k) + sourceI[1](i, j+1, k)) * 0.5;
                  Fz = (sourceI[2](i,j,k) + sourceI[2](i, j, k+1)) * 0.5;

                  // Force de repulsion :
                  Frx = (repuls[0](i,j,k) + repuls[0](i+1, j, k)) * 0.5;
                  Fry = (repuls[1](i,j,k) + repuls[1](i, j+1, k)) * 0.5;
                  Frz = (repuls[2](i,j,k) + repuls[2](i, j, k+1)) * 0.5;

                  // Force de repulsion :
                  Frax = (absrepuls[0](i,j,k) + absrepuls[0](i+1, j, k)) * 0.5;
                  Fray = (absrepuls[1](i,j,k) + absrepuls[1](i, j+1, k)) * 0.5;
                  Fraz = (absrepuls[2](i,j,k) + absrepuls[2](i, j, k+1)) * 0.5;

                  // Gradient de l'indicatrice au centre de la maille i,j,k
                  dIdx = (gradI[0](i,j,k) + gradI[0](i+1, j, k)) * 0.5;
                  dIdy = (gradI[1](i,j,k) + gradI[1](i, j+1, k)) * 0.5;
                  dIdz = (gradI[2](i,j,k) + gradI[2](i, j, k+1)) * 0.5;
                }

              int idx =0;
              for (auto& itr : cas.thermique_)
                {
                  const IJK_Field_double& temperature = itr.get_temperature();
                  const double T = temperature(i,j,k);

                  double T_adim_bulles;
                  if (liste_post_instantanes.contient_("TEMPERATURE_ADIM_BULLES"))
                    {
                      const IJK_Field_double& temperature_adim_bulles = itr.get_temperature_adim_bulles();
                      T_adim_bulles = temperature_adim_bulles(i,j,k);
                    }
                  else
                    {
                      T_adim_bulles = 0;
                    }


                  // Derivee seconde de la temperature :
                  FixedVector<IJK_Field_double, 3>& gradT = itr.get_gradient_temperature();
                  double ddTdxdx = 0.;
                  double ddTdxdy = 0.;
                  double ddTdxdz = 0.;
                  double ddTdydx = 0.;
                  double ddTdydy = 0.;
                  double ddTdydz = 0.;
                  double ddTdzdx = 0.;
                  double ddTdzdy = 0.;
                  double ddTdzdz = 0.;
                  face_to_cell_gradient(gradT[0],gradT[1],gradT[2],
                                        i,j,k,
                                        dz,
                                        ddTdxdx,ddTdydx,ddTdzdx,
                                        ddTdxdy,ddTdydy,ddTdzdy,
                                        ddTdxdz,ddTdydz,ddTdzdz,
                                        on_the_first_cell, on_the_last_cell,
                                        1 /* bc_type for gradT at the wall :
									    assumed equal to faces values */);
#define AJOUTT(somme,val) moyv(somme,idx) += val
                  AJOUTT(TI_MOY,chi*T);
                  AJOUTT(TTI_MOY, chi*T*T);
                  AJOUTT(ITU_MOY, chi*T*u);
                  AJOUTT(ITV_MOY, chi*T*v);
                  AJOUTT(ITW_MOY, chi*T*w);
                  AJOUTT(ITUU_MOY, chi*T*u*u);
                  AJOUTT(ITUV_MOY, chi*T*u*v);
                  AJOUTT(ITUW_MOY, chi*T*u*w);
                  AJOUTT(ITVV_MOY, chi*T*v*v);
                  AJOUTT(ITVW_MOY, chi*T*v*w);
                  AJOUTT(ITWW_MOY, chi*T*w*w);
                  AJOUTT(ITdPdx_MOY, chi*T*dPdx);
                  AJOUTT(ITdPdy_MOY, chi*T*dPdy);
                  AJOUTT(ITdPdz_MOY, chi*T*dPdz);
                  AJOUTT(ITddUdxdx_MOY, chi*T*ddUdxdx);
                  AJOUTT(ITddUdydy_MOY, chi*T*ddUdydy);
                  AJOUTT(ITddUdzdz_MOY, chi*T*ddUdzdz);
                  AJOUTT(ITddVdxdx_MOY, chi*T*ddVdxdx);
                  AJOUTT(ITddVdydy_MOY, chi*T*ddVdydy);
                  AJOUTT(ITddVdzdz_MOY, chi*T*ddVdzdz);
                  AJOUTT(ITddWdxdx_MOY, chi*T*ddWdxdx);
                  AJOUTT(ITddWdydy_MOY, chi*T*ddWdydy);
                  AJOUTT(ITddWdzdz_MOY, chi*T*ddWdzdz);
                  AJOUTT(IUddTdxdx_MOY, chi*u*ddTdxdx);
                  AJOUTT(IUddTdydy_MOY, chi*u*ddTdydy);
                  AJOUTT(IUddTdzdz_MOY, chi*u*ddTdzdz);
                  AJOUTT(IVddTdxdx_MOY, chi*v*ddTdxdx);
                  AJOUTT(IVddTdydy_MOY, chi*v*ddTdydy);
                  AJOUTT(IVddTdzdz_MOY, chi*v*ddTdzdz);
                  AJOUTT(IWddTdxdx_MOY, chi*w*ddTdxdx);
                  AJOUTT(IWddTdydy_MOY, chi*w*ddTdydy);
                  AJOUTT(IWddTdzdz_MOY, chi*w*ddTdzdz);
                  AJOUTT(ITbulles_MOY, chi*T_adim_bulles);
#undef AJOUTT
                  idx++;
                }

#define AJOUT(somme,val) moy[somme] += val
              // moyennes
              AJOUT(I_MOY,chi);
              AJOUT(UI_MOY,u*chi);
              AJOUT(VI_MOY,v*chi);
              AJOUT(WI_MOY,w*chi);
              AJOUT(PI_MOY,p*chi);
              AJOUT(UIv_MOY,u*chiv);
              AJOUT(VIv_MOY,v*chiv);
              AJOUT(WIv_MOY,w*chiv);
              AJOUT(PIv_MOY,p*chiv);
              // correlations 2 eme ordre vitesse :
              AJOUT(UUI_MOY,u*u*chi);
              AJOUT(VVI_MOY,v*v*chi);
              AJOUT(WWI_MOY,w*w*chi);
              AJOUT(UVI_MOY,u*v*chi);
              AJOUT(VWI_MOY,v*w*chi);
              AJOUT(UWI_MOY,u*w*chi);
              // 3 eme ordre (10 termes)
              AJOUT(UUUI_MOY,u*u*u*chi);
              AJOUT(UUVI_MOY,u*u*v*chi);
              AJOUT(UUWI_MOY,u*u*w*chi);
              AJOUT(UVVI_MOY,u*v*v*chi);
              AJOUT(UVWI_MOY,u*v*w*chi);
              AJOUT(UWWI_MOY,u*w*w*chi);
              AJOUT(VVVI_MOY,v*v*v*chi);
              AJOUT(VVWI_MOY,v*v*w*chi);
              AJOUT(VWWI_MOY,v*w*w*chi);
              AJOUT(WWWI_MOY,w*w*w*chi);
              // Correlation vitesse/pression :
              AJOUT(UPI_MOY,u*p*chi);
              AJOUT(VPI_MOY,v*p*chi);
              AJOUT(WPI_MOY,w*p*chi);
              // Maintenant, ce qui fait intervenir les derivees :
              // 1. de l'indicatrice  :
              //^^^^^^^^^^^^^^^^^^^^^^^
              // ui*dI/dxb : (9)
              AJOUT(UdIdx_MOY,u*dIdx);
              AJOUT(VdIdx_MOY,v*dIdx);
              AJOUT(WdIdx_MOY,w*dIdx);
              AJOUT(UdIdy_MOY,u*dIdy);
              AJOUT(VdIdy_MOY,v*dIdy);
              AJOUT(WdIdy_MOY,w*dIdy);
              AJOUT(UdIdz_MOY,u*dIdz);
              AJOUT(VdIdz_MOY,v*dIdz);
              AJOUT(WdIdz_MOY,w*dIdz);
              //^^^^^^^^^^^^^^^^^^^^^^^
              // ui*p*dI/dxb : (9)
              AJOUT(UPdIdx_MOY,u*p*dIdx);
              AJOUT(VPdIdx_MOY,v*p*dIdx);
              AJOUT(WPdIdx_MOY,w*p*dIdx);
              AJOUT(UPdIdy_MOY,u*p*dIdy);
              AJOUT(VPdIdy_MOY,v*p*dIdy);
              AJOUT(WPdIdy_MOY,w*p*dIdy);
              AJOUT(UPdIdz_MOY,u*p*dIdz);
              AJOUT(VPdIdz_MOY,v*p*dIdz);
              AJOUT(WPdIdz_MOY,w*p*dIdz);
              //^^^^^^^^^^^^^^^^^^^^^^^
              // ui*uj*dI/dxb : (18)
              AJOUT(UUdIdx_MOY,u*u*dIdx);
              AJOUT(VVdIdx_MOY,v*v*dIdx);
              AJOUT(WWdIdx_MOY,w*w*dIdx);
              AJOUT(UUdIdy_MOY,u*u*dIdy);
              AJOUT(VVdIdy_MOY,v*v*dIdy);
              AJOUT(WWdIdy_MOY,w*w*dIdy);
              AJOUT(UUdIdz_MOY,u*u*dIdz);
              AJOUT(VVdIdz_MOY,v*v*dIdz);
              AJOUT(WWdIdz_MOY,w*w*dIdz);
              AJOUT(UVdIdx_MOY,u*v*dIdx);
              AJOUT(UWdIdx_MOY,u*w*dIdx);
              AJOUT(VWdIdx_MOY,v*w*dIdx);
              AJOUT(UVdIdy_MOY,u*v*dIdy);
              AJOUT(UWdIdy_MOY,u*w*dIdy);
              AJOUT(VWdIdy_MOY,v*w*dIdy);
              AJOUT(UVdIdz_MOY,u*v*dIdz);
              AJOUT(UWdIdz_MOY,u*w*dIdz);
              AJOUT(VWdIdz_MOY,v*w*dIdz);
              //^^^^^^^^^^^^^^^^^^^^^^
              // p*dI/dxb : (3)
              AJOUT(PdIdx_MOY,p*dIdx);
              AJOUT(PdIdy_MOY,p*dIdy);
              AJOUT(PdIdz_MOY,p*dIdz);
              //^^^^^^^^^^^^^^^^^^^^^^
              // I*dI/dxb : (3)
              AJOUT(IdIdx_MOY,chi*dIdx);
              AJOUT(IdIdy_MOY,chi*dIdy);
              AJOUT(IdIdz_MOY,chi*dIdz);
              // 2. de la vitesse  :
              //^^^^^^^^^^^^^^^^^^^^^^^
              // I*dUidxb : (9)
              AJOUT(IdUdx_MOY,chi*dUdx);
              AJOUT(IdUdy_MOY,chi*dUdy);
              AJOUT(IdUdz_MOY,chi*dUdz);
              AJOUT(IdVdx_MOY,chi*dVdx);
              AJOUT(IdVdy_MOY,chi*dVdy);
              AJOUT(IdVdz_MOY,chi*dVdz);
              AJOUT(IdWdx_MOY,chi*dWdx);
              AJOUT(IdWdy_MOY,chi*dWdy);
              AJOUT(IdWdz_MOY,chi*dWdz);
              // I*P*dUidxb : (9)
              AJOUT(IPdUdx_MOY,chi*p*dUdx);
              AJOUT(IPdUdy_MOY,chi*p*dUdy);
              AJOUT(IPdUdz_MOY,chi*p*dUdz);
              AJOUT(IPdVdx_MOY,chi*p*dVdx);
              AJOUT(IPdVdy_MOY,chi*p*dVdy);
              AJOUT(IPdVdz_MOY,chi*p*dVdz);
              AJOUT(IPdWdx_MOY,chi*p*dWdx);
              AJOUT(IPdWdy_MOY,chi*p*dWdy);
              AJOUT(IPdWdz_MOY,chi*p*dWdz);
              // 3. de l'indicatrice et de la vitesse  :
              //^^^^^^^^^^^^^^^^^^^^^^^
              // dUidxb*dIdxb : (9)
              AJOUT(dUdxdIdx_MOY,dUdx*dIdx);
              AJOUT(dUdydIdy_MOY,dUdy*dIdy);
              AJOUT(dUdzdIdz_MOY,dUdz*dIdz);
              AJOUT(dVdxdIdx_MOY,dVdx*dIdx);
              AJOUT(dVdydIdy_MOY,dVdy*dIdy);
              AJOUT(dVdzdIdz_MOY,dVdz*dIdz);
              AJOUT(dWdxdIdx_MOY,dWdx*dIdx);
              AJOUT(dWdydIdy_MOY,dWdy*dIdy);
              AJOUT(dWdzdIdz_MOY,dWdz*dIdz);
              // Ui*dUjdxb*dIdxb : (27)
              AJOUT(UdUdxdIdx_MOY,u*dUdx*dIdx);
              AJOUT(UdVdxdIdx_MOY,u*dVdx*dIdx);
              AJOUT(UdWdxdIdx_MOY,u*dWdx*dIdx);
              AJOUT(UdUdydIdy_MOY,u*dUdy*dIdy);
              AJOUT(UdVdydIdy_MOY,u*dVdy*dIdy);
              AJOUT(UdWdydIdy_MOY,u*dWdy*dIdy);
              AJOUT(UdUdzdIdz_MOY,u*dUdz*dIdz);
              AJOUT(UdVdzdIdz_MOY,u*dVdz*dIdz);
              AJOUT(UdWdzdIdz_MOY,u*dWdz*dIdz);
              AJOUT(VdUdxdIdx_MOY,v*dUdx*dIdx);
              AJOUT(VdVdxdIdx_MOY,v*dVdx*dIdx);
              AJOUT(VdWdxdIdx_MOY,v*dWdx*dIdx);
              AJOUT(VdUdydIdy_MOY,v*dUdy*dIdy);
              AJOUT(VdVdydIdy_MOY,v*dVdy*dIdy);
              AJOUT(VdWdydIdy_MOY,v*dWdy*dIdy);
              AJOUT(VdUdzdIdz_MOY,v*dUdz*dIdz);
              AJOUT(VdVdzdIdz_MOY,v*dVdz*dIdz);
              AJOUT(VdWdzdIdz_MOY,v*dWdz*dIdz);
              AJOUT(WdUdxdIdx_MOY,w*dUdx*dIdx);
              AJOUT(WdVdxdIdx_MOY,w*dVdx*dIdx);
              AJOUT(WdWdxdIdx_MOY,w*dWdx*dIdx);
              AJOUT(WdUdydIdy_MOY,w*dUdy*dIdy);
              AJOUT(WdVdydIdy_MOY,w*dVdy*dIdy);
              AJOUT(WdWdydIdy_MOY,w*dWdy*dIdy);
              AJOUT(WdUdzdIdz_MOY,w*dUdz*dIdz);
              AJOUT(WdVdzdIdz_MOY,w*dVdz*dIdz);
              AJOUT(WdWdzdIdz_MOY,w*dWdz*dIdz);
              // 4. de la vitesse et de la vitesse (fois l'indic) :
              //^^^^^^^^^^^^^^^^^^^^^^^
              AJOUT(IdUdxdUdx_MOY,chi*dUdx*dUdx);
              AJOUT(IdUdxdUdy_MOY,chi*dUdx*dUdy);
              AJOUT(IdUdxdUdz_MOY,chi*dUdx*dUdz);
              AJOUT(IdUdxdVdx_MOY,chi*dUdx*dVdx);
              AJOUT(IdUdxdVdy_MOY,chi*dUdx*dVdy);
              AJOUT(IdUdxdVdz_MOY,chi*dUdx*dVdz);
              AJOUT(IdUdxdWdx_MOY,chi*dUdx*dWdx);
              AJOUT(IdUdxdWdy_MOY,chi*dUdx*dWdy);
              AJOUT(IdUdxdWdz_MOY,chi*dUdx*dWdz);
              AJOUT(IdUdydUdy_MOY,chi*dUdy*dUdy);
              AJOUT(IdUdydUdz_MOY,chi*dUdy*dUdz);
              AJOUT(IdUdydVdx_MOY,chi*dUdy*dVdx);
              AJOUT(IdUdydVdy_MOY,chi*dUdy*dVdy);
              AJOUT(IdUdydVdz_MOY,chi*dUdy*dVdz);
              AJOUT(IdUdydWdx_MOY,chi*dUdy*dWdx);
              AJOUT(IdUdydWdy_MOY,chi*dUdy*dWdy);
              AJOUT(IdUdydWdz_MOY,chi*dUdy*dWdz);
              AJOUT(IdUdzdUdz_MOY,chi*dUdz*dUdz);
              AJOUT(IdUdzdVdx_MOY,chi*dUdz*dVdx);
              AJOUT(IdUdzdVdy_MOY,chi*dUdz*dVdy);
              AJOUT(IdUdzdVdz_MOY,chi*dUdz*dVdz);
              AJOUT(IdUdzdWdx_MOY,chi*dUdz*dWdx);
              AJOUT(IdUdzdWdy_MOY,chi*dUdz*dWdy);
              AJOUT(IdUdzdWdz_MOY,chi*dUdz*dWdz);
              AJOUT(IdVdxdVdx_MOY,chi*dVdx*dVdx);
              AJOUT(IdVdxdVdy_MOY,chi*dVdx*dVdy);
              AJOUT(IdVdxdVdz_MOY,chi*dVdx*dVdz);
              AJOUT(IdVdxdWdx_MOY,chi*dVdx*dWdx);
              AJOUT(IdVdxdWdy_MOY,chi*dVdx*dWdy);
              AJOUT(IdVdxdWdz_MOY,chi*dVdx*dWdz);
              AJOUT(IdVdydVdy_MOY,chi*dVdy*dVdy);
              AJOUT(IdVdydVdz_MOY,chi*dVdy*dVdz);
              AJOUT(IdVdydWdx_MOY,chi*dVdy*dWdx);
              AJOUT(IdVdydWdy_MOY,chi*dVdy*dWdy);
              AJOUT(IdVdydWdz_MOY,chi*dVdy*dWdz);
              AJOUT(IdVdzdVdz_MOY,chi*dVdz*dVdz);
              AJOUT(IdVdzdWdx_MOY,chi*dVdz*dWdx);
              AJOUT(IdVdzdWdy_MOY,chi*dVdz*dWdy);
              AJOUT(IdVdzdWdz_MOY,chi*dVdz*dWdz);
              AJOUT(IdWdxdWdx_MOY,chi*dWdx*dWdx);
              AJOUT(IdWdxdWdy_MOY,chi*dWdx*dWdy);
              AJOUT(IdWdxdWdz_MOY,chi*dWdx*dWdz);
              AJOUT(IdWdydWdy_MOY,chi*dWdy*dWdy);
              AJOUT(IdWdydWdz_MOY,chi*dWdy*dWdz);
              AJOUT(IdWdzdWdz_MOY,chi*dWdz*dWdz);
              //^^^^^^^^^^^^^^^^^^^^^^^
              // Other correlations to avoid using gradI :
              AJOUT(IdPdx_MOY, chi*dPdx);
              AJOUT(IdPdy_MOY, chi*dPdy);
              AJOUT(IdPdz_MOY, chi*dPdz);
              AJOUT(IUdPdx_MOY, chi*UdPdx);
              AJOUT(IUdPdy_MOY, chi*u*dPdy);
              AJOUT(IUdPdz_MOY, chi*u*dPdz);
              AJOUT(IVdPdx_MOY, chi*v*dPdx);
              AJOUT(IVdPdy_MOY, chi*VdPdy);
              AJOUT(IVdPdz_MOY, chi*v*dPdz);
              AJOUT(IWdPdx_MOY, chi*w*dPdx);
              AJOUT(IWdPdy_MOY, chi*w*dPdy);
              AJOUT(IWdPdz_MOY, chi*WdPdz);
              //
              AJOUT(IddUdxx_MOY, chi*ddUdxdx);
              AJOUT(IddVdxx_MOY, chi*ddVdxdx);
              AJOUT(IddWdxx_MOY, chi*ddWdxdx);
              AJOUT(IddUdyy_MOY, chi*ddUdydy);
              AJOUT(IddVdyy_MOY, chi*ddVdydy);
              AJOUT(IddWdyy_MOY, chi*ddWdydy);
              AJOUT(IddUdzz_MOY, chi*ddUdzdz);
              AJOUT(IddVdzz_MOY, chi*ddVdzdz);
              AJOUT(IddWdzz_MOY, chi*ddWdzdz);
              //
              AJOUT(IUdUdx_MOY, chi*u*dUdx);
              AJOUT(IUdVdx_MOY, chi*u*dVdx);
              AJOUT(IUdWdx_MOY, chi*u*dWdx);
              AJOUT(IUdUdy_MOY, chi*u*dUdy);
              AJOUT(IUdVdy_MOY, chi*u*dVdy);
              AJOUT(IUdWdy_MOY, chi*u*dWdy);
              AJOUT(IUdUdz_MOY, chi*u*dUdz);
              AJOUT(IUdVdz_MOY, chi*u*dVdz);
              AJOUT(IUdWdz_MOY, chi*u*dWdz);
              //
              AJOUT(IVdUdx_MOY, chi*v*dUdx);
              AJOUT(IVdVdx_MOY, chi*v*dVdx);
              AJOUT(IVdWdx_MOY, chi*v*dWdx);
              AJOUT(IVdUdy_MOY, chi*v*dUdy);
              AJOUT(IVdVdy_MOY, chi*v*dVdy);
              AJOUT(IVdWdy_MOY, chi*v*dWdy);
              AJOUT(IVdUdz_MOY, chi*v*dUdz);
              AJOUT(IVdVdz_MOY, chi*v*dVdz);
              AJOUT(IVdWdz_MOY, chi*v*dWdz);
              //
              AJOUT(IWdUdx_MOY, chi*w*dUdx);
              AJOUT(IWdVdx_MOY, chi*w*dVdx);
              AJOUT(IWdWdx_MOY, chi*w*dWdx);
              AJOUT(IWdUdy_MOY, chi*w*dUdy);
              AJOUT(IWdVdy_MOY, chi*w*dVdy);
              AJOUT(IWdWdy_MOY, chi*w*dWdy);
              AJOUT(IWdUdz_MOY, chi*w*dUdz);
              AJOUT(IWdVdz_MOY, chi*w*dVdz);
              AJOUT(IWdWdz_MOY, chi*w*dWdz);
              //
              //^^^^^^^^^^^^^^^^^^^^^^^
              // Other quantities :
              //
              // Force : pot * dIdxb
              AJOUT(Fx_MOY,Fx);
              AJOUT(Fy_MOY,Fy);
              AJOUT(Fz_MOY,Fz);
              //
              // Repulsion : pot_rep * dIdxb
              AJOUT(Frx_MOY,Frx);
              AJOUT(Fry_MOY,Fry);
              AJOUT(Frz_MOY,Frz);
              //
              // Repulsion : pot_rep * dIdxb
              AJOUT(Frax_MOY,Frax);
              AJOUT(Fray_MOY,Fray);
              AJOUT(Fraz_MOY,Fraz);
              //
              // La dissipation :
              AJOUT(DISSIP_MOY,chi*pseudo_dissip);
              AJOUT(DISSIP_VAP_MOY,chiv*pseudo_dissip);
              AJOUT(TRUE_DISSIP_MOY,chi*true_dissip);
              AJOUT(TRUE_DISSIP_VAP_MOY,chiv*true_dissip);
              //
              // correlations 2 eme ordre vitesse (cote vapeur) :
              AJOUT(UUIv_MOY,u*u*chiv);
              AJOUT(VVIv_MOY,v*v*chiv);
              AJOUT(WWIv_MOY,w*w*chiv);
              AJOUT(UVIv_MOY,u*v*chiv);
              AJOUT(VWIv_MOY,v*w*chiv);
              AJOUT(UWIv_MOY,u*w*chiv);
              // b : bis (autre interpolation)
              AJOUT(UUIbv_MOY,uu*chiv);
              AJOUT(VVIbv_MOY,vv*chiv);
              AJOUT(WWIbv_MOY,ww*chiv);
              // correlations 2 eme ordre vitesse (liquide, autre interpolation)  :
              AJOUT(UUIb_MOY,uu*chi);
              AJOUT(VVIb_MOY,vv*chi);
              AJOUT(WWIb_MOY,ww*chi);
              // Il manque une seule correl triple qui aurait la meme interp :
              // u*v*w car il faut interpoler chaque compo avant le produit...
              AJOUT(UUUIb_MOY,uuu*chi);
              AJOUT(UUVIb_MOY,uu*v*chi);
              AJOUT(UUWIb_MOY,uu*w*chi);
              AJOUT(UVVIb_MOY,u*vv*chi);
              AJOUT(UWWIb_MOY,u*ww*chi);
              AJOUT(VVVIb_MOY,vvv*chi);
              AJOUT(VVWIb_MOY,vv*w*chi);
              AJOUT(VWWIb_MOY,v*ww*chi);
              AJOUT(WWWIb_MOY,www*chi);
              //
              // ai (et nk?):
              AJOUT(AI_MOY, ai);
              AJOUT(Nx_MOY, Nx);
              AJOUT(Ny_MOY, Ny);
              AJOUT(Nz_MOY, Nz);
              //
              // AJOUT(NK_MOY, ???);
              //
              // Les stats manquantes pour l'equation d'epsilon :
              // Attention, elles ne sont pas toutes implementees!!!
              AJOUT(IdUdxdUdxdUdx_MOY,chi * dUdx * dUdx * dUdx);
              AJOUT(IdUdydUdydUdx_MOY,chi * dUdy * dUdy * dUdx);
              AJOUT(IdUdzdUdzdUdx_MOY,chi * dUdz * dUdz * dUdx);
              AJOUT(IdUdxdVdxdUdy_MOY,chi * dUdx * dVdx * dUdy);
              AJOUT(IdUdydVdydUdy_MOY,chi * dUdy * dVdy * dUdy);
              AJOUT(IdUdzdVdzdUdy_MOY,chi * dUdz * dVdz * dUdy);
              AJOUT(IdUdxdWdxdUdz_MOY,chi * dUdx * dWdx * dUdz);
              AJOUT(IdUdydWdydUdz_MOY,chi * dUdy * dWdy * dUdz);
              AJOUT(IdUdzdWdzdUdz_MOY,chi * dUdz * dWdz * dUdz);
              AJOUT(IdVdxdUdxdVdx_MOY,chi * dVdx * dUdx * dVdx);
              AJOUT(IdVdydUdydVdx_MOY,chi * dVdy * dUdy * dVdx);
              AJOUT(IdVdzdUdzdVdx_MOY,chi * dVdz * dUdz * dVdx);
              AJOUT(IdVdxdVdxdVdy_MOY,chi * dVdx * dVdx * dVdy);
              AJOUT(IdVdydVdydVdy_MOY,chi * dVdy * dVdy * dVdy);
              AJOUT(IdVdzdVdzdVdy_MOY,chi * dVdz * dVdz * dVdy);
              AJOUT(IdVdxdWdxdVdz_MOY,chi * dVdx * dWdx * dVdz);
              AJOUT(IdVdydWdydVdz_MOY,chi * dVdy * dWdy * dVdz);
              AJOUT(IdVdzdWdzdVdz_MOY,chi * dVdz * dWdz * dVdz);
              AJOUT(IdWdxdUdxdWdx_MOY,chi * dWdx * dUdx * dWdx);
              AJOUT(IdWdydUdydWdx_MOY,chi * dWdy * dUdy * dWdx);
              AJOUT(IdWdzdUdzdWdx_MOY,chi * dWdz * dUdz * dWdx);
              AJOUT(IdWdxdVdxdWdy_MOY,chi * dWdx * dVdx * dWdy);
              AJOUT(IdWdydVdydWdy_MOY,chi * dWdy * dVdy * dWdy);
              AJOUT(IdWdzdVdzdWdy_MOY,chi * dWdz * dVdz * dWdy);
              AJOUT(IdWdxdWdxdWdz_MOY,chi * dWdx * dWdx * dWdz);
              AJOUT(IdWdydWdydWdz_MOY,chi * dWdy * dWdy * dWdz);
              AJOUT(IdWdzdWdzdWdz_MOY,chi * dWdz * dWdz * dWdz);
              AJOUT(IdUdxdUdxW_MOY,chi * dUdx * dUdx * w);
              AJOUT(IdUdydUdyW_MOY,chi * dUdy * dUdy * w);
              AJOUT(IdUdzdUdzW_MOY,chi * dUdz * dUdz * w);
              AJOUT(IdVdxdVdxW_MOY,chi * dVdx * dVdx * w);
              AJOUT(IdVdydVdyW_MOY,chi * dVdy * dVdy * w);
              AJOUT(IdVdzdVdzW_MOY,chi * dVdz * dVdz * w);
              AJOUT(IdWdxdWdxW_MOY,chi * dWdx * dWdx * w);
              AJOUT(IdWdydWdyW_MOY,chi * dWdy * dWdy * w);
              AJOUT(IdWdzdWdzW_MOY,chi * dWdz * dWdz * w);
              AJOUT(IdUdxddPdxdx_MOY,chi * dUdx * ddPdxdx);
              AJOUT(IdUdyddPdxdy_MOY,chi * dUdy * ddPdxdy);
              AJOUT(IdUdzddPdxdz_MOY,chi * dUdz * ddPdxdz);
              AJOUT(IdVdxddPdydx_MOY,chi * dVdx * ddPdydx);
              AJOUT(IdVdyddPdydy_MOY,chi * dVdy * ddPdydy);
              AJOUT(IdVdzddPdydz_MOY,chi * dVdz * ddPdydz);
              AJOUT(IdWdxddPdzdx_MOY,chi * dWdx * ddPdzdx);
              AJOUT(IdWdyddPdzdy_MOY,chi * dWdy * ddPdzdy);
              AJOUT(IdWdzddPdzdz_MOY,chi * dWdz * ddPdzdz);
              AJOUT(IddUdxdxddUdxdx_MOY,chi * ddUdxdx * ddUdxdx);
              AJOUT(IddUdxdyddUdxdy_MOY,chi * ddUdxdy * ddUdxdy);
              AJOUT(IddUdxdzddUdxdz_MOY,chi * ddUdxdz * ddUdxdz);
              AJOUT(IddUdydxddUdydx_MOY,chi * ddUdydx * ddUdydx);
              AJOUT(IddUdydyddUdydy_MOY,chi * ddUdydy * ddUdydy);
              AJOUT(IddUdydzddUdydz_MOY,chi * ddUdydz * ddUdydz);
              AJOUT(IddUdzdxddUdzdx_MOY,chi * ddUdzdx * ddUdzdx);
              AJOUT(IddUdzdyddUdzdy_MOY,chi * ddUdzdy * ddUdzdy);
              AJOUT(IddUdzdzddUdzdz_MOY,chi * ddUdzdz * ddUdzdz);
              AJOUT(IddVdxdxddVdxdx_MOY,chi * ddVdxdx * ddVdxdx);
              AJOUT(IddVdxdyddVdxdy_MOY,chi * ddVdxdy * ddVdxdy);
              AJOUT(IddVdxdzddVdxdz_MOY,chi * ddVdxdz * ddVdxdz);
              AJOUT(IddVdydxddVdydx_MOY,chi * ddVdydx * ddVdydx);
              AJOUT(IddVdydyddVdydy_MOY,chi * ddVdydy * ddVdydy);
              AJOUT(IddVdydzddVdydz_MOY,chi * ddVdydz * ddVdydz);
              AJOUT(IddVdzdxddVdzdx_MOY,chi * ddVdzdx * ddVdzdx);
              AJOUT(IddVdzdyddVdzdy_MOY,chi * ddVdzdy * ddVdzdy);
              AJOUT(IddVdzdzddVdzdz_MOY,chi * ddVdzdz * ddVdzdz);
              AJOUT(IddWdxdxddWdxdx_MOY,chi * ddWdxdx * ddWdxdx);
              AJOUT(IddWdxdyddWdxdy_MOY,chi * ddWdxdy * ddWdxdy);
              AJOUT(IddWdxdzddWdxdz_MOY,chi * ddWdxdz * ddWdxdz);
              AJOUT(IddWdydxddWdydx_MOY,chi * ddWdydx * ddWdydx);
              AJOUT(IddWdydyddWdydy_MOY,chi * ddWdydy * ddWdydy);
              AJOUT(IddWdydzddWdydz_MOY,chi * ddWdydz * ddWdydz);
              AJOUT(IddWdzdxddWdzdx_MOY,chi * ddWdzdx * ddWdzdx);
              AJOUT(IddWdzdyddWdzdy_MOY,chi * ddWdzdy * ddWdzdy);
              AJOUT(IddWdzdzddWdzdz_MOY,chi * ddWdzdz * ddWdzdz);
              AJOUT(dIdxddUdxdxdUdx_MOY,dIdx * ddUdxdx * dUdx);
              AJOUT(dIdxddUdxdydUdy_MOY,dIdx * ddUdxdy * dUdy);
              AJOUT(dIdxddUdxdzdUdz_MOY,dIdx * ddUdxdz * dUdz);
              AJOUT(dIdyddUdydxdUdx_MOY,dIdy * ddUdydx * dUdx);
              AJOUT(dIdyddUdydydUdy_MOY,dIdy * ddUdydy * dUdy);
              AJOUT(dIdyddUdydzdUdz_MOY,dIdy * ddUdydz * dUdz);
              AJOUT(dIdzddUdzdxdUdx_MOY,dIdz * ddUdzdx * dUdx);
              AJOUT(dIdzddUdzdydUdy_MOY,dIdz * ddUdzdy * dUdy);
              AJOUT(dIdzddUdzdzdUdz_MOY,dIdz * ddUdzdz * dUdz);
              AJOUT(dIdxddVdxdxdVdx_MOY,dIdx * ddVdxdx * dVdx);
              AJOUT(dIdxddVdxdydVdy_MOY,dIdx * ddVdxdy * dVdy);
              AJOUT(dIdxddVdxdzdVdz_MOY,dIdx * ddVdxdz * dVdz);
              AJOUT(dIdyddVdydxdVdx_MOY,dIdy * ddVdydx * dVdx);
              AJOUT(dIdyddVdydydVdy_MOY,dIdy * ddVdydy * dVdy);
              AJOUT(dIdyddVdydzdVdz_MOY,dIdy * ddVdydz * dVdz);
              AJOUT(dIdzddVdzdxdVdx_MOY,dIdz * ddVdzdx * dVdx);
              AJOUT(dIdzddVdzdydVdy_MOY,dIdz * ddVdzdy * dVdy);
              AJOUT(dIdzddVdzdzdVdz_MOY,dIdz * ddVdzdz * dVdz);
              AJOUT(dIdxddWdxdxdWdx_MOY,dIdx * ddWdxdx * dWdx);
              AJOUT(dIdxddWdxdydWdy_MOY,dIdx * ddWdxdy * dWdy);
              AJOUT(dIdxddWdxdzdWdz_MOY,dIdx * ddWdxdz * dWdz);
              AJOUT(dIdyddWdydxdWdx_MOY,dIdy * ddWdydx * dWdx);
              AJOUT(dIdyddWdydydWdy_MOY,dIdy * ddWdydy * dWdy);
              AJOUT(dIdyddWdydzdWdz_MOY,dIdy * ddWdydz * dWdz);
              AJOUT(dIdzddWdzdxdWdx_MOY,dIdz * ddWdzdx * dWdx);
              AJOUT(dIdzddWdzdydWdy_MOY,dIdz * ddWdzdy * dWdy);
              AJOUT(dIdzddWdzdzdWdz_MOY,dIdz * ddWdzdz * dWdz);
              AJOUT(dIdxddUdxdz_MOY,dIdx * ddUdxdz);
              AJOUT(dIdyddUdydz_MOY,dIdy * ddUdydz);
              AJOUT(dIdzddUdzdz_MOY,dIdz * ddUdzdz);
              AJOUT(dIdzdUdxdUdx_MOY,dIdz * dUdx * dUdx);
              AJOUT(dIdzdUdydUdy_MOY,dIdz * dUdy * dUdy);
              AJOUT(dIdzdUdzdUdz_MOY,dIdz * dUdz * dUdz);
              AJOUT(dIdzdVdxdVdx_MOY,dIdz * dVdx * dVdx);
              AJOUT(dIdzdVdydVdy_MOY,dIdz * dVdy * dVdy);
              AJOUT(dIdzdVdzdVdz_MOY,dIdz * dVdz * dVdz);
              AJOUT(dIdzdWdxdWdx_MOY,dIdz * dWdx * dWdx);
              AJOUT(dIdzdWdydWdy_MOY,dIdz * dWdy * dWdy);
              AJOUT(dIdzdWdzdWdz_MOY,dIdz * dWdz * dWdz);
              // Ajout de tout ce qui avait du dIdx qui est mainetenant aiNx...
              AJOUT(UaiNx_MOY,u*aiNx);
              AJOUT(VaiNx_MOY,v*aiNx);
              AJOUT(WaiNx_MOY,w*aiNx);
              AJOUT(UaiNy_MOY,u*aiNy);
              AJOUT(VaiNy_MOY,v*aiNy);
              AJOUT(WaiNy_MOY,w*aiNy);
              AJOUT(UaiNz_MOY,u*aiNz);
              AJOUT(VaiNz_MOY,v*aiNz);
              AJOUT(WaiNz_MOY,w*aiNz);
              AJOUT(UPaiNx_MOY,u*p*aiNx);
              AJOUT(VPaiNx_MOY,v*p*aiNx);
              AJOUT(WPaiNx_MOY,w*p*aiNx);
              AJOUT(UPaiNy_MOY,u*p*aiNy);
              AJOUT(VPaiNy_MOY,v*p*aiNy);
              AJOUT(WPaiNy_MOY,w*p*aiNy);
              AJOUT(UPaiNz_MOY,u*p*aiNz);
              AJOUT(VPaiNz_MOY,v*p*aiNz);
              AJOUT(WPaiNz_MOY,w*p*aiNz);
              AJOUT(UUaiNx_MOY,uu*aiNx);
              AJOUT(VVaiNx_MOY,vv*aiNx);
              AJOUT(WWaiNx_MOY,ww*aiNx);
              AJOUT(UUaiNy_MOY,uu*aiNy);
              AJOUT(VVaiNy_MOY,vv*aiNy);
              AJOUT(WWaiNy_MOY,ww*aiNy);
              AJOUT(UUaiNz_MOY,uu*aiNz);
              AJOUT(VVaiNz_MOY,vv*aiNz);
              AJOUT(WWaiNz_MOY,ww*aiNz);
              AJOUT(UVaiNx_MOY,u*v*aiNx);
              AJOUT(UWaiNx_MOY,u*w*aiNx);
              AJOUT(VWaiNx_MOY,v*w*aiNx);
              AJOUT(UVaiNy_MOY,u*v*aiNy);
              AJOUT(UWaiNy_MOY,u*w*aiNy);
              AJOUT(VWaiNy_MOY,v*w*aiNy);
              AJOUT(UVaiNz_MOY,u*v*aiNz);
              AJOUT(UWaiNz_MOY,u*w*aiNz);
              AJOUT(VWaiNz_MOY,v*w*aiNz);
              AJOUT(PaiNx_MOY,p*aiNx);
              AJOUT(PaiNy_MOY,p*aiNy);
              AJOUT(PaiNz_MOY,p*aiNz);
              AJOUT(IaiNx_MOY,chi*aiNx);
              AJOUT(IaiNy_MOY,chi*aiNy);
              AJOUT(IaiNz_MOY,chi*aiNz);
              AJOUT(dUdxaiNx_MOY,dUdx*aiNx);
              AJOUT(dVdxaiNx_MOY,dVdx*aiNx);
              AJOUT(dWdxaiNx_MOY,dWdx*aiNx);
              AJOUT(dUdyaiNy_MOY,dUdy*aiNy);
              AJOUT(dVdyaiNy_MOY,dVdy*aiNy);
              AJOUT(dWdyaiNy_MOY,dWdy*aiNy);
              AJOUT(dUdzaiNz_MOY,dUdz*aiNz);
              AJOUT(dVdzaiNz_MOY,dVdz*aiNz);
              AJOUT(dWdzaiNz_MOY,dWdz*aiNz);
              AJOUT(UdUdxaiNx_MOY,u*dUdx*aiNx);
              AJOUT(UdVdxaiNx_MOY,u*dVdx*aiNx);
              AJOUT(UdWdxaiNx_MOY,u*dWdx*aiNx);
              AJOUT(UdUdyaiNy_MOY,u*dUdy*aiNy);
              AJOUT(UdVdyaiNy_MOY,u*dVdy*aiNy);
              AJOUT(UdWdyaiNy_MOY,u*dWdy*aiNy);
              AJOUT(UdUdzaiNz_MOY,u*dUdz*aiNz);
              AJOUT(UdVdzaiNz_MOY,u*dVdz*aiNz);
              AJOUT(UdWdzaiNz_MOY,u*dWdz*aiNz);
              AJOUT(VdUdxaiNx_MOY,v*dUdx*aiNx);
              AJOUT(VdVdxaiNx_MOY,v*dVdx*aiNx);
              AJOUT(VdWdxaiNx_MOY,v*dWdx*aiNx);
              AJOUT(VdUdyaiNy_MOY,v*dUdy*aiNy);
              AJOUT(VdVdyaiNy_MOY,v*dVdy*aiNy);
              AJOUT(VdWdyaiNy_MOY,v*dWdy*aiNy);
              AJOUT(VdUdzaiNz_MOY,v*dUdz*aiNz);
              AJOUT(VdVdzaiNz_MOY,v*dVdz*aiNz);
              AJOUT(VdWdzaiNz_MOY,v*dWdz*aiNz);
              AJOUT(WdUdxaiNx_MOY,w*dUdx*aiNx);
              AJOUT(WdVdxaiNx_MOY,w*dVdx*aiNx);
              AJOUT(WdWdxaiNx_MOY,w*dWdx*aiNx);
              AJOUT(WdUdyaiNy_MOY,w*dUdy*aiNy);
              AJOUT(WdVdyaiNy_MOY,w*dVdy*aiNy);
              AJOUT(WdWdyaiNy_MOY,w*dWdy*aiNy);
              AJOUT(WdUdzaiNz_MOY,w*dUdz*aiNz);
              AJOUT(WdVdzaiNz_MOY,w*dVdz*aiNz);
              AJOUT(WdWdzaiNz_MOY,w*dWdz*aiNz);
              AJOUT(aiNxddUdxdxdUdx_MOY,aiNx * ddUdxdx * dUdx);
              AJOUT(aiNxddUdxdydUdy_MOY,aiNx * ddUdxdy * dUdy);
              AJOUT(aiNxddUdxdzdUdz_MOY,aiNx * ddUdxdz * dUdz);
              AJOUT(aiNyddUdydydUdy_MOY,aiNy * ddUdydy * dUdy);
              AJOUT(aiNyddUdydzdUdz_MOY,aiNy * ddUdydz * dUdz);
              AJOUT(aiNzddUdzdzdUdz_MOY,aiNz * ddUdzdz * dUdz);
              AJOUT(aiNxddVdxdxdVdx_MOY,aiNx * ddVdxdx * dVdx);
              AJOUT(aiNxddVdxdydVdy_MOY,aiNx * ddVdxdy * dVdy);
              AJOUT(aiNxddVdxdzdVdz_MOY,aiNx * ddVdxdz * dVdz);
              AJOUT(aiNyddVdydydVdy_MOY,aiNy * ddVdydy * dVdy);
              AJOUT(aiNyddVdydzdVdz_MOY,aiNy * ddVdydz * dVdz);
              AJOUT(aiNzddVdzdzdVdz_MOY,aiNz * ddVdzdz * dVdz);
              AJOUT(aiNxddWdxdxdWdx_MOY,aiNx * ddWdxdx * dWdx);
              AJOUT(aiNxddWdxdydWdy_MOY,aiNx * ddWdxdy * dWdy);
              AJOUT(aiNxddWdxdzdWdz_MOY,aiNx * ddWdxdz * dWdz);
              AJOUT(aiNyddWdydydWdy_MOY,aiNy * ddWdydy * dWdy);
              AJOUT(aiNyddWdydzdWdz_MOY,aiNy * ddWdydz * dWdz);
              AJOUT(aiNzddWdzdzdWdz_MOY,aiNz * ddWdzdz * dWdz);
              AJOUT(aiNxddUdxdz_MOY,aiNx * ddUdxdz);
              AJOUT(aiNyddUdydz_MOY,aiNy * ddUdydz);
              AJOUT(aiNzddUdzdz_MOY,aiNz * ddUdzdz);
              AJOUT(aiNzdUdxdUdx_MOY,aiNz * dUdx * dUdx);
              AJOUT(aiNzdUdydUdy_MOY,aiNz * dUdy * dUdy);
              AJOUT(aiNzdUdzdUdz_MOY,aiNz * dUdz * dUdz);
              AJOUT(aiNzdVdxdVdx_MOY,aiNz * dVdx * dVdx);
              AJOUT(aiNzdVdydVdy_MOY,aiNz * dVdy * dVdy);
              AJOUT(aiNzdVdzdVdz_MOY,aiNz * dVdz * dVdz);
              AJOUT(aiNzdWdxdWdx_MOY,aiNz * dWdx * dWdx);
              AJOUT(aiNzdWdydWdy_MOY,aiNz * dWdy * dWdy);
              AJOUT(aiNzdWdzdWdz_MOY,aiNz * dWdz * dWdz);

              AJOUT(aiNx_MOY,aiNx);
              AJOUT(aiNy_MOY,aiNy);
              AJOUT(aiNz_MOY,aiNz);
              AJOUT(kaiNx_MOY,kai*Nx);
              AJOUT(kaiNy_MOY,kai*Ny);
              AJOUT(kaiNz_MOY,kai*Nz);

              AJOUT(aiNxx_MOY,ai*Nxx);
              AJOUT(aiNxy_MOY,ai*Nxy);
              AJOUT(aiNxz_MOY,ai*Nxz);
              AJOUT(aiNyy_MOY,ai*Nyy);
              AJOUT(aiNyz_MOY,ai*Nyz);
              AJOUT(aiNzz_MOY,ai*Nzz);

              AJOUT(aiNNxxxx_MOY,ai*Nxx*Nxx);
              AJOUT(aiNNxxxy_MOY,ai*Nxx*Nxy);
              AJOUT(aiNNxxxz_MOY,ai*Nxx*Nxz);
              AJOUT(aiNNxxyy_MOY,ai*Nxx*Nyy);
              AJOUT(aiNNxxyz_MOY,ai*Nxx*Nyz);
              AJOUT(aiNNxxzz_MOY,ai*Nxx*Nzz);

              AJOUT(aiNNxyyy_MOY,ai*Nxy*Nyy);
              AJOUT(aiNNxyyz_MOY,ai*Nxy*Nyz);
              AJOUT(aiNNxyzz_MOY,ai*Nxy*Nzz);
              AJOUT(aiNNxzyy_MOY,ai*Nxz*Nyy);
              AJOUT(aiNNxzyz_MOY,ai*Nxz*Nyz);
              AJOUT(aiNNxzzz_MOY,ai*Nxz*Nzz);

              AJOUT(aiNNyyyy_MOY,ai*Nyy*Nyy);
              AJOUT(aiNNyyyz_MOY,ai*Nyy*Nyz);
              AJOUT(aiNNyyzz_MOY,ai*Nyy*Nzz);
              AJOUT(aiNNyzzz_MOY,ai*Nyz*Nzz);
              AJOUT(aiNNzzzz_MOY,ai*Nzz*Nzz);

              AJOUT(kai_MOY,kai);





              if(coef_immobilisation_>1e-16)
                {
                  double intu = integral_vitesse_i(i,j,k);
                  double intv = integral_vitesse_j(i,j,k);
                  double intw = integral_vitesse_k(i,j,k);
                  double intp = integral_pression(i,j,k);
                  double rappelx=force_rappel[0](i,j,k);
                  double rappely=force_rappel[1](i,j,k);
                  double rappelz=force_rappel[2](i,j,k);
                  double p_seuil =0., p_seuil_isup=0., p_seuil_iinf=0., p_seuil_jsup=0., p_seuil_jinf=0. , p_seuil_ksup=0. , p_seuil_kinf=0.;
                  p_seuil = Calculer_valeur_seuil(p_seuil, pression(i,j,k), p_seuil_max, p_seuil_min);
                  p_seuil_isup = Calculer_valeur_seuil(p_seuil_isup, pression(i+1,j,k), p_seuil_max, p_seuil_min);
                  p_seuil_jsup = Calculer_valeur_seuil(p_seuil_jsup, pression(i,j+1,k), p_seuil_max, p_seuil_min);
                  p_seuil_ksup = Calculer_valeur_seuil(p_seuil_ksup, pression(i,j,k+1), p_seuil_max, p_seuil_min);
                  p_seuil_iinf = Calculer_valeur_seuil(p_seuil_iinf, pression(i-1,j,k), p_seuil_max, p_seuil_min);
                  p_seuil_jinf = Calculer_valeur_seuil(p_seuil_jinf, pression(i,j-1,k), p_seuil_max, p_seuil_min);
                  p_seuil_kinf = Calculer_valeur_seuil(p_seuil_kinf, pression(i,j,k-1), p_seuil_max, p_seuil_min);
                  double dPdx_seuil;
                  double dPdy_seuil;
                  double dPdz_seuil;
                  dPdx_seuil=((p_seuil_isup-p_seuil)*(1./dx_)+(p_seuil-p_seuil_iinf)*(1./dx_))*0.5;
                  dPdy_seuil=((p_seuil_jsup-p_seuil)*(1./dy_)+(p_seuil-p_seuil_jinf)*(1./dy_))*0.5;
                  dPdz_seuil=((p_seuil_ksup-p_seuil)*(1./dz)+(p_seuil-p_seuil_kinf)*(1./dz))*0.5;
                  chi = floor(indicatrice(i,j,k));
                  chiv = floor(1.-indicatrice(i,j,k));
                  double chinp = 1.0-indicatrice_np(i,j,k);

                  double  dUintdx=0., dVintdx=0., dWintdx=0., dUintdy=0., dVintdy=0., dWintdy=0., dUintdz=0., dVintdz=0., dWintdz=0.;
                  double  ddUintdxdx=0., ddVintdxdx=0., ddWintdxdx=0., ddUintdydy=0., ddVintdydy=0., ddWintdydy=0., ddUintdzdz=0., ddVintdzdz=0., ddWintdzdz=0.;
                  double pseudo_dissipint = calculer_gradients_vitesse(integral_vitesse_i, integral_vitesse_j, integral_vitesse_k,
                                                                       i,j,k,
                                                                       dz,
                                                                       dUintdx, dVintdx, dWintdx,
                                                                       dUintdy, dVintdy, dWintdy,
                                                                       dUintdz, dVintdz, dWintdz,
                                                                       ddUintdxdx, ddVintdxdx, ddWintdxdx,
                                                                       ddUintdydy, ddVintdydy, ddWintdydy,
                                                                       ddUintdzdz, ddVintdzdz, ddWintdzdz,
                                                                       on_the_first_cell, on_the_last_cell);
//                  double true_dissipint = calculer_vraie_dissipation(pseudo_dissipint,
//                                                                     dUintdx, dUintdy, dUintdz,
//                                                                     dVintdx, dVintdy, dVintdz,
//                                                                     dWintdx, dWintdy, dWintdz);

                  // ajout des stats pour bulle fixe et moyenne tmp
                  AJOUT(UI_INT,intu*chi);
                  AJOUT(VI_INT,intv*chi);
                  AJOUT(WI_INT,intw*chi);
                  AJOUT(UI_INTUI_INT, intu*intu*chi);
                  AJOUT(UI_INTVI_INT, intu*intv*chi);
                  AJOUT(UI_INTWI_INT, intu*intw*chi);
                  AJOUT(VI_INTVI_INT, intv*intv*chi );
                  AJOUT(VI_INTWI_INT, intv*intw*chi);
                  AJOUT(WI_INTWI_INT, intw*intw*chi);
                  AJOUT(UI_INTUI_INTUI_INT, intu*intu*intu*chi);
                  AJOUT(UI_INTUI_INTVI_INT, intu*intu*intv*chi );
                  AJOUT(UI_INTUI_INTWI_INT, intu*intu*intw*chi );
                  AJOUT(UI_INTVI_INTUI_INT, intu*intv*intu*chi );
                  AJOUT(UI_INTVI_INTVI_INT, intu*intv*intv*chi );
                  AJOUT(UI_INTVI_INTWI_INT, intu*intv*intw*chi );
                  AJOUT(UI_INTWI_INTUI_INT, intu*intw*intu*chi );
                  AJOUT(UI_INTWI_INTVI_INT, intu*intw*intv*chi );
                  AJOUT(UI_INTWI_INTWI_INT, intu*intw*intw*chi );
                  AJOUT(VI_INTUI_INTUI_INT, intv*intu*intu*chi);
                  AJOUT(VI_INTUI_INTVI_INT, intv*intu*intv*chi );
                  AJOUT(VI_INTUI_INTWI_INT, intv*intu*intw*chi );
                  AJOUT(VI_INTVI_INTUI_INT, intv*intv*intu*chi );
                  AJOUT(VI_INTVI_INTVI_INT, intv*intv*intv*chi );
                  AJOUT(VI_INTVI_INTWI_INT, intv*intv*intw*chi );
                  AJOUT(VI_INTWI_INTUI_INT, intv*intw*intu*chi );
                  AJOUT(VI_INTWI_INTVI_INT, intv*intw*intv*chi );
                  AJOUT(VI_INTWI_INTWI_INT, intv*intw*intw*chi );
                  AJOUT(WI_INTUI_INTUI_INT, intw*intu*intu*chi);
                  AJOUT(WI_INTUI_INTVI_INT, intw*intu*intv*chi );
                  AJOUT(WI_INTUI_INTWI_INT, intw*intu*intw*chi );
                  AJOUT(WI_INTVI_INTUI_INT, intw*intv*intu*chi );
                  AJOUT(WI_INTVI_INTVI_INT, intw*intv*intv*chi );
                  AJOUT(WI_INTVI_INTWI_INT, intw*intv*intw*chi );
                  AJOUT(WI_INTWI_INTUI_INT, intw*intw*intu*chi );
                  AJOUT(WI_INTWI_INTVI_INT, intw*intw*intv*chi );
                  AJOUT(WI_INTWI_INTWI_INT, intw*intw*intw*chi );
                  AJOUT( UI_INTUIUI, intu*u*u*chi );
                  AJOUT( UI_INTUIVI, intu*u*v*chi );
                  AJOUT( UI_INTUIWI, intu*u*w*chi );
                  AJOUT( UI_INTVIUI, intu*v*u*chi );
                  AJOUT( UI_INTVIVI, intu*v*v*chi );
                  AJOUT( UI_INTVIWI, intu*v*w*chi );
                  AJOUT( UI_INTWIUI, intu*w*u*chi );
                  AJOUT( UI_INTWIVI, intu*w*v*chi );
                  AJOUT( UI_INTWIWI, intu*w*w*chi );
                  AJOUT( VI_INTUIUI, intv*u*u*chi );
                  AJOUT( VI_INTUIVI, intv*u*v*chi );
                  AJOUT( VI_INTUIWI, intv*u*w*chi );
                  AJOUT( VI_INTVIUI, intv*v*u*chi );
                  AJOUT( VI_INTVIVI, intv*v*v*chi );
                  AJOUT( VI_INTVIWI, intv*v*w*chi );
                  AJOUT( VI_INTWIUI, intv*w*u*chi );
                  AJOUT( VI_INTWIVI, intv*w*v*chi );
                  AJOUT( VI_INTWIWI, intv*w*w*chi );
                  AJOUT( WI_INTUIUI, intw*u*u*chi );
                  AJOUT( WI_INTUIVI, intw*u*v*chi );
                  AJOUT( WI_INTUIWI, intw*u*w*chi );
                  AJOUT( WI_INTVIUI, intw*v*u*chi );
                  AJOUT( WI_INTVIVI, intw*v*v*chi );
                  AJOUT( WI_INTVIWI, intw*v*w*chi );
                  AJOUT( WI_INTWIUI, intw*w*u*chi );
                  AJOUT( WI_INTWIVI, intw*w*v*chi );
                  AJOUT( WI_INTWIWI, intw*w*w*chi );
                  AJOUT(UI_INTP_INT, intu*intp*chi );
                  AJOUT(VI_INTP_INT, intv*intp*chi  );
                  AJOUT(WI_INTP_INT, intw*intp*chi  );

                  AJOUT(dUINTdxdUINTdx, dUintdx*dUintdx*chi );
                  AJOUT(dUINTdxdUINTdy, dUintdx*dUintdy*chi );
                  AJOUT(dUINTdxdUINTdz, dUintdx*dUintdz*chi );
                  AJOUT(dUINTdydUINTdy, dUintdy*dUintdy*chi );
                  AJOUT(dUINTdydUINTdz, dUintdy*dUintdz*chi );
                  AJOUT(dUINTdzdUINTdz, dUintdz*dUintdz*chi );
                  AJOUT(dUINTdxdVINTdx, dUintdx*dVintdx*chi );
                  AJOUT(dUINTdxdVINTdy, dUintdx*dVintdy*chi );
                  AJOUT(dUINTdxdVINTdz, dUintdx*dVintdz*chi );
                  AJOUT(dUINTdydVINTdy, dUintdy*dVintdy*chi );
                  AJOUT(dUINTdydVINTdz, dUintdy*dVintdz*chi );
                  AJOUT(dUINTdzdVINTdz, dUintdz*dVintdz*chi );
                  AJOUT(dUINTdxdWINTdx, dUintdx*dWintdx*chi );
                  AJOUT(dUINTdxdWINTdy, dUintdx*dWintdy*chi );
                  AJOUT(dUINTdxdWINTdz, dUintdx*dWintdz*chi );
                  AJOUT(dUINTdydWINTdy, dUintdy*dWintdy*chi );
                  AJOUT(dUINTdydWINTdz, dUintdy*dWintdz*chi );
                  AJOUT(dUINTdzdWINTdz, dUintdz*dWintdz*chi );


                  AJOUT(dVINTdxdUINTdx, dVintdx*dUintdx*chi );
                  AJOUT(dVINTdxdUINTdy, dVintdx*dUintdy*chi );
                  AJOUT(dVINTdxdUINTdz, dVintdx*dUintdz*chi );
                  AJOUT(dVINTdydUINTdy, dVintdy*dUintdy*chi );
                  AJOUT(dVINTdydUINTdz, dVintdy*dUintdz*chi );
                  AJOUT(dVINTdzdUINTdz, dVintdz*dUintdz*chi );
                  AJOUT(dVINTdxdVINTdx, dVintdx*dVintdx*chi );
                  AJOUT(dVINTdxdVINTdy, dVintdx*dVintdy*chi );
                  AJOUT(dVINTdxdVINTdz, dVintdx*dVintdz*chi );
                  AJOUT(dVINTdydVINTdy, dVintdy*dVintdy*chi );
                  AJOUT(dVINTdydVINTdz, dVintdy*dVintdz*chi );
                  AJOUT(dVINTdzdVINTdz, dVintdz*dVintdz*chi );
                  AJOUT(dVINTdxdWINTdx, dVintdx*dWintdx*chi );
                  AJOUT(dVINTdxdWINTdy, dVintdx*dWintdy*chi );
                  AJOUT(dVINTdxdWINTdz, dVintdx*dWintdz*chi );
                  AJOUT(dVINTdydWINTdy, dVintdy*dWintdy*chi );
                  AJOUT(dVINTdydWINTdz, dVintdy*dWintdz*chi );
                  AJOUT(dVINTdzdWINTdz, dVintdz*dWintdz*chi );

                  AJOUT(dWINTdxdUINTdx, dWintdx*dUintdx*chi );
                  AJOUT(dWINTdxdUINTdy, dWintdx*dUintdy*chi );
                  AJOUT(dWINTdxdUINTdz, dWintdx*dUintdz*chi );
                  AJOUT(dWINTdydUINTdy, dWintdy*dUintdy*chi );
                  AJOUT(dWINTdydUINTdz, dWintdy*dUintdz*chi );
                  AJOUT(dWINTdzdUINTdz, dWintdz*dUintdz*chi );
                  AJOUT(dWINTdxdVINTdx, dWintdx*dVintdx*chi );
                  AJOUT(dWINTdxdVINTdy, dWintdx*dVintdy*chi );
                  AJOUT(dWINTdxdVINTdz, dWintdx*dVintdz*chi );
                  AJOUT(dWINTdydVINTdy, dWintdy*dVintdy*chi );
                  AJOUT(dWINTdydVINTdz, dWintdy*dVintdz*chi );
                  AJOUT(dWINTdzdVINTdz, dWintdz*dVintdz*chi );
                  AJOUT(dWINTdxdWINTdx, dWintdx*dWintdx*chi );
                  AJOUT(dWINTdxdWINTdy, dWintdx*dWintdy*chi );
                  AJOUT(dWINTdxdWINTdz, dWintdx*dWintdz*chi );
                  AJOUT(dWINTdydWINTdy, dWintdy*dWintdy*chi );
                  AJOUT(dWINTdydWINTdz, dWintdy*dWintdz*chi );
                  AJOUT(dWINTdzdWINTdz, dWintdz*dWintdz*chi );
                  AJOUT(P_INTdUINTdx, intp*dUintdx*chi);
                  AJOUT(P_INTdUINTdy, intp*dUintdy*chi);
                  AJOUT(P_INTdUINTdz, intp*dUintdz*chi);
                  AJOUT(P_INTdVINTdx, intp*dVintdx*chi);
                  AJOUT(P_INTdVINTdy, intp*dVintdy*chi);
                  AJOUT(P_INTdVINTdz, intp*dVintdz*chi);
                  AJOUT(P_INTdWINTdx, intp*dWintdx*chi);
                  AJOUT(P_INTdWINTdy, intp*dWintdy*chi);
                  AJOUT(P_INTdWINTdz, intp*dWintdz*chi);
                  AJOUT(DISSIP_INT,chi*pseudo_dissipint);
//                  AJOUT(DISSIP_VAP_INT,chiv*pseudo_dissipint);
//                  AJOUT(TRUE_DISSIP_INT,chi*true_dissipint);
//                  AJOUT(TRUE_DISSIP_VAP_INT,chiv*true_dissipint);
                  AJOUT(U_NOPERTURBE,u*chinp);
                  AJOUT(V_NOPERTURBE,v*chinp);
                  AJOUT(W_NOPERTURBE,w*chinp);
                  AJOUT(P_NOPERTURBE,p*chinp);
                  AJOUT(I_NP,chinp);
                  AJOUT(P_NP, p*chinp);
                  AJOUT(PIU_NP, u*p*chinp);
                  AJOUT(PIV_NP, v*p*chinp);
                  AJOUT(PIW_NP, w*p*chinp);
                  AJOUT(IPdUdx_NP,chinp*dUdx*p );
                  AJOUT(IPdVdx_NP,chinp*dVdx*p );
                  AJOUT(IPdWdx_NP,chinp*dWdx*p );
                  AJOUT(IPdUdy_NP,chinp*dUdy*p );
                  AJOUT(IPdVdy_NP,chinp*dVdy*p );
                  AJOUT(IPdWdy_NP,chinp*dWdy*p );
                  AJOUT(IPdUdz_NP,chinp*dUdz*p );
                  AJOUT(IPdVdz_NP,chinp*dVdz*p );
                  AJOUT(IPdWdz_NP,chinp*dWdz*p );
                  AJOUT(IdPdx_NP, chinp*dPdx );
                  AJOUT(IdPdy_NP, chinp*dPdy);
                  AJOUT(IdPdz_NP, chinp*dPdz);
                  AJOUT(PI_seuil, p_seuil*chi);
                  AJOUT(PIU_seuil,u*p_seuil*chi );
                  AJOUT(PIV_seuil,v*p_seuil*chi );
                  AJOUT(PIW_seuil,w*p_seuil*chi );
                  AJOUT(IPdUdx_seuil, chi*dUdx*p_seuil);
                  AJOUT(IPdVdx_seuil, chi*dVdx*p_seuil);
                  AJOUT(IPdWdx_seuil, chi*dWdx*p_seuil);
                  AJOUT(IPdUdy_seuil, chi*dUdy*p_seuil);
                  AJOUT(IPdVdy_seuil, chi*dVdy*p_seuil);
                  AJOUT(IPdWdy_seuil, chi*dWdy*p_seuil);
                  AJOUT(IPdUdz_seuil, chi*dUdz*p_seuil);
                  AJOUT(IPdVdz_seuil, chi*dVdz*p_seuil);
                  AJOUT(IPdWdz_seuil, chi*dWdz*p_seuil);
                  AJOUT(IdPdx_seuil, dPdx_seuil*chi );
                  AJOUT(IdPdy_seuil, dPdy_seuil*chi );
                  AJOUT(IdPdz_seuil, dPdz_seuil*chi );
                  AJOUT(P_NP_seuil, p_seuil*chinp);
                  AJOUT(PIU_NP_seuil,u*p_seuil*chinp );
                  AJOUT(PIV_NP_seuil,v*p_seuil*chinp );
                  AJOUT(PIW_NP_seuil,w*p_seuil*chinp );
                  AJOUT(IPdUdx_NP_seuil, chinp*dUdx*p_seuil);
                  AJOUT(IPdVdx_NP_seuil, chinp*dVdx*p_seuil);
                  AJOUT(IPdWdx_NP_seuil, chinp*dWdx*p_seuil);
                  AJOUT(IPdUdy_NP_seuil, chinp*dUdy*p_seuil);
                  AJOUT(IPdVdy_NP_seuil, chinp*dVdy*p_seuil);
                  AJOUT(IPdWdy_NP_seuil, chinp*dWdy*p_seuil);
                  AJOUT(IPdUdz_NP_seuil, chinp*dUdz*p_seuil);
                  AJOUT(IPdVdz_NP_seuil, chinp*dVdz*p_seuil);
                  AJOUT(IPdWdz_NP_seuil, chinp*dWdz*p_seuil);
                  AJOUT(IdPdx_NP_seuil, dPdx_seuil*chinp );
                  AJOUT(IdPdy_NP_seuil, dPdy_seuil*chinp );
                  AJOUT(IdPdz_NP_seuil, dPdz_seuil*chinp );
                  AJOUT(Ivrappelx, chiv*rappelx);
                  AJOUT(Ivrappely, chiv*rappely);
                  AJOUT(Ivrappelz, chiv*rappelz);
                  AJOUT(IP_INT, chi*intp );
                  AJOUT(UIvrappelx, chiv*rappelx*u);
                  AJOUT(VIvrappelx, chiv*rappelx*v);
                  AJOUT(WIvrappelx, chiv*rappelx*w);
                  AJOUT(UIvrappely, chiv*rappely*u);
                  AJOUT(VIvrappely, chiv*rappely*v);
                  AJOUT(WIvrappely, chiv*rappely*w);
                  AJOUT(UIvrappelz, chiv*rappelz*u);
                  AJOUT(VIvrappelz, chiv*rappelz*v);
                  AJOUT(WIvrappelz, chiv*rappelz*w);
                }

              //NEW TERMS FOR THE DEVIATORIC PART OF THE VISCOUS STRESS TENSOR
              AJOUT(dUdyaiNx_MOY, dUdy*aiNx);
              AJOUT(dUdzaiNx_MOY, dUdz*aiNx);
              AJOUT(dVdxaiNy_MOY, dVdx*aiNy);
              AJOUT(dVdzaiNy_MOY, dVdz*aiNy);
              AJOUT(dWdxaiNz_MOY, dWdx*aiNz);
              AJOUT(dWdyaiNz_MOY, dWdy*aiNz);

              // Iv*dUidxb : (9)
              AJOUT(dUdxIv_MOY,chiv*dUdx);
              AJOUT(dUdyIv_MOY,chiv*dUdy);
              AJOUT(dUdzIv_MOY,chiv*dUdz);
              AJOUT(dVdxIv_MOY,chiv*dVdx);
              AJOUT(dVdyIv_MOY,chiv*dVdy);
              AJOUT(dVdzIv_MOY,chiv*dVdz);
              AJOUT(dWdxIv_MOY,chiv*dWdx);
              AJOUT(dWdyIv_MOY,chiv*dWdy);
              AJOUT(dWdzIv_MOY,chiv*dWdz);
              // Attempt for viscous diffusion
              AJOUT(IddUdxdxU_MOY, chi*ddUdxdx*u);
              AJOUT(IddVdxdxU_MOY, chi*ddVdxdx*u);
              AJOUT(IddWdxdxU_MOY, chi*ddWdxdx*u);
              AJOUT(IddUdxdxV_MOY, chi*ddUdxdx*v);
              AJOUT(IddVdxdxV_MOY, chi*ddVdxdx*v);
              AJOUT(IddWdxdxV_MOY, chi*ddWdxdx*v);
              AJOUT(IddUdxdxW_MOY, chi*ddUdxdx*w);
              AJOUT(IddVdxdxW_MOY, chi*ddVdxdx*w);
              AJOUT(IddWdxdxW_MOY, chi*ddWdxdx*w);

              AJOUT(IddUdydyU_MOY, chi*ddUdydy*u);
              AJOUT(IddVdydyU_MOY, chi*ddVdydy*u);
              AJOUT(IddWdydyU_MOY, chi*ddWdydy*u);
              AJOUT(IddUdydyV_MOY, chi*ddUdydy*v);
              AJOUT(IddVdydyV_MOY, chi*ddVdydy*v);
              AJOUT(IddWdydyV_MOY, chi*ddWdydy*v);
              AJOUT(IddUdydyW_MOY, chi*ddUdydy*w);
              AJOUT(IddVdydyW_MOY, chi*ddVdydy*w);
              AJOUT(IddWdydyW_MOY, chi*ddWdydy*w);

              AJOUT(IddUdzdzU_MOY, chi*ddUdzdz*u);
              AJOUT(IddVdzdzU_MOY, chi*ddVdzdz*u);
              AJOUT(IddWdzdzU_MOY, chi*ddWdzdz*u);
              AJOUT(IddUdzdzV_MOY, chi*ddUdzdz*v);
              AJOUT(IddVdzdzV_MOY, chi*ddVdzdz*v);
              AJOUT(IddWdzdzV_MOY, chi*ddWdzdz*v);
              AJOUT(IddUdzdzW_MOY, chi*ddUdzdz*w);
              AJOUT(IddVdzdzW_MOY, chi*ddVdzdz*w);
              AJOUT(IddWdzdzW_MOY, chi*ddWdzdz*w);

              AJOUT(Ix_MOY , x*chi);
              AJOUT(xaiNx_MOY, x*aiNx);
              AJOUT(xaiNy_MOY, x*aiNy);
              AJOUT(xaiNz_MOY, x*aiNz);

              AJOUT(dUdxIx_MOY, x*chi*dUdx);
              AJOUT(dUdyIx_MOY, x*chi*dUdy);
              AJOUT(dUdzIx_MOY, x*chi*dUdz);
              AJOUT(dVdxIx_MOY, x*chi*dVdx);
              AJOUT(dVdyIx_MOY, x*chi*dVdy);
              AJOUT(dVdzIx_MOY, x*chi*dVdz);
              AJOUT(dWdxIx_MOY, x*chi*dWdx);
              AJOUT(dWdyIx_MOY, x*chi*dWdy);
              AJOUT(dWdzIx_MOY, x*chi*dWdz);

              AJOUT(xUaiNx_MOY, x*u*aiNx);
              AJOUT(xVaiNx_MOY, x*v*aiNx);
              AJOUT(xWaiNx_MOY, x*w*aiNx);
              AJOUT(xUaiNy_MOY, x*u*aiNy);
              AJOUT(xVaiNy_MOY, x*v*aiNy);
              AJOUT(xWaiNy_MOY, x*w*aiNy);
              AJOUT(xUaiNz_MOY, x*u*aiNz);
              AJOUT(xVaiNz_MOY, x*v*aiNz);
              AJOUT(xWaiNz_MOY, x*w*aiNz);

              //NEW TERMS FOR THE EXTENDED PRESSURE FIELDS ((Need to avoid the invalid cells in this esteem))

              //if (chi != 0 && chi != 1){ Inserire tuti quelli definiti solo sull'interfaccia`
              //nel caso di pl vanno comunque escluse dal calcolo le celle invalide
              //}
              AJOUT(P_VAP_Iv_MOY,pv*chiv);
              //Momentum equation for the vapor phase
              AJOUT(P_VAP_aiNx_MOY,pv*aiNx); //definito solo sull'interfaccia
              AJOUT(P_VAP_aiNy_MOY,pv*aiNy); //definito solo sull'interfaccia
              AJOUT(P_VAP_aiNz_MOY,pv*aiNz); //definito solo sull'interfaccia

              AJOUT(P_LIQ_I_MOY,pl*chi);
              //Momentum equation for the liquid phase
              AJOUT(P_LIQ_aiNx_MOY,pl*aiNx); //definito solo sull'interfaccia
              AJOUT(P_LIQ_aiNy_MOY,pl*aiNy); //definito solo sull'interfaccia
              AJOUT(P_LIQ_aiNz_MOY,pl*aiNz); //definito solo sull'interfaccia
              // Reynolds stress equation (only liquid phase)
              AJOUT(UP_LIQ_aiNx_MOY, u*pl*aiNx); //definito solo sull'interfaccia
              AJOUT(VP_LIQ_aiNx_MOY, v*pl*aiNx); //definito solo sull'interfaccia
              AJOUT(WP_LIQ_aiNx_MOY, w*pl*aiNx); //definito solo sull'interfaccia
              AJOUT(UP_LIQ_aiNy_MOY, u*pl*aiNy);  //definito solo sull'interfaccia
              AJOUT(VP_LIQ_aiNy_MOY, v*pl*aiNy);  //definito solo sull'interfaccia
              AJOUT(WP_LIQ_aiNy_MOY, w*pl*aiNy); //definito solo sull'interfaccia
              AJOUT(UP_LIQ_aiNz_MOY, u*pl*aiNz); //definito solo sull'interfaccia
              AJOUT(VP_LIQ_aiNz_MOY, v*pl*aiNz); //definito solo sull'interfaccia
              AJOUT(WP_LIQ_aiNz_MOY, w*pl*aiNz); //definito solo sull'interfaccia
              //Puo essere utile definirli su tutto il dominio, come correzione ai valori precedenti??
              AJOUT(IP_LIQ_dUdx_MOY, chi*pl*dUdx);
              AJOUT(IP_LIQ_dUdy_MOY, chi*pl*dUdy);
              AJOUT(IP_LIQ_dUdz_MOY, chi*pl*dUdz);
              AJOUT(IP_LIQ_dVdx_MOY, chi*pl*dVdx);
              AJOUT(IP_LIQ_dVdy_MOY, chi*pl*dVdy);
              AJOUT(IP_LIQ_dVdz_MOY, chi*pl*dVdz);
              AJOUT(IP_LIQ_dWdx_MOY, chi*pl*dWdx);
              AJOUT(IP_LIQ_dWdy_MOY, chi*pl*dWdy);
              AJOUT(IP_LIQ_dWdz_MOY, chi*pl*dWdz);

              AJOUT(UP_LIQ_I_MOY, chi*u*pl);
              AJOUT(VP_LIQ_I_MOY, chi*v*pl);
              AJOUT(WP_LIQ_I_MOY, chi*w*pl);

              //New stats to compare all pressure terms
              //AJOUT(IdP_LIQ_dx_MOY, chi*pl*);

#undef AJOUT
            }
        }
      //rang_begin_liq_ = BEG_LIQ;
      for (int i = 0; i < nval_; i++)
        tmp(kglob, i) = moy[i] * facteur;

      // For thermal variables :
      for (int i = 0; i < nvalt_; i++)
        for (int j = 0; j < nb_thermal_fields_; j++)
          tmpt(kglob, i,j) = moyv(i,j) * facteur;
    } // end loop k

// Somme sur tous les processeurs:
  mp_sum_for_each_item(occ); // contains the number of invalid cells within the full slice
  //                            The second index is 0 for liq, 1 for vap.
  mp_sum_for_each_item(tmp);
  mp_sum_for_each_item(tmpt);

// Sur processeur 0, ajout de la contribution a l'integrale temporelle:
  if (Process::je_suis_maitre())
    {
      Nom str_occ = Nom(" ");
      Nom str_occ_vap = Nom(" ");
      for (int k = 0; k < nktot; k++)
        {
          double oc_liq = occ(k,0);
          double oc_vap = occ(k,1);
          str_occ += Nom(oc_liq)+Nom(" ");
          str_occ_vap += Nom(oc_vap)+Nom(" ");
          double facteur_liq = nijtot/(double)((nijtot)-oc_liq);
          double facteur_vap = nijtot/(double)((nijtot)-oc_vap);
          // Correct the value in tmp to account for invalid extended pressure cells :
          // BEG and END are included
          for (int i = BEG_LIQ; i <= END_LIQ; i++)
            tmp(k, i) *= facteur_liq;

          for (int i = BEG_VAP; i <= END_VAP; i++)
            tmp(k, i) *= facteur_vap;

          for (int i = 0; i < nval_; i++)
            {
              integrale_temporelle_[i][k] += tmp(k, i) * dt;  //Numero di righe e di colonne e invertito riapetto a tmp
              moyenne_spatiale_instantanee_[i][k] = tmp(k, i);
            }
          for (int i = 0; i < nvalt_; i++)
            for (int j = 0; j < nb_thermal_fields_; j++)
              {
                integrale_temporelle_temperature_(i,k,j) += tmpt(k, i,j) * dt;
                moyenne_spatiale_instantanee_temperature_(i,k,j) = tmpt(k, i,j);
              }
        }
      Cout<< "Number of invalid cells for liquid: "  << str_occ << finl;
      Cout<< "Number of invalid cells for vapor: "  << str_occ_vap << finl;
      t_integration_ += dt;
#if 0
      ArrOfDouble sumTI(nvalt_);
      sumTI = 0.;
      double sumI = 0.;
      for (int k = 0; k < nktot; k++)
        {
          sumI  +=moyenne_spatiale_instantanee_[I_MOY][k];
          for (int j = 0; j < nb_thermal_fields_; j++)
            sumTI(j) +=moyenne_spatiale_instantanee_temperature_(TI_MOY, k, j);
        }

      Cerr<< "Tl_test (bis) " << sumI/nktot << " " ;
      for (int j = 0; j < nb_thermal_fields_; j++)
        Cerr << sumTI(j)/nktot << " " ;
      Cerr << finl;
#endif
    }

#ifdef STAT_VERBOSE
// Pour verifications : impressions des normes :
  Journal() << "VERBOSE update_state: " << finl;
  Journal() << "         erreur, erreur_ddP, divU, laplP, d_divU_dx, d_divU_dy, d_divU_dz" << finl;
  Journal() << "   L1    "
            << L1_erreur << " "
            << L1_erreur_ddP << " "
            << L1_divU << " ";
  Journal() << L1_laplP << " "
            << L1_d_divU_dx << " "
            << L1_d_divU_dy << " "
            << L1_d_divU_dz << finl;
  Journal() << "   L2    "
            << L2_erreur << " "
            << L2_erreur_ddP << " "
            << L2_divU << " ";
  Journal() << L2_laplP << " "
            << L2_d_divU_dx << " "
            << L2_d_divU_dy << " "
            << L2_d_divU_dz << finl;
#endif

}

// This method is called only on master proc.
void Statistiques_dns_ijk_FT::postraiter_thermique(const double current_time) const
{
  if ((!Process::je_suis_maitre()) && (nb_thermal_fields_ !=0))
    return;

  const int nz = elem_coord_.size_array(); // nombre de points en K
  for (int ith = 0; ith < nb_thermal_fields_; ith++)
    {
      Nom n("");
      for (int flag_valeur_instantanee=0; flag_valeur_instantanee<2; flag_valeur_instantanee++)
        {
          if (ref_ijk_ft_->disable_diphasique_)
            n="monophasique_";
          else
            n="diphasique_";

          if (flag_valeur_instantanee == 0)
            n+="statistiques_thermique_";
          else
            n+="moyenne_spatiale_thermique_";

          n += Nom(ith)+Nom("_")+Nom(current_time)+Nom(".txt");
          SFichier os(n);
          os.setf(ios::scientific);
          os.precision(15);
          if (flag_valeur_instantanee == 0)
            {
              os << "# temps_integration " <<  t_integration_ << finl;
              os << "# Impression des moyennes temporelles" << finl;
            }
          else
            os << "# Impression des moyennes spatiales instantanee" << finl;

          os << "# colonne 1 : coordonnee_K" << finl;
          for (int i = 0; i < nvalt_; i++)
            os << "# colonne " << i+2 << " : " << noms_moyennes_temperature_[i] << finl;

          for (int j = 0; j < nz; j++)
            {
              char s[100];
              snprintf(s, 100, "%16.16e ", elem_coord_[j]);
              os << s;
              for (int i = 0; i < nvalt_; i++)
                {
                  double x;
                  if (flag_valeur_instantanee == 0)
                    x = integrale_temporelle_temperature_(i,j,ith) / t_integration_;
                  else
                    x = moyenne_spatiale_instantanee_temperature_(i,j,ith);
                  snprintf(s, 100, "%16.16e ", x);
                  os << s;
                }
              os << finl;
            }
        }
    }
}

void Statistiques_dns_ijk_FT::postraiter(Sortie& os, int flag_valeur_instantanee) const
{
  // Mother-class method
  Statistiques_dns_ijk::postraiter(os,flag_valeur_instantanee);
  // Completer par la thermique
  //if (nb_thermal_fields_)
  // posttraiter_thermique a besoin du temps mais sinon on aurrait pu l'appeler la.
}

// Impression des donnees pour reprise
Sortie& Statistiques_dns_ijk_FT::printOn(Sortie& os) const
{
  Statistiques_dns_ijk::printOn(os);
  return os;
}

Sortie& Statistiques_dns_ijk_FT::completer_print(Sortie& os) const
{
  if (Process::je_suis_maitre() && (nb_thermal_fields_>0))
    {
      const int nz = elem_coord_.size_array(); // nombre de points en K
      ArrOfDouble tmp_integrale_temporelle(nz);
      for (int ith = 0; ith < nb_thermal_fields_; ith++)
        {
          os  << "  thermal_fields " << ith<< " { " << "\n";
          for (int i = 0; i < nvalt_; i++)
            {
              for (int k=0; k< nz; k++)
                tmp_integrale_temporelle[k] = integrale_temporelle_temperature_(i,k,ith);

              os << noms_moyennes_temperature_[i] << " " << tmp_integrale_temporelle << finl;
            }
          os << "  }" << finl;
        }
      Cout << "Statistical thermal fields written for restart: t_integration=" << t_integration_ << finl;
    }
  return os;
}

// Reprise des donnees stat dans un fichier reprise
// Attention, cette methode peut etre appelee avant initialize() !
Entree& Statistiques_dns_ijk_FT::readOn(Entree& is)
{
  Statistiques_dns_ijk::readOn(is);
  return is;
}

int Statistiques_dns_ijk_FT::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  Cerr << "[Statistiques_dns_ijk_FT] Motcle : " << mot << finl;
  Param param(que_suis_je()+Nom("_fields"));
  int ret_val = 1;
  if (mot=="thermal_fields")
    {
      VECT(ArrOfDouble) tmp;
      tmp.dimensionner(nval_);
      int ith = -1;
      is >> ith;

      for (int i = 0; i < nvalt_; i++)
        {
          param.ajouter(noms_moyennes_temperature_[i], &tmp[i]); // que mettre.. Ne rien faire ??? il faudrait un VECT(VECT(A
        }
      param.lire_avec_accolades(is);

      // Copy the tmp into the good cell of the target :
      const int nz = tmp[0].size_array();
      // At last, we know the sizes and can dimension the tables correctly :
      if (ith ==0)
        {
          // no need to repeat it several times, and I'm afraid the resize
          // might erase the first field if no option like NO_INIT is used...
          integrale_temporelle_temperature_.resize(nvalt_,nz,nb_thermal_fields_);
          moyenne_spatiale_instantanee_temperature_.resize(nvalt_,nz,nb_thermal_fields_);
        }
      for (int i = 0; i < nvalt_; i++)
        for (int k = 0; k< nz; k++)
          integrale_temporelle_temperature_(i,k,ith) = tmp[i][k] ; // (nvalt_,nb_elem_k_tot,nb_thermal_fields_)

      Cerr << "End of reading 'thermal_fields' in Statistiques_dns_ijk_FT" << finl;
    }
  else
    {
      Cerr << "Keyword " << mot <<  " does not belong to Statistiques_dns_ijk_FT, "
           << "resorting to the mother-class for resolution" << finl;
      ret_val = Statistiques_dns_ijk::lire_motcle_non_standard(mot, is);
    }
  return ret_val;
}

void Statistiques_dns_ijk_FT::completer_read(Param& param)
{
  // At this stage, nz is 0, we haven't read the initialize.
  // Thus, we cannot resize integrale_temporelle_temperature_ or moyenne_spatiale_instantanee_temperature_ correctly.
  // However, we can do something for nb_thermal_fields_ from the problem !
  {
    int idx =0;
    for ( auto&& itr = ref_ijk_ft_->thermique_.begin(); itr != ref_ijk_ft_->thermique_.end(); ++itr )
      {
        idx++;
      }
    //curseur = ref_ijk_ft_->thermique_; //RAZ
    nb_thermal_fields_ = idx;
    Journal() <<"ooooo "<<Process::me() << finl;
    //envoyer_broadcast(nb_thermal_fields_, 0);
    Cerr << "nb_thermal_fields=" << nb_thermal_fields_<< finl;
    Journal() <<"nb_thermal_fields=" << nb_thermal_fields_<< finl;
  }

  param.ajouter_non_std("thermal_fields",(this));
  // ou pour lire quand meme les anciens directement sans accolade :
  // Ce format n'a normalement jamais existe pour la temperature.
  /* for (int i = 0; i < nvalt_; i++)
    {
      param.ajouter(noms_moyennes_temperature_[i], &integrale_temporelle_temperature_[i]); // que mettre.. Ne rien faire ??? il faudrait un VECT(VECT(A
    } */
  return;
}

// Attention, cette methode est appelee apres readOn(),
// il ne faut pas casser les donnees lues
void Statistiques_dns_ijk_FT::initialize(const IJK_FT_double& ijk_ft, const IJK_Grid_Geometry& geom)
{
  ref_ijk_ft_=ijk_ft;
  dx_=geom.get_constant_delta(0); //modif AT 20/06/2013
  dy_=geom.get_constant_delta(1); //modif AT 20/06/2013
  tab_dz_=geom.get_delta(2); //modif AT 20/06/2013
  const int n = geom.get_nb_elem_tot(DIRECTION_K);
  elem_coord_.resize_array(n);
  const ArrOfDouble& coord_z = geom.get_node_coordinates(DIRECTION_K);
  for (int i = 0; i < n; i++)
    elem_coord_[i] = (coord_z[i] + coord_z[i+1]) * 0.5;

  int idx =0;
  for ( auto&& itr = ref_ijk_ft_->thermique_.begin(); itr != ref_ijk_ft_->thermique_.end(); ++itr )
    {
      idx++;
    }
  //curseur = ref_ijk_ft_->thermique_; //RAZ
  nb_thermal_fields_ = idx;
  Journal() <<"ooooo "<< Process::me() << " with nb_thermal_fields=" << nb_thermal_fields_ << finl;
  Cerr <<"ooooo "<< Process::me() << " with nb_thermal_fields=" << nb_thermal_fields_ << finl;
  if (Process::je_suis_maitre())
    {
      moyenne_spatiale_instantanee_.dimensionner(nval_);
      for (int i = 0; i < nval_; i++)
        {
          moyenne_spatiale_instantanee_[i].resize_array(n);
        }
      //   moyenne_spatiale_instantanee_[WFACE_MOY].resize_array(n+1);//En effet, il y a une face de plus que d'elem et ce tableau est au face
      integrale_temporelle_temperature_.resize(nvalt_,n,nb_thermal_fields_);
      moyenne_spatiale_instantanee_temperature_.resize(nvalt_,n,nb_thermal_fields_);
    }

  int flag = 0;
  if (integrale_temporelle_.size() == 0)
    {
      integrale_temporelle_.dimensionner(nval_);
      t_integration_ = 0.;
    }
  else
    flag = 1;
  for (int i = 0; i < nval_; i++)
    {
      if (flag)
        {
          if (integrale_temporelle_[i].size_array() != n)
            {
              Cerr << "Erreur dans Statistiques_dns_ijk_FT::initialize: reprise avec mauvais nombre de mailles en z" << finl;
              Cerr << "integrale_temporelle_[" << i << "].size_array() = " << integrale_temporelle_[i].size_array() << finl;
              Process::exit();
            }
        }
      else
        {
          integrale_temporelle_[i].resize_array(n);

        }
    }
//modif AT 20/06/2013
  if (check_converge_)
    {
      //on recupere les vitesses moyennes aux faces pour calculer les fluctuations de vitesse.
      vit_moy_.dimensionner(3);
      for (int i = 0; i < 3; i++)
        vit_moy_[i].resize_array(n);

      vit_moy_[0]=integrale_temporelle_[0];
      vit_moy_[1]=integrale_temporelle_[1];
      vit_moy_[2]=integrale_temporelle_[40];//car on veut w aux faces !
    }
  check_stats_ = 0;
}

int Statistiques_dns_ijk_FT::initialize(const IJK_FT_double& ijk_ft, const IJK_Splitting& splitting,
                                        const int check_stats)
{
  initialize(ijk_ft,splitting.get_grid_geometry());
  // Pour les deriv de U, V et W :
  allocate_cell_vector(gradU_, splitting, 0);
  allocate_cell_vector(gradV_, splitting, 0);
  allocate_cell_vector(gradW_, splitting, 0);
  int nalloc=9;
  // veut-on verifier les derivees :
  check_stats_ = check_stats;
  if (check_stats_)
    {
      // Pour verification des stats :
      // Pour stocker les deriv secondes :
      allocate_cell_vector(grad2Pi_, splitting, 0);
      allocate_cell_vector(grad2Pc_, splitting, 0);
      allocate_cell_vector(grad2Ui_, splitting, 0);
      allocate_cell_vector(grad2Uc_, splitting, 0);
      allocate_cell_vector(grad2Vi_, splitting, 0);
      allocate_cell_vector(grad2Vc_, splitting, 0);
      allocate_cell_vector(grad2Wi_, splitting, 0);
      allocate_cell_vector(grad2Wc_, splitting, 0);
      nalloc +=24;
    }
  return nalloc;
}

const FixedVector<IJK_Field_double, 3>& Statistiques_dns_ijk_FT::get_IJK_vector_field(const Nom& nom) const
{

  if (nom== "gradU")
    return gradU_;

  if (nom== "gradV")
    return gradV_;

  if (nom== "gradW")
    return gradW_;

  if (nom== "grad2Pi")
    return grad2Pi_;

  if (nom == "grad2Pc")
    return grad2Pc_;

  if (nom== "grad2Ui")
    return grad2Ui_;

  if (nom == "grad2Uc")
    return grad2Uc_;

  if (nom== "grad2Vi")
    return grad2Vi_;

  if (nom == "grad2Vc")
    return grad2Vc_;

  if (nom== "grad2Wi")
    return grad2Wi_;

  if (nom == "grad2Wc")
    return grad2Wc_;

  Cerr << "Erreur dans Statistiques_dns_ijk_FT::get_IJK_vector_field : "
       << "Champ demande : " << nom
       << "Liste des champs possibles : "  << finl;
  Process::exit();
  throw;
}

double Statistiques_dns_ijk_FT::compute_desequil_alpha(const IJK_Grid_Geometry& geom_NS,
                                                       const double portee_wall_repulsion) const
{
  if (Process::je_suis_maitre())
    {
      const double zmin = geom_NS.get_origin(DIRECTION_K);
      const double zmax = zmin+geom_NS.get_domain_length(DIRECTION_K);
      const int n = elem_coord_.size_array(); // nombre de points en K
      const int i =0 ; // LA CASE DU TAUX DE VIDE...
      double almoy_g = 0., almoy_d = 0.;
      int ng = 0, nd = 0;
      for (int j = 0; j < n; j++)
        {
          const double z = elem_coord_[j];
          if (z - zmin < portee_wall_repulsion)
            {
              // Dans le domaine de repulsion gauche
              almoy_g +=moyenne_spatiale_instantanee_[i][j];
              ng++;
            }
          else if (zmax - z < portee_wall_repulsion )
            {
              // Dans le domaine de repulsion droite
              almoy_d +=moyenne_spatiale_instantanee_[i][j];
              nd++;
            }
        }
      if (ng>0)
        almoy_g /= ng;
      if (nd>0)
        almoy_d /= nd;
      return (almoy_g - almoy_d);

    }
  else
    {
      Cerr << " Statistiques_dns_ijk_ft::compute_desequil_alpha() should be called on master" << finl;
      Process::exit();
      return -1;
    }
}
