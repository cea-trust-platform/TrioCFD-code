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
// File:        FTd_tools.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/3
//
//////////////////////////////////////////////////////////////////////////////
#ifndef FTd_tools_inclus
#define FTd_tools_inclus

class DoubleTabFT;
class DoubleTab;

typedef double FTd_matrice22[2][2];
typedef double FTd_vecteur2[2];
typedef double FTd_matrice33[3][3];
typedef double FTd_vecteur3[3];


double FTd_calculer_max_norme_vecteurs(const DoubleTab& vecteur);
double FTd_calculer_aire_triangle(const FTd_vecteur2& coord_som0,const FTd_vecteur2& coord_som1,const FTd_vecteur2& coord_som2);
double FTd_calculer_volume_tetraedre(const FTd_vecteur3& coord_som0,const FTd_vecteur3& coord_som1,const FTd_vecteur3& coord_som2,const FTd_vecteur3& coord_som3);
int FTd_compare_sommets_global(int pe0, int numOwner0, int pe1, int numOwner1);

#endif
