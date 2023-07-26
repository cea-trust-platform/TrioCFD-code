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

#ifndef Problemes_rayo_included
#define Problemes_rayo_included

#include <Pb_Thermohydraulique_Turbulent_QC.h>
#include <Pb_Thermohydraulique_Turbulent.h>
#include <Pb_Hydraulique_Turbulent.h>
#include <Pb_Thermohydraulique_QC.h>
#include <Pb_Thermohydraulique.h>
#include <Pb_Hydraulique.h>
#include <Pb_Conduction.h>

//class Problemes_rayo
//{ };

class Pb_Rayo_Conduction : public Pb_Conduction
{
  Declare_instanciable(Pb_Rayo_Conduction);
};

class Pb_Rayo_Hydraulique : public Pb_Hydraulique
{
  Declare_instanciable(Pb_Rayo_Hydraulique);
};

class Pb_Rayo_Thermohydraulique : public Pb_Thermohydraulique
{
  Declare_instanciable(Pb_Rayo_Thermohydraulique);
};

class Pb_Rayo_Thermohydraulique_QC : public Pb_Thermohydraulique_QC
{
  Declare_instanciable(Pb_Rayo_Thermohydraulique_QC);
};

class Pb_Rayo_Hydraulique_Turbulent : public Pb_Hydraulique_Turbulent
{
  Declare_instanciable(Pb_Rayo_Hydraulique_Turbulent);
};

class Pb_Rayo_Thermohydraulique_Turbulent : public Pb_Thermohydraulique_Turbulent
{
  Declare_instanciable(Pb_Rayo_Thermohydraulique_Turbulent);
};

class Pb_Rayo_Thermohydraulique_Turbulent_QC : public Pb_Thermohydraulique_Turbulent_QC
{
  Declare_instanciable(Pb_Rayo_Thermohydraulique_Turbulent_QC);
};

#endif /* Problemes_rayo_included */
