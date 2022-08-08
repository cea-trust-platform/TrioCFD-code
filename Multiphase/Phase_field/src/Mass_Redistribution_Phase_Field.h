/*
 * Mass_Redistribution_Phase_Field.h
 *
 *  Created on: May 6, 2022
 *      Author: Shambhavi Nandan
 */

#ifndef Mass_Redistribution_Phase_Field_included
#define Mass_Redistribution_Phase_Field_included

#include <Zone_VDF.h>
#include <Zone_Cl_VDF.h>
#include <Ref_Zone_Cl_VDF.h>
#include <Ref_Zone_VDF.h>
#include <Equation.h>
#include <Process.h>
#include <vector>
#include <MD_Vector_tools.h>

class Mass_Redistribution_Phase_Field
{

public:
   Mass_Redistribution_Phase_Field() = delete;
  ~Mass_Redistribution_Phase_Field() = delete;

  static void impose_mass_redistribution(const Zone_VDF&, DoubleTab&, DoubleVect, DoubleVect);
  static DoubleTab c_ini;


private:
  REF(Zone_VDF) la_zone_vdf;
  REF(Zone_Cl_VDF) la_zcl_vdf;

public:
  static double epsilon_mass_redistribute;
};



#endif /* MASS_REDISTRIBUTION_PHASE_FIELD_H_ */
