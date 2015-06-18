//
// Warning : DO NOT EDIT !
// To update this file, run: make depend
//
#include <verifie_pere.h>
#include <Paroi_DWF_hyd_VDF.h>
#include <Source_DWF_VDF.h>
void instancie_src_Zoom_VDF_Turbulence() {
Cerr << "src_Zoom_VDF_Turbulence" << finl;
Paroi_DWF_hyd_VDF inst1;verifie_pere(inst1);
Source_DWF_VDF inst2;verifie_pere(inst2);
}
