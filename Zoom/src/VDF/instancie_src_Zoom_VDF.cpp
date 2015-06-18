//
// Warning : DO NOT EDIT !
// To update this file, run: make depend
//
#include <verifie_pere.h>
#include <Source_DC_VDF.h>
#include <Source_DC_VDF_NS.h>
#include <Source_LDC_VDF.h>
#include <Source_LDC_VDF_NS.h>
void instancie_src_Zoom_VDF() {
Cerr << "src_Zoom_VDF" << finl;
Source_DC_VDF inst1;verifie_pere(inst1);
Source_DC_VDF_NS inst2;verifie_pere(inst2);
Source_LDC_VDF inst3;verifie_pere(inst3);
Source_LDC_VDF_NS inst4;verifie_pere(inst4);
}
