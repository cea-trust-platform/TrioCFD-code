//
// Warning : DO NOT EDIT !
// To update this file, run: make depend
//
#include <verifie_pere.h>
#include <Convection_Diffusion_Phase_field.h>
#include <Navier_Stokes_phase_field.h>
#include <Pb_Phase_field.h>
#include <Schema_Phase_field.h>
void instancie_src_Phase_field() {
Cerr << "src_Phase_field" << finl;
Convection_Diffusion_Phase_field inst1;verifie_pere(inst1);
Navier_Stokes_phase_field inst2;verifie_pere(inst2);
Pb_Phase_field inst3;verifie_pere(inst3);
Schema_Phase_field inst4;verifie_pere(inst4);
}
