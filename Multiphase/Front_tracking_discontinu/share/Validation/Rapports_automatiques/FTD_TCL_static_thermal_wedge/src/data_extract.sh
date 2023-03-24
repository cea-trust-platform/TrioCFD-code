#!/bin/bash
err=`ls *err` # PAR_source.err
grep "instant-Q" $err | awk '{print $5, $6, $7}' > instant_Qmic_meso.txt
grep "Volume_phase_0" $err | awk '{print $4, $2}' > vap_vol.txt
awk '{print $1*1000, 2*((3.0*$2)/(4.0*(22/7)))^(1/3)*1000}' vap_vol.txt > dia_full.txt ; awk 'NR % 100 == 0' dia_full.txt > dia.txt
grep "vevap_int" $err | awk '{print $6, $3}' > micro_vap.txt
awk '{print $1*1000, 2*((3.0*$2)/(4*(22/7)))^(1/3)*1000}' micro_vap.txt > dia_micro.txt
grep "time_sommet_tcl=" $err | awk '{print $6*1000, 2*$4*1000}' > base_dia.txt 
grep "v_cl_test=" $err | awk '{print $6, $4}' > CL_vel.txt
grep "instantaneous meso-evaporation =" $err | awk '{print $3, -$7}' > meso.txt
grep "Instantaneous tcl-evaporation =" $err | awk '{print $3, $7}' > tcl.txt
awk '{print $1, 2*((3.0*$2)/(4*(22/7)))^(1/3)}' tcl.txt > dia_tcl.txt
grep "theta_app_degree=" $err | awk '{print $9, $11}' > theta_app.txt

