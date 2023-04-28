#!/bin/bash

echo 1 | awk -v Pr=0.025 '{for (i=1; i<390;i++) print i " " Pr*i*exp(-1e-2*(Pr*i)^4/(1+5*Pr^3*i))+(2.12*log((1+i)*(1))+(3.85*Pr^0.333-1.3)^2+2.12*log(Pr))*exp(-(1+5*Pr^3*i)/1e-2/(Pr*i)^4);}'
