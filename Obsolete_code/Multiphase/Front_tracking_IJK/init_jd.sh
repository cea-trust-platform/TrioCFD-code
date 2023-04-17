if test "$HOSTNAME" == "is221705"
then
  . /home/jd235850/Trio_U/v1.6.6/Trio_U/bin/Init_Trio_U
elif test "$HOSTNAME" == "gre058509"
then
  .  /work/jd235850/v166/vobs/Trio_U/bin/Init_Trio_U
elif test "$HOSTNAME" == "eris-ib.cluster"
then
  source ~triou/env_triou_1.6.6.sh 
else
  echo HOSTNAME inconnu
  # ne pas faire exit car le script est lance avec . espace, donc le script est  #lance dans le shell ou on a lance le script. avec ./ le script est lance dans un sous shell, donc le exit permet de sortir de ce sous shell
fi
export execdb=$PWD/New_algo_qc
export execopt=$PWD/New_algo_qc_opt
export exec=$execdb
if test "$HOSTNAME" == "eris-ib.cluster"
then
  export exec=$execopt
fi
