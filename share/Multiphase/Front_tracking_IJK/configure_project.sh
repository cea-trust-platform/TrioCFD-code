if [ "$1" = "" ]
then
    echo Usage: 
    echo  First time, to compile project:
    echo    . share/configure_project pcbenoit\|jade\|titane|env
    echo    . share/configure_project config
    echo    make debug\|optim
    echo  Next times:
    echo    . share/configure_project pcbenoit\|jade\|titane
elif [ "$1" = "config" ]
then
    $TRIO_U_ROOT/bin/baltik/bin/baltik_build_configure
    ./configure
else
    case $1 in
	pcbenoit)
	TRIO_U_ROOT=/work/mathieu/VUES/165/vobs/Trio_U
	;;
	jade)
	TRIO_U_ROOT=/work/triou/toto/Trio_U
	;;
	titane)
	TRIO_U_ROOT=/work/cont002/triou/toto/Trio_U
	;;
	env)
	;;
    esac
    pushd $TRIO_U_ROOT
    . bin/Init_Trio_U
    popd
    
    echo TRIO_U_ROOT=$TRIO_U_ROOT
    export exec=`pwd`/New_algo_qc_opt
    echo Changing exec=$exec
    echo Setting BALTIKHOSTNAME $1
    export BALTIKHOSTNAME=$1
    export BALTIKROOT=$PWD
fi
