
function replace {
    file=$1
    var=$2
    text=$3
    sed "s/$2/$3/" $file >${file}tmp
    mv -f ${file}tmp $file
}

function replaceline {
    file=$1
    tag=$2
    replacementtext="$3"
    echo "$replacementtext" >${file}tmptext
    sed "/$tag/{z
r ${file}tmptext
}" $file >${file}tmp
    mv -f ${file}tmp $file
}

function writecoord {
	export N=$1
	file=$2
	export domsize=$3
	if test -z "$coord_k_formula"; then
	    coord_k_formula='z'
	fi
	expr='BEGIN{N=ENVIRON["N"];zmax=ENVIRON["domsize"]}{z=$1/N*zmax;printf "%20.16g\n",'${coord_k_formula}'}'
	#seq 0 $N | awk 'BEGIN{N=ENVIRON["N"];zmax=ENVIRON["domsize"]}{printf "%20.16g\n",$1*zmax/N}' >$file 
	seq 0 $N | awk "$expr" >$file
}

function gen_data {
    src=$1
    dest=$2
    if test -z "$nodes"; then
	nodes="17 33 17"
    fi
    if test -z "$split"; then
	split="1 1 1"
    fi
    if test -z "$splitijk"; then
	splitijk=$split
    fi
    if test -z "$ghost"; then
	ghost=2
    fi
    if test -z "$precision"; then
	precision="mixed"
    fi
    if test -z "$seuil"; then
	seuil="1e-10"
    fi
    if test -z "$iterjacobi"; then
	iterjacobi="2 5 5"
    fi
    if test -z "$iterjacobipre"; then
	iterjacobipre=$iterjacobi
    fi
    if test -z "$vx"; then
	# build an initial velocity field
	# parabolic
	vx=`echo 0.1 0.1 0.1 | awk '{x=$1
y=$2
z=$3
printf "50*y*(%f-y)*%f",y,4/(y*y)
}'`
	# add some high freq components
	for i in 1  4  16; do
	    for j in  2  8 ; do
		for k in 3 9; do
		    vx=${vx}`echo 0.1 0.1 0.1 $i $j $k | awk '{x=$1
y=$2
z=$3
fx=$4*6.2831853071795862/x
fy=$5*6.2831853071795862/y
fz=$6*6.2831853071795862/z
a=20
printf "+%f*sin(x*%f)*sin(y*%f)*sin(z*%f)",a,fx,fy,fz
}'`
		done
	    done
	done
	echo vx=$vx
    fi
    if test -z "$vy"; then
	vy="0"
    fi
    if test -z "$vz"; then
	vz="0"
    fi


    if test -z "$coarsen"; then
	coarsen="                        coarsen_operators 4
                        Coarsen_Operator_Uniform { }
                        Coarsen_Operator_Uniform { }
                        Coarsen_Operator_Uniform { coarsen_i 1 coarsen_j 1 coarsen_k 2 }"
    fi
    if test -z "$pressuresolver"; then
	pressuresolver="                solveur_pression multigrille_adrien {
			solver_precision SOLVER_PRECISION
			COARSEN_OPERATOR
			ghost_size GHOST_SIZE
			pre_smooth_steps ITERJACOBIPRE
			smooth_steps ITERJACOBI
			relax_jacobi 3 0.7 0.65 0.66
			solveur_grossier GCP { impr seuil 1e-11 precond ssor { omega 1.5 } }
			check_residu 1
			seuil SEUIL
		        impr		
			nb_full_mg_steps 2 3 1 
		}"
    fi
    cp $src $dest
    
    replaceline $dest PRESSURE_SOLVER "$pressuresolver"
    replaceline $dest COARSEN_OPERATOR "$coarsen"
    replace $dest ITERJACOBIPRE "$iterjacobipre"
    replace $dest ITERJACOBI "$iterjacobi"
    replace $dest SOLVER_PRECISION "$precision"
    replace $dest SEUIL "$seuil"
    replace $dest GHOST_SIZE "$ghost"
    replace $dest NNODES "$nodes"
    replace $dest SPLITTING_VDF "$split"
    replace $dest SPLITTING_IJK "$splitijk"
    replace $dest DOM_SIZE_X 0.1
    replace $dest DOM_SIZE_Z 0.1
    replace $dest INITIALVELOCITYX "$vx"
    replace $dest OPCONV "$opconv"
    n=`echo "$nodes" |awk '{print $2-1}'`
    echo n=$n
    writecoord $n coord_level0.txt 0.1
    writecoord $[$n/2] coord_level1.txt 0.1
    writecoord $[$n/4] coord_level2.txt 0.1
}

function extract_norm_residu {
    awk '$0~"Norme de Ax"{a=$6+1e-20;print $6,$7,$7/a}' $1 >$2
}

# for uniform mesh:
coord_k_formula='z'
# for non uniform mesh:
coord_k_formula='(z+z*z*0.1)*0.99009900990099009'
nodes="33 65 33"
#nodes="49 97 49"
if test "$coord_k_formula" != 'z'
then
  coarsen="                        coarsen_operators 3
                        Coarsen_Operator_Uniform { }
                        Coarsen_Operator_Uniform { }
                        Coarsen_Operator_K { file_z_coord coord_level1.txt }"
fi
split="2 2 2"
#split="1 1 1"
opconv=centre4ijk
mkdir tmp
mkdir tmp/lata
cd tmp
gen_data ../model.data a.data
if test "$split" == "1 1 1"; then
    echo Running sequential
    triou a 1>out 2>err
else
    n=`echo $split | awk '{print $1*$2*$3}'`
    echo Running parallel $n
    triou a $n 1>out 2>err
fi
