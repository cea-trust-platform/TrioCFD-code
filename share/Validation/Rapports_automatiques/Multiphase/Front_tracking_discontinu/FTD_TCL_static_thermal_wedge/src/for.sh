source /volatile/catC/linkai/TRIOCFD/env_TrioCFD.sh

for angle in 5 10 15 20 30 50 80 90 
do
    cd CAS_ANGLE$angle
    cp ../postpro_flux.py .
    offset=$(grep offset info.txt | sed "s/offset=//")
    Nx=$(grep Nx info.txt | sed "s/Nx=//")
    angle=$(grep theta info.txt | sed "s/theta=//")
    python postpro_flux.py $offset $Nx $angle
    cd ..
done
