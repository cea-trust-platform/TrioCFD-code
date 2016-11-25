rm -f *png

visit -cli -s Visus.py -nowin

montage mesh_hexa_*.png  -tile 3x2 -geometry 1024x1024+1+1 mesh_hexa.png
montage mesh_tetra_*.png  -tile 3x3 -geometry 1024x1024+1+1 mesh_tetra.png
