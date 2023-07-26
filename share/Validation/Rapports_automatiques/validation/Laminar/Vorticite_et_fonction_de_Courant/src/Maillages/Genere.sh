#visit -cli -s Visus.py -nowin

montage mesh_cart_*.png  -tile 3x3 -geometry 1024x1024+1+1 mesh_cart.png
montage mesh_quad_*.png  -tile 3x3 -geometry 1024x1024+1+1 mesh_quad.png
montage mesh_tri_*.png  -tile 3x2 -geometry 1024x1024+1+1 mesh_tri.png
