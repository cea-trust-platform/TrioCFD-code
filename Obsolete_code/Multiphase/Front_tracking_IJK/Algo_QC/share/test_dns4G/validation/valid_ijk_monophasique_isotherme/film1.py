from math import sin, cos, pi
from render import *
from wall_tools import *
import os

snorky_root=os.environ.get("SNORKY_ROOT")

init_global_thread_pool(7)

set_Journal_level(1)

npartmax=200000
timefact=0.5
engine=RenderTool(0)
camera=Camera()

camera.set_pos(0.2,0.01,-0.7);
camera.set_rot(0,0);
camera.set_scale(0.1);

screensize=0.3
camera.set_stereo_params(4.*screensize/3., screensize, 0.06, 2.*screensize/3., 10.*screensize/3.);
camera.set_pos(-0.4,0.75*0.2,-0.7);
camera.set_rot(0,90);
camera.set_scale(0.2);
#camera.set_swap_yzrotation(1)
engine.set_camera(camera)

db=LataFilterTool()
db.open('toto/titi.lata','')

list=RenderObjetsReferences()

cred  =[0.4,0,0,1,1]
cgreen=[0.4,1,1,1,0]
cblue =[1,1,0,0,0]
colormap=FloatTab(5, 4)
for i in range(5):
    colormap.set(i,0,cred[i])
    colormap.set(i,1,cgreen[i])
    colormap.set(i,2,cblue[i])
    colormap.set(i,3,1)
    pass
font=FontTool()
font.load_pgm_file(snorky_root+"/courier_bold.pgm")
partmode=1
if (1):
    part=ParticlesTool()
    part.init(db, "DOM", "VELOCITY", "FACES", npartmax)
    #part.setScalarColormap("TEMPERATURE", "ELEM", 0.6, 0.9, colormap);
    part.setPointSize(0.003)
    part.setLightFactor(2)
    part.set_cache_all_fields(0)
    part.set_repeat_database(1)
    if (partmode==1):
        part.setMaxlifetime(0.001,100,1e-6)
        part.set_particles_per_second(1500*npartmax)
    else:
        part.setMaxlifetime(3e-4,1e-6,0)
        part.set_particles_per_second(1e9)
        pass
    part.set_integration_step(6e-6*timefact)
    part.set_motion_blur_duration(6e-6)
    
    srcs=[]
    srcs.append(ParticleSourceBox())
    if (partmode==1):
        #srcs[0].setBox(0, 0, 0.75, 0.05, 1, 0.75)
        srcs[0].setBox(0, 0.1, 0., 0.05, 0.1, 1.5)
    else:
        srcs[0].setBox(0, 0.1, 0., 4., 0.1, 1.5)
        pass
    srcs[0].setWeight(1)
    # srcs[0].setGrid(1000,1000,5)
    part.add_source_ref(srcs[0])
    #srcs.append(ParticleSourceBox())
    #srcs[1].setBox(0,0,1, 0.05,1,1)
    #srcs[1].setWeight(0.06)
    #part.add_source_ref(srcs[1])
        
    
    list.push_back(part)
##     stats=TimeStatDisplay()
##     stats.set_particles(part);
##     stats.set_text(font, "TOTO")
##     stats.set_display_mode(TextRenderer.DISPLAY_MODE_2D_UNIT)
##     stats.set_size(0.02)
##     stats.set_color(1,1,1,1)
##     stats.set_translation(0.,0.)
##     list.push_back(stats);
    pass

if 0:
    vdfmesh = StructuredMeshRenderer()
    vdfmesh.set_geometry(db, 'dom_IJK')
    list.push_back(vdfmesh)
    pass

engine.set_speed_ratio(1e-6*50*timefact)

engine.run_interactive(list)

##streamencode='ffmpeg'
##if streamencode=='ffmpeg':
##    target=screenfilm
##    tilesize_x=6720
##    tilesize_y=3776
##    stereo_side_by_side=1
##else:
##    target=cea_wall_saclay
##    tilesize_x = 10080/2
##    tilesize_y = 5688
##    cea_wall_saclay.supersampling_=3
##    stereo_side_by_side=0
##    pass
##filename='film'
##engine.set_srgb_brightness(1)
##engine.set_offline_tile_size(tilesize_x, tilesize_y, ColorFloat16)
##filename=prepare_film(engine, target, filename, streamencode)

##engine.set_offline_render_size(target.xres_, target.yres_, target.supersampling_)
##engine.set_movie_prefix(filename)
##engine.set_movie_framerate(25)
##engine.set_display_offline(1)
##engine.set_offline_stereo(1,stereo_side_by_side)
##engine.run_offline_render(60*25,list)
