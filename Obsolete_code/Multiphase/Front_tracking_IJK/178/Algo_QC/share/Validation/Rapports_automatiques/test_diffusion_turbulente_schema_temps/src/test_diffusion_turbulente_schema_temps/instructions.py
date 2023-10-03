# Instructions visit
DefineScalarExpression("x","coord(grid_geom2)[0]")
DefineScalarExpression("y","coord(grid_geom2)[1]")
DefineScalarExpression("z","coord(grid_geom2)[2]")
DefineScalarExpression("simu_vz","VELOCITY_Z_FACES_grid_geom2_dual")
DefineScalarExpression("ana_vz_tinit","2.*0.001*sin(x*6.28318530717958647*8.)*sin(z*6.28318530717958647*8.)")
DefineScalarExpression("ana_vz_tfinal","0.58653381*2.*0.001*sin(x*6.28318530717958647*8.)*sin(z*6.28318530717958647*8.)")
DefineScalarExpression("error_vz","VELOCITY_Z_FACES_grid_geom2_dual-0.58653381*2.*0.001*sin(x*6.28318530717958647*8.)*sin(z*6.28318530717958647*8.)")

