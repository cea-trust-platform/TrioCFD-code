# Instructions visit
DefineScalarExpression("x","coord(grid_geom2)[0]")
DefineScalarExpression("y","coord(grid_geom2)[1]")
DefineScalarExpression("z","coord(grid_geom2)[2]")
DefineScalarExpression("simu_dv_x","D_VELOCITY_X_FACES_grid_geom2_dual*point_constant(grid_geom2, 1.7828944017115125e+00)")
DefineScalarExpression("simu_dv_y","D_VELOCITY_Y_FACES_grid_geom2_dual*point_constant(grid_geom2, 1.7828944017115125e+00)")
DefineScalarExpression("simu_dv_z","D_VELOCITY_Z_FACES_grid_geom2_dual*point_constant(grid_geom2, 1.7828944017115125e+00)")
DefineScalarExpression("ana_dv_x","0.*x*y*z")
DefineScalarExpression("ana_dv_y","(-5.)*2.*0.001*6.28318530717958647*6.28318530717958647*1.8101668014267859e-05*sin(x*6.28318530717958647)*sin(z*6.28318530717958647*2)")
DefineScalarExpression("ana_dv_z","0.*x*y*z")
DefineScalarExpression("error_dv_x","D_VELOCITY_X_FACES_grid_geom2_dual*point_constant(grid_geom2, 1.7828944017115125e+00)-0.*x*y*z")
DefineScalarExpression("error_dv_y","D_VELOCITY_Y_FACES_grid_geom2_dual*point_constant(grid_geom2, 1.7828944017115125e+00)-(-5.)*2.*0.001*6.28318530717958647*6.28318530717958647*1.8101668014267859e-05*sin(x*6.28318530717958647)*sin(z*6.28318530717958647*2)")
DefineScalarExpression("error_dv_z","D_VELOCITY_Z_FACES_grid_geom2_dual*point_constant(grid_geom2, 1.7828944017115125e+00)-0.*x*y*z")

