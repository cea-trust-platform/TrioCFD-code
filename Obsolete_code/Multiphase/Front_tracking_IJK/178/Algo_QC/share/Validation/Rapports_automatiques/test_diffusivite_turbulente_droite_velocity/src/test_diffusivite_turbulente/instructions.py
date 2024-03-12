# Instructions visit
DefineScalarExpression("x","coord(grid_geom2)[0]")
DefineScalarExpression("y","coord(grid_geom2)[1]")
DefineScalarExpression("z","coord(grid_geom2)[2]")
DefineScalarExpression("simu_unst","point_constant(grid_geom2, 1.0000000000000000e+00)/TEMPERATURE_ELEM_grid_geom2")
DefineScalarExpression("ana_unst","(1./293.-.5/293.*z)")
DefineScalarExpression("error_unst","point_constant(grid_geom2, 1.0000000000000000e+00)/TEMPERATURE_ELEM_grid_geom2-(1./293.-.5/293.*z)")
