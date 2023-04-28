# Instructions visit
DefineScalarExpression("x","coord(grid_geom2)[0]")
DefineScalarExpression("y","coord(grid_geom2)[1]")
DefineScalarExpression("z","coord(grid_geom2)[2]")
DefineScalarExpression("simu_t","TEMPERATURE_ELEM_grid_geom2")
DefineScalarExpression("ana_t","(293.+293.*z)")
DefineScalarExpression("error_t","TEMPERATURE_ELEM_grid_geom2-(293.+293.*z)")
