function run_calibration()

%[panel, geometry, phantom, objOffset, optimization, data, spAttrb] = calibration_setup_upenngtr4_RADA()
[panel, geometry, phantom, objOffset, optimization, data, spAttrb] = calibration_setup_upenngtr4_RADB()
[measurements,geometry] = cbct_geometric_calibration(panel, geometry, phantom, objOffset, optimization, data, spAttrb)