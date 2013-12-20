/*
 * ********+*********+*********+*********+*********+*********+*********+*******
 * By Andrew Shao
 * 9 April 2013
 * Description: Here the main tracer array
 * ********+*********+*********+*********+*********+*********+*********+*******
*/


void initialize_netcdf_vardesc ( )
{
	struct vardesc vars[NOVARS] =
	{
			{
					"D", "Basin Depth", 'h', '1', '1', "meter", 'f', 0
			}, // 0
			{
					"mn_u",
					"Zonal Velocity",
					'u',
					'L',
					's',
					"meter second-1",
					'f',
					0
			}, // 1
			{
					"mn_v",
					"Meridional Velocity",
					'v',
					'L',
					's',
					"meter second-1",
					'f',
					0
			}, // 2
			{
					"mn_h", "Layer Thickness", 'h', 'L', 's', "meter", 'd', 9.e-10
			}, // 3
			{
					"mn_uhtm",
					"zonal vel - cumul",
					'u',
					'L',
					's',
					"meter second-1",
					'f',
					0
			}, // 4
			{
					"mn_vhtm",
					"merid vel - cumul",
					'v',
					'L',
					's',
					"meter second-1",
					'f',
					0
			}, // 5
			{
					"mn_ea",
					"downward entrainment",
					'h',
					'L',
					's',
					"meters",
					'd',
					0
			}, // 6
			{
					"mn_eb", "upward entrainment", 'h', 'L', 's', "meters", 'd', 0
			}, // 7
			{
					"mn_eaml",
					"downward ML detrain",
					'h',
					'1',
					's',
					"meters",
					'd',
					0
			}, // 8
			{
					"mn_age",
					"age tracer",
					'h',
					'L',
					's',
					"years",
					'd',
					-3.17097930758233E-14
			}// 9

	};
}
