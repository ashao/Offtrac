void read_var3d( char inpath[200], char varname[200], int imon, double*** data);
void read_fieldave(int imon);
void read_D();
void read_uvw(int imon,int itts);
void read_h(int imon,int inxt);
void read_fields(int imon,int itts);
void read_biotic_bc(int imon,int itts);
void read_ts(int imon,int itts);
# ifdef REST_ARCTIC
void read_restoring(int imon,int itts);
# endif
void read_clim(int imon,int inxt,int ilst);
#if defined JO2_SPEC || defined UNIVENT
void read_misc_fields();
#endif
#ifdef RESTART
void read_tracer_init(int imon, char *run_name);
#else
void read_tracer_init(int imon);
#endif
void read_buoy(int imon);
void read_sponge(void);
void read_grid();
#ifdef LEV_OXY
void read_oxy_ic(void);
#endif
//begin ashao
# if defined(CFC) || defined(SF6)
int read_trac(double sc_coeffs[4], double sol_coeffs[7], int tracyear[MAXRECS], double Ncon[MAXRECS], double Scon[MAXRECS]);
# endif
