/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__glutamate
#define _nrn_initial _nrn_initial__glutamate
#define nrn_cur _nrn_cur__glutamate
#define _nrn_current _nrn_current__glutamate
#define nrn_jacob _nrn_jacob__glutamate
#define nrn_state _nrn_state__glutamate
#define _net_receive _net_receive__glutamate 
#define state state__glutamate 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gNMDAmax _p[0]
#define gNMDAmax_columnindex 0
#define gAMPAmax _p[1]
#define gAMPAmax_columnindex 1
#define e _p[2]
#define e_columnindex 2
#define dend _p[3]
#define dend_columnindex 3
#define pos _p[4]
#define pos_columnindex 4
#define locx _p[5]
#define locx_columnindex 5
#define locy _p[6]
#define locy_columnindex 6
#define stim _p[7]
#define stim_columnindex 7
#define inmda _p[8]
#define inmda_columnindex 8
#define iampa _p[9]
#define iampa_columnindex 9
#define gnmda _p[10]
#define gnmda_columnindex 10
#define local_v _p[11]
#define local_v_columnindex 11
#define tt _p[12]
#define tt_columnindex 12
#define A _p[13]
#define A_columnindex 13
#define B _p[14]
#define B_columnindex 14
#define gampa _p[15]
#define gampa_columnindex 15
#define DA _p[16]
#define DA_columnindex 16
#define DB _p[17]
#define DB_columnindex 17
#define Dgampa _p[18]
#define Dgampa_columnindex 18
#define _g _p[19]
#define _g_columnindex 19
#define _nd_area  *_ppvar[0]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 /* declaration of user functions */
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(Object* _ho) { void* create_point_process(int, Object*);
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt(void*);
 static double _hoc_loc_pnt(void* _vptr) {double loc_point_process(int, void*);
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(void* _vptr) {double has_loc_point(void*);
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(void* _vptr) {
 double get_loc_point_process(void*); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 0,0
};
 static Member_func _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 0, 0
};
 /* declare global and static user variables */
#define Vset Vset_glutamate
 double Vset = -60;
#define Voff Voff_glutamate
 double Voff = 0;
#define gama gama_glutamate
 double gama = 0.08;
#define n n_glutamate
 double n = 0.25;
#define tau2 tau2_glutamate
 double tau2 = 2;
#define tau1 tau1_glutamate
 double tau1 = 50;
#define tau_ampa tau_ampa_glutamate
 double tau_ampa = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "tau1_glutamate", "ms",
 "tau2_glutamate", "ms",
 "tau_ampa_glutamate", "ms",
 "n_glutamate", "/mM",
 "gama_glutamate", "/mV",
 "gNMDAmax", "nS",
 "gAMPAmax", "nS",
 "e", "mV",
 "A", "nS",
 "B", "nS",
 "gampa", "nS",
 "inmda", "nA",
 "iampa", "nA",
 "gnmda", "nS",
 0,0
};
 static double A0 = 0;
 static double B0 = 0;
 static double delta_t = 0.01;
 static double gampa0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "tau1_glutamate", &tau1_glutamate,
 "tau2_glutamate", &tau2_glutamate,
 "tau_ampa_glutamate", &tau_ampa_glutamate,
 "n_glutamate", &n_glutamate,
 "gama_glutamate", &gama_glutamate,
 "Voff_glutamate", &Voff_glutamate,
 "Vset_glutamate", &Vset_glutamate,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void nrn_cur(NrnThread*, _Memb_list*, int);
static void  nrn_jacob(NrnThread*, _Memb_list*, int);
 static void _hoc_destroy_pnt(void* _vptr) {
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(NrnThread*, _Memb_list*, int);
static void _ode_matsol(NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[2]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"glutamate",
 "gNMDAmax",
 "gAMPAmax",
 "e",
 "dend",
 "pos",
 "locx",
 "locy",
 "stim",
 0,
 "inmda",
 "iampa",
 "gnmda",
 "local_v",
 "tt",
 0,
 "A",
 "B",
 "gampa",
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 20, _prop);
 	/*initialize range parameters*/
 	gNMDAmax = 0;
 	gAMPAmax = 1;
 	e = 0;
 	dend = 0;
 	pos = 0;
 	locx = 0;
 	locy = 0;
 	stim = 0;
  }
 	_prop->param = _p;
 	_prop->param_size = 20;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _glutamate_reg() {
	int _vectorized = 0;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 0,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 20, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 glutamate glutamate.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double F = 96480.0;
 static double R = 8.314;
 
#define PI _nrnunit_PI[_nrnunit_use_legacy_]
static double _nrnunit_PI[2] = {0x1.921fb54442d18p+1, 3.14159}; /* 3.14159265358979312 */
static int _reset;
static char *modelname = "Glutamatergic synapse with network activation";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[3], _dlist1[3];
 static int state(_threadargsproto_);
 extern int state_discon_flag_;
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   DA = - A / tau1 ;
   DB = - B / tau2 ;
   Dgampa = - gampa / tau_ampa ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 DA = DA  / (1. - dt*( ( - 1.0 ) / tau1 )) ;
 DB = DB  / (1. - dt*( ( - 1.0 ) / tau2 )) ;
 Dgampa = Dgampa  / (1. - dt*( ( - 1.0 ) / tau_ampa )) ;
  return 0;
}
 /*END CVODE*/
 static int state () {_reset=0;
 {
    A = A + (1. - exp(dt*(( - 1.0 ) / tau1)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau1 ) - A) ;
    B = B + (1. - exp(dt*(( - 1.0 ) / tau2)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau2 ) - B) ;
    gampa = gampa + (1. - exp(dt*(( - 1.0 ) / tau_ampa)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau_ampa ) - gampa) ;
   }
  return 0;
}
 
static int _ode_count(int _type){ return 3;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 ();
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 3; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  A = A0;
  B = B0;
  gampa = gampa0;
 {
   gnmda = 0.0 ;
   gampa = 0.0 ;
   A = 0.0 ;
   B = 0.0 ;
   stim = 0.0 ;
   tt = 0.0 ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel();
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   if ( ( stim  == 1.0 )  && ( tt < t ) ) {
     state_discontinuity ( _cvode_ieq + 0, & A , A + gNMDAmax ) ;
     state_discontinuity ( _cvode_ieq + 1, & B , B + gNMDAmax ) ;
     state_discontinuity ( _cvode_ieq + 2, & gampa , gampa + gAMPAmax ) ;
     tt = t + 2.0 ;
     }
   local_v = v * ( 1.0 - Voff ) + Vset * Voff ;
   gnmda = ( A - B ) / ( 1.0 + n * exp ( - gama * local_v ) ) ;
   inmda = ( 1e-3 ) * gnmda * ( v - e ) ;
   iampa = ( 1e-3 ) * gampa * ( v - e ) ;
   }
 _current += inmda;
 _current += iampa;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _g = _nrn_current(_v + .001);
 	{ state_discon_flag_ = 1; _rhs = _nrn_current(_v); state_discon_flag_ = 0;
 	}
 _g = (_g - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
 { error =  state();
 if(error){fprintf(stderr,"at line 92 in file glutamate.mod:\n	SOLVE state METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 }}}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = A_columnindex;  _dlist1[0] = DA_columnindex;
 _slist1[1] = B_columnindex;  _dlist1[1] = DB_columnindex;
 _slist1[2] = gampa_columnindex;  _dlist1[2] = Dgampa_columnindex;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "glutamate.mod";
static const char* nmodl_file_text = 
  "COMMENT\n"
  "//******************************************//\n"
  "// Created by Alon Poleg-Polsky 			//\n"
  "//    alon.poleg-polsky@ucdenver.edu		//\n"
  "//		2018								//\n"
  "//******************************************//\n"
  "ENDCOMMENT\n"
  "\n"
  "TITLE Glutamatergic synapse with network activation\n"
  "\n"
  "NEURON {\n"
  "	POINT_PROCESS glutamate\n"
  "	NONSPECIFIC_CURRENT inmda, iampa\n"
  "	RANGE e ,gAMPAmax, gNMDAmax, inmda, iampa\n"
  "\n"
  "	RANGE gnmda, gampa, dend, pos, locx, locy, local_v\n"
  "	RANGE stim, tt\n"
  "\n"
  "	GLOBAL n, gama, tau_ampa\n"
  "	:GLOBAL Pr\n"
  "	GLOBAL tau1, tau2\n"
  "	GLOBAL Voff, Vset\n"
  "\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(nA) 	= (nanoamp)\n"
  "	(mV)	= (millivolt)\n"
  "	(nS) 	= (nanomho)\n"
  "	(mM)    = (milli/liter)\n"
  "    F		= 96480 (coul)\n"
  "    R       = 8.314 (volt-coul/degC)\n"
  " 	PI = (pi) (1)\n"
  "	(mA) = (milliamp)\n"
  "	(um) = (micron)\n"
  "\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	gNMDAmax=	0	(nS)\n"
  "	gAMPAmax=	1	(nS)\n"
  "	e= 			0	(mV)\n"
  "	tau1=		50	(ms)	\n"
  "	tau2=		2	(ms)	\n"
  "	tau_ampa=	1	(ms)	\n"
  "	n=			0.25 	(/mM)	\n"
  "	gama=		0.08 	(/mV) \n"
  "	dt (ms)\n"
  "	v		(mV)\n"
  "	dend=		0\n"
  "	pos=		0\n"
  "	locx=		0\n"
  "	locy=		0\n"
  "	:Pr=		1\n"
  "	Voff=		0		:0 - voltage dependent 1- voltage independent\n"
  "	Vset=		-60		:set voltage when voltage independent		\n"
  "	stim=		0	\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	inmda		(nA)  \n"
  "	iampa		(nA)  \n"
  "	gnmda		(nS)\n"
  "	local_v\n"
  "	tt\n"
  "}\n"
  "STATE {\n"
  "	A 		(nS)\n"
  "	B 		(nS)\n"
  "	gampa 	(nS)\n"
  "\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "    gnmda=	0 \n"
  "    gampa=	0 \n"
  "	A=		0\n"
  "	B=		0\n"
  "	stim=	0\n"
  "	tt=		0\n"
  "}    \n"
  "\n"
  "BREAKPOINT {  \n"
  "	if((stim == 1) && (tt < t)){\n"
  "		:if(scop_random()<=Pr){\n"
  "			state_discontinuity( A, A + gNMDAmax)\n"
  "			state_discontinuity( B, B + gNMDAmax)\n"
  "			state_discontinuity( gampa, gampa + gAMPAmax)\n"
  "			tt= t + 2\n"
  "		:}\n"
  "	}\n"
  "	SOLVE state METHOD cnexp\n"
  "	local_v= v * (1 - Voff) + Vset * Voff	:temp voltage\n"
  "	gnmda= (A - B) / (1 + n * exp(-gama * local_v) )\n"
  "	inmda= (1e-3) * gnmda * (v-e)\n"
  "	iampa= (1e-3) * gampa * (v- e)\n"
  "}\n"
  "\n"
  "DERIVATIVE state {\n"
  "	A'= -A / tau1\n"
  "	B'= -B / tau2\n"
  "	gampa'= -gampa / tau_ampa\n"
  "}\n"
  ;
#endif
