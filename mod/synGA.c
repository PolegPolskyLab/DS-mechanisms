/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
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
 
#define nrn_init _nrn_init__synGA
#define _nrn_initial _nrn_initial__synGA
#define nrn_cur _nrn_cur__synGA
#define _nrn_current _nrn_current__synGA
#define nrn_jacob _nrn_jacob__synGA
#define nrn_state _nrn_state__synGA
#define _net_receive _net_receive__synGA 
#define states states__synGA 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define g_linMax _p[0]
#define g_linMax_columnindex 0
#define g _p[1]
#define g_columnindex 1
#define g_vgMax _p[2]
#define g_vgMax_columnindex 2
#define i_lin _p[3]
#define i_lin_columnindex 3
#define i_vg _p[4]
#define i_vg_columnindex 4
#define g_lin _p[5]
#define g_lin_columnindex 5
#define g_vg _p[6]
#define g_vg_columnindex 6
#define ninf _p[7]
#define ninf_columnindex 7
#define A _p[8]
#define A_columnindex 8
#define B _p[9]
#define B_columnindex 9
#define DA _p[10]
#define DA_columnindex 10
#define DB _p[11]
#define DB_columnindex 11
#define v _p[12]
#define v_columnindex 12
#define _g _p[13]
#define _g_columnindex 13
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
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 /* declaration of user functions */
 static double _hoc_calc_ninf(void*);
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
 _extcall_prop = _prop;
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
 "calc_ninf", _hoc_calc_ninf,
 0, 0
};
#define calc_ninf calc_ninf_synGA
 extern double calc_ninf( _threadargsprotocomma_ double , double , double , double );
 /* declare global and static user variables */
#define e e_synGA
 double e = 0;
#define ninf15 ninf15_synGA
 double ninf15 = 0;
#define ninf25 ninf25_synGA
 double ninf25 = 0;
#define ninf35 ninf35_synGA
 double ninf35 = 0;
#define ninf45 ninf45_synGA
 double ninf45 = 0;
#define ninf55 ninf55_synGA
 double ninf55 = 0;
#define ninf65 ninf65_synGA
 double ninf65 = 0;
#define ninf75 ninf75_synGA
 double ninf75 = 0;
#define ninf10 ninf10_synGA
 double ninf10 = 0;
#define ninf20 ninf20_synGA
 double ninf20 = 0;
#define ninf30 ninf30_synGA
 double ninf30 = 0;
#define ninf40 ninf40_synGA
 double ninf40 = 0;
#define ninf50 ninf50_synGA
 double ninf50 = 0;
#define ninf60 ninf60_synGA
 double ninf60 = 0;
#define ninf70 ninf70_synGA
 double ninf70 = 0;
#define ninf80 ninf80_synGA
 double ninf80 = 0;
#define tauB tauB_synGA
 double tauB = 100;
#define tauA tauA_synGA
 double tauA = 10;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "tauA_synGA", "ms",
 "tauB_synGA", "ms",
 "e_synGA", "mV",
 "g_linMax", "nS",
 "g", "nS",
 "g_vgMax", "nS",
 "A", "nS",
 "B", "nS",
 "i_lin", "nA",
 "i_vg", "nA",
 0,0
};
 static double A0 = 0;
 static double B0 = 0;
 static double delta_t = 0.01;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "tauA_synGA", &tauA_synGA,
 "tauB_synGA", &tauB_synGA,
 "ninf80_synGA", &ninf80_synGA,
 "ninf70_synGA", &ninf70_synGA,
 "ninf60_synGA", &ninf60_synGA,
 "ninf50_synGA", &ninf50_synGA,
 "ninf40_synGA", &ninf40_synGA,
 "ninf30_synGA", &ninf30_synGA,
 "ninf20_synGA", &ninf20_synGA,
 "ninf10_synGA", &ninf10_synGA,
 "ninf75_synGA", &ninf75_synGA,
 "ninf65_synGA", &ninf65_synGA,
 "ninf55_synGA", &ninf55_synGA,
 "ninf45_synGA", &ninf45_synGA,
 "ninf35_synGA", &ninf35_synGA,
 "ninf25_synGA", &ninf25_synGA,
 "ninf15_synGA", &ninf15_synGA,
 "e_synGA", &e_synGA,
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
"synGA",
 "g_linMax",
 "g",
 "g_vgMax",
 0,
 "i_lin",
 "i_vg",
 "g_lin",
 "g_vg",
 "ninf",
 0,
 "A",
 "B",
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
 	_p = nrn_prop_data_alloc(_mechtype, 14, _prop);
 	/*initialize range parameters*/
 	g_linMax = 0.1;
 	g = 0;
 	g_vgMax = 0.1;
  }
 	_prop->param = _p;
 	_prop->param_size = 14;
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

 void _synGA_reg() {
	int _vectorized = 1;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 1,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 14, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 synGA synGA.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[2], _dlist1[2];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset = 0; {
   DA = - A / tauA ;
   DB = - B / tauB ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
 DA = DA  / (1. - dt*( ( - 1.0 ) / tauA )) ;
 DB = DB  / (1. - dt*( ( - 1.0 ) / tauB )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) { {
    A = A + (1. - exp(dt*(( - 1.0 ) / tauA)))*(- ( 0.0 ) / ( ( - 1.0 ) / tauA ) - A) ;
    B = B + (1. - exp(dt*(( - 1.0 ) / tauB)))*(- ( 0.0 ) / ( ( - 1.0 ) / tauB ) - B) ;
   }
  return 0;
}
 
double calc_ninf ( _threadargsprotocomma_ double _ln_low , double _ln_high , double _lv_low , double _lv_high ) {
   double _lcalc_ninf;
 _lcalc_ninf = ( ( v - _lv_low ) * _ln_high - ( v - _lv_high ) * _ln_low ) / ( _lv_high - _lv_low ) ;
   
return _lcalc_ninf;
 }
 
static double _hoc_calc_ninf(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  calc_ninf ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) );
 return(_r);
}
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 (_p, _ppvar, _thread, _nt);
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  A = A0;
  B = B0;
 
}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 initmodel(_p, _ppvar, _thread, _nt);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   g_lin = g * g_linMax ;
   i_lin = ( 1e-3 ) * g_lin * ( v - e ) ;
   A = A + g ;
   B = B + g ;
   if ( v <= - 80.0 ) {
     ninf = ninf80 ;
     }
   else {
     if ( v <= - 75.0 ) {
       ninf = calc_ninf ( _threadargscomma_ ninf80 , ninf75 , - 80.0 , - 75.0 ) ;
       }
     else {
       if ( v <= - 70.0 ) {
         ninf = calc_ninf ( _threadargscomma_ ninf75 , ninf70 , - 75.0 , - 70.0 ) ;
         }
       else {
         if ( v <= - 65.0 ) {
           ninf = calc_ninf ( _threadargscomma_ ninf70 , ninf65 , - 70.0 , - 65.0 ) ;
           }
         else {
           if ( v <= - 60.0 ) {
             ninf = calc_ninf ( _threadargscomma_ ninf65 , ninf60 , - 65.0 , - 60.0 ) ;
             }
           else {
             if ( v <= - 55.0 ) {
               ninf = calc_ninf ( _threadargscomma_ ninf60 , ninf55 , - 60.0 , - 55.0 ) ;
               }
             else {
               if ( v <= - 50.0 ) {
                 ninf = calc_ninf ( _threadargscomma_ ninf55 , ninf50 , - 55.0 , - 50.0 ) ;
                 }
               else {
                 if ( v <= - 45.0 ) {
                   ninf = calc_ninf ( _threadargscomma_ ninf50 , ninf45 , - 50.0 , - 45.0 ) ;
                   }
                 else {
                   if ( v <= - 40.0 ) {
                     ninf = calc_ninf ( _threadargscomma_ ninf45 , ninf40 , - 45.0 , - 40.0 ) ;
                     }
                   else {
                     if ( v <= - 35.0 ) {
                       ninf = calc_ninf ( _threadargscomma_ ninf40 , ninf35 , - 40.0 , - 35.0 ) ;
                       }
                     else {
                       if ( v <= - 30.0 ) {
                         ninf = calc_ninf ( _threadargscomma_ ninf35 , ninf30 , - 35.0 , - 30.0 ) ;
                         }
                       else {
                         if ( v <= - 25.0 ) {
                           ninf = calc_ninf ( _threadargscomma_ ninf30 , ninf25 , - 30.0 , - 25.0 ) ;
                           }
                         else {
                           if ( v <= - 20.0 ) {
                             ninf = calc_ninf ( _threadargscomma_ ninf25 , ninf20 , - 25.0 , - 20.0 ) ;
                             }
                           else {
                             if ( v <= - 15.0 ) {
                               ninf = calc_ninf ( _threadargscomma_ ninf20 , ninf15 , - 20.0 , - 15.0 ) ;
                               }
                             else {
                               if ( v <= - 10.0 ) {
                                 ninf = calc_ninf ( _threadargscomma_ ninf15 , ninf10 , - 15.0 , - 10.0 ) ;
                                 }
                               else {
                                 ninf = ninf10 ;
                                 }
                               }
                             }
                           }
                         }
                       }
                     }
                   }
                 }
               }
             }
           }
         }
       }
     }
   g_vg = ( B - A ) * ninf * g_vgMax ;
   if ( g_vg < 0.0 ) {
     g_vg = 0.0 ;
     }
   i_vg = ( 1e-3 ) * g_vg * ( v - e ) ;
   }
 _current += i_lin;
 _current += i_vg;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
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
 
}
 
}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 
}
 
}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 {   states(_p, _ppvar, _thread, _nt);
  }}}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = A_columnindex;  _dlist1[0] = DA_columnindex;
 _slist1[1] = B_columnindex;  _dlist1[1] = DB_columnindex;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "synGA.mod";
static const char* nmodl_file_text = 
  "\n"
  "NEURON {\n"
  "	POINT_PROCESS synGA\n"
  "	:SUFFIX synGA\n"
  "	NONSPECIFIC_CURRENT i_lin, i_vg\n"
  "	RANGE g_linMax, g_vg, g_vgMax, g_vg, g, g_lin\n"
  "	RANGE ninf \n"
  "	GLOBAL ninf80, ninf70,ninf60, ninf50, ninf40, ninf30, ninf20, ninf10\n"
  "	GLOBAL ninf75, ninf65, ninf55, ninf45, ninf35, ninf25, ninf15\n"
  "	GLOBAL tauA, tauB\n"
  "	GLOBAL e\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "		:linear\n"
  "	g_linMax= 0.1	(nS)\n"
  "	g= 0		(nS)\n"
  "		:voltage gated\n"
  "	g_vgMax= 0.1	(nS)\n"
  "	tauA= 10		(ms)\n"
  "	tauB= 100		(ms)\n"
  "	ninf80= 0\n"
  "	ninf70= 0\n"
  "	ninf60= 0\n"
  "	ninf50= 0\n"
  "	ninf40= 0\n"
  "	ninf30= 0\n"
  "	ninf20= 0\n"
  "	ninf10= 0\n"
  "	ninf75= 0\n"
  "	ninf65= 0\n"
  "	ninf55= 0\n"
  "	ninf45= 0\n"
  "	ninf35= 0\n"
  "	ninf25= 0\n"
  "	ninf15= 0\n"
  "\n"
  "		:params\n"
  "	v 				(mV)\n"
  "	e= 0			(mV)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	i_lin 			(nA)		\n"
  "	i_vg			(nA)\n"
  "	g_lin\n"
  "	g_vg	\n"
  "	ninf\n"
  "	\n"
  "}\n"
  "\n"
  "STATE { 	\n"
  "	A 		(nS)\n"
  "	B 		(nS) \n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "		:linear\n"
  "	g_lin= g * g_linMax\n"
  "	\n"
  "	i_lin= (1e-3) * g_lin * (v - e)\n"
  "\n"
  "		:voltage gated\n"
  "	:state_discontinuity(A, A + g_vgMax)\n"
  "	:state_discontinuity(B, B + g_vgMax)\n"
  "	A= A + g\n"
  "	B= B + g\n"
  "    SOLVE states METHOD cnexp\n"
  "	if(v <= -80){\n"
  "		ninf= ninf80\n"
  "	}else{\n"
  "		if(v <= -75){\n"
  "			ninf= calc_ninf(ninf80,ninf75,-80,-75)\n"
  "		}else{\n"
  "			if(v <= -70){\n"
  "				ninf= calc_ninf(ninf75,ninf70,-75,-70)\n"
  "			}else{\n"
  "				if(v <= -65){\n"
  "					ninf= calc_ninf(ninf70,ninf65,-70,-65)\n"
  "				}else{				\n"
  "					if(v <= -60){\n"
  "						ninf= calc_ninf(ninf65,ninf60,-65,-60)\n"
  "					}else{\n"
  "						if(v <= -55){\n"
  "							ninf= calc_ninf(ninf60,ninf55,-60,-55)\n"
  "						}else{\n"
  "							if(v <= -50){\n"
  "								ninf= calc_ninf(ninf55,ninf50,-55,-50)\n"
  "							}else{\n"
  "								if(v <= -45){\n"
  "									ninf= calc_ninf(ninf50,ninf45,-50,-45)\n"
  "								}else{\n"
  "									if(v <= -40){\n"
  "										ninf= calc_ninf(ninf45,ninf40,-45,-40)									\n"
  "									}else{\n"
  "										if(v <= -35){\n"
  "											ninf= calc_ninf(ninf40,ninf35,-40,-35)									\n"
  "										}else{\n"
  "											if(v <= -30){\n"
  "												ninf= calc_ninf(ninf35,ninf30,-35,-30)									\n"
  "											}else{\n"
  "												if(v <= -25){\n"
  "													ninf= calc_ninf(ninf30,ninf25,-30,-25)									\n"
  "												}else{\n"
  "													if(v <= -20){\n"
  "														ninf= calc_ninf(ninf25,ninf20,-25,-20)									\n"
  "													}else{\n"
  "														if(v <= -15){\n"
  "															ninf= calc_ninf(ninf20,ninf15,-20,-15)									\n"
  "														}else{\n"
  "															if(v <= -10){\n"
  "																ninf= calc_ninf(ninf15,ninf10,-15,-10)									\n"
  "															}else{													\n"
  "																ninf= ninf10\n"
  "															}\n"
  "														}\n"
  "													}\n"
  "												}\n"
  "											}\n"
  "										}\n"
  "									}\n"
  "								}\n"
  "							}\n"
  "						}\n"
  "					}\n"
  "				}\n"
  "			}\n"
  "		}\n"
  "	}\n"
  "	\n"
  "	\n"
  "	g_vg= (B - A) * ninf * g_vgMax\n"
  "	if(g_vg < 0) {g_vg= 0}\n"
  "	i_vg= (1e-3) * g_vg * (v - e)\n"
  "} \n"
  "\n"
  "DERIVATIVE states {\n"
  "	A'= -A / tauA\n"
  "	B'= -B / tauB\n"
  "}\n"
  "\n"
  "FUNCTION calc_ninf(n_low,n_high,v_low,v_high) {\n"
  "	calc_ninf= ((v - v_low) * n_high - (v - v_high) * n_low) / (v_high - v_low)\n"
  "}\n"
  ;
#endif
