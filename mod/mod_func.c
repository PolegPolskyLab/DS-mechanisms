#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _BK_reg();
extern void _SynPointer_reg();
extern void _SynVec_reg();
extern void _caGA_reg();
extern void _ca_lin_reg();
extern void _cadiff_reg();
extern void _calRGC_reg();
extern void _calRGCfix_reg();
extern void _canrgc_reg();
extern void _glutamate_reg();
extern void _ih_reg();
extern void _kGA_reg();
extern void _kap_reg();
extern void _kca_reg();
extern void _km_reg();
extern void _kslow_reg();
extern void _kv_reg();
extern void _nav12_reg();
extern void _nav16_reg();
extern void _spike_reg();
extern void _synGA_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," BK.mod");
fprintf(stderr," SynPointer.mod");
fprintf(stderr," SynVec.mod");
fprintf(stderr," caGA.mod");
fprintf(stderr," ca_lin.mod");
fprintf(stderr," cadiff.mod");
fprintf(stderr," calRGC.mod");
fprintf(stderr," calRGCfix.mod");
fprintf(stderr," canrgc.mod");
fprintf(stderr," glutamate.mod");
fprintf(stderr," ih.mod");
fprintf(stderr," kGA.mod");
fprintf(stderr," kap.mod");
fprintf(stderr," kca.mod");
fprintf(stderr," km.mod");
fprintf(stderr," kslow.mod");
fprintf(stderr," kv.mod");
fprintf(stderr," nav12.mod");
fprintf(stderr," nav16.mod");
fprintf(stderr," spike.mod");
fprintf(stderr," synGA.mod");
fprintf(stderr, "\n");
    }
_BK_reg();
_SynPointer_reg();
_SynVec_reg();
_caGA_reg();
_ca_lin_reg();
_cadiff_reg();
_calRGC_reg();
_calRGCfix_reg();
_canrgc_reg();
_glutamate_reg();
_ih_reg();
_kGA_reg();
_kap_reg();
_kca_reg();
_km_reg();
_kslow_reg();
_kv_reg();
_nav12_reg();
_nav16_reg();
_spike_reg();
_synGA_reg();
}
