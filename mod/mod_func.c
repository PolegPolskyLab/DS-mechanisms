#include <stdio.h>
extern void _SynPointer_reg();
ORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _SynVec_reg();
extern void _caGA_reg();
extern void _cadiff_reg();
extern void _calRGC_reg();
extern void _canrgc_reg();
extern void _glutamate_reg();
extern void _ih_reg();
extern void _kGA_reg();
extern void _ka_reg();
extern void _kap_reg();
extern void _kca_reg();
extern void _kdf_reg();
extern void _kdr2_reg();
extern void _kfast_reg();
extern void _km_reg();
extern void _kslow_reg();
extern void _kv_reg();
extern void _mGluR_reg();
extern void _synGA_reg();

fprintf(stderr," SynPointer.mod");
dio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," SynVec.mod");
fprintf(stderr," caGA.mod");
fprintf(stderr," cadiff.mod");
fprintf(stderr," calRGC.mod");
fprintf(stderr," canrgc.mod");
fprintf(stderr," glutamate.mod");
fprintf(stderr," ih.mod");
fprintf(stderr," kGA.mod");
fprintf(stderr," ka.mod");
fprintf(stderr," kap.mod");
fprintf(stderr," kca.mod");
fprintf(stderr," kdf.mod");
fprintf(stderr," kdr2.mod");
fprintf(stderr," kfast.mod");
fprintf(stderr," km.mod");
fprintf(stderr," kslow.mod");
fprintf(stderr," kv.mod");
fprintf(stderr," mGluR.mod");
fprintf(stderr," synGA.mod");
fprintf(stderr, "\n");
    }
_SynPointer_reg();
_SynVec_reg();
_caGA_reg();
_cadiff_reg();
_calRGC_reg();
_canrgc_reg();
_glutamate_reg();
_ih_reg();
_kGA_reg();
_ka_reg();
_kap_reg();
_kca_reg();
_kdf_reg();
_kdr2_reg();
_kfast_reg();
_km_reg();
_kslow_reg();
_kv_reg();
_mGluR_reg();
_synGA_reg();
}
