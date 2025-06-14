
NEURON {
	SUFFIX caGA
	USEION ca READ eca WRITE ica
	:GLOBAL mN
	RANGE m, h, mN, gMax,gca
}

PARAMETER {
	gMax = 0.1   	(pS/um2)
	mN= 1	
	v 		(mV)
}

ASSIGNED {
	ica 				(mA/cm2)
	gca				(pS/um2)
	eca 				(mV)
}

STATE { m h }

BREAKPOINT {
    SOLVE states METHOD cnexp
    gca = gMax * (m)^mN * h
	ica = (1e-4) * gca * (v - eca)
} 

DERIVATIVE states {
    m' =  (minf(v)-m)/(mtau(v))
    h' =  (hinf(v)-h)/(htau(v))
}

FUNCTION_TABLE mtau(v(mV)) (ms)
FUNCTION_TABLE htau(v(mV)) (ms)
FUNCTION_TABLE minf(v(mV))
FUNCTION_TABLE hinf(v(mV))