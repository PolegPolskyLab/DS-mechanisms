
NEURON {
	SUFFIX kGA
	USEION k READ ek WRITE ik
	RANGE n, h, nN, gMax,gk
}

PARAMETER {
	gMax = 0.1   	(pS/um2)
	nN= 1	
	v 		(mV)
}

ASSIGNED {
	ik 				(mA/cm2)
	gk				(pS/um2)
	ek 				(mV)
}

STATE { n h }

BREAKPOINT {
    SOLVE states METHOD cnexp
    gk = gMax * (n)^nN * h
	ik = (1e-4) * gk * (v - ek)
} 

DERIVATIVE states {
    n' =  (ninf(v)-n)/(ntau(v))
    h' =  (hinf(v)-h)/(htau(v))
}

FUNCTION_TABLE ntau(v(mV)) (ms)
FUNCTION_TABLE htau(v(mV)) (ms)
FUNCTION_TABLE ninf(v(mV))
FUNCTION_TABLE hinf(v(mV))