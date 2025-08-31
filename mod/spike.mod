TITLE HH channel


NEURON {
	SUFFIX spike
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	RANGE gnabar, gkbar
	GLOBAL taus,taun,taum,tauh,tausb
	GLOBAL tausv,tausd,mN,nN
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
	v (mV)


	gnabar=.20 (mho/cm2)
	gkbar=.12 (mho/cm2)
	ena = 40 (mV)
	ek = -80 (mV)
	taum=0.05
	tauh=0.5
	taus=50
	tausv=30
	tausd=1
	taun=1
	mN=3
	nN=3
	tausb=0.5
}
STATE {
	m h n s
}
ASSIGNED {
	ina (mA/cm2)
	ik (mA/cm2)

}

BREAKPOINT {
	SOLVE states
	ina = gnabar*h*s*(v - ena)*m^mN
	ik = gkbar*(v - ek)*n^nN
}

PROCEDURE states() {	: exact when v held constant
	LOCAL sigmas
	sigmas=1/(1+exp((v+tausv)/tausd))
	m = m + (1 - exp(-dt/taum))*(1 / (1 + exp((v + 40)/(-3)))  - m)
	h = h + (1 - exp(-dt/tauh))*(1 / (1 + exp((v + 45)/3))  - h)
	s = s + (1 - exp(-dt/(taus*sigmas+tausb)))*(1 / (1 + exp((v + 44)/3))  - s)
	n = n + (1 - exp(-dt/taun))*(1 / (1 + exp((v + 40)/(-3)))  - n)
	VERBATIM
	return 0;
	ENDVERBATIM
}

