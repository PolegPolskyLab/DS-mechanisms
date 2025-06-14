

NEURON {
	SUFFIX cadiff
	USEION ca READ ica, cai WRITE cai
	GLOBAL taur,cainf,factor
}

CONSTANT {
	FARADAY = 96489				
}

PARAMETER {
	taur	= 50				
	cainf	= 0.0001
	factor=100000	
}
ASSIGNED {
	ica
}
STATE {
	cai	 
}

INITIAL {
	cai = cainf
}

BREAKPOINT {
	SOLVE state METHOD euler
}

DERIVATIVE state {
	cai'=-(factor)*ica/(2 * FARADAY )-(cai-cainf)/taur
}
