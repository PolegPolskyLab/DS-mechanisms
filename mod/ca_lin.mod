
: "Linear" calcium channel
NEURON {
	SUFFIX caLinear
	USEION ca WRITE ica
	RANGE  gbar, ica
	
}

PARAMETER {
	gbar = .000001   	(mho/cm2)	
}

ASSIGNED {
	ica 		(mA/cm2)
	v
}
 
BREAKPOINT {
	if(v>-60){
		ica = -gbar*(v + 60)
	}else{
		ica = 0
	}
} 
