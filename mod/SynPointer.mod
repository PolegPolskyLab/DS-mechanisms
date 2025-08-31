: Graded Synapse guided by presynaptic signals

NEURON {
POINT_PROCESS SynPointer
	RANGE gain,g_syn 
	RANGE e, baseline
	NONSPECIFIC_CURRENT i
	POINTER pre
}

PARAMETER {
	e=			-60				:reversal potential
	gain=		1				:gain factor
	baseline=	0.0001			: make -60 for voltage
	pre=		0				:presynsptic volatage/calcium
}

INITIAL {
	pre= baseline				:100nM, set by presynaptic release
}
ASSIGNED {
	v (millivolt)
	i (nanoamp)
	g_syn
}
 
BREAKPOINT {
	g_syn= (pre-baseline) * gain
	if(g_syn < 0){
		g_syn= 0
	}
	
	i = (1e-3) * g_syn * (v - e)
}
 
 
