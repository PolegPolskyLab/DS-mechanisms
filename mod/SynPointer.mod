: Graded Synapse guided by presynaptic signals

NEURON {
POINT_PROCESS SynPointer
	RANGE gain,g_syn 
	GLOBAL global_gain
	RANGE e
	RANGE postX,postY,preX,preY,cellNum
	NONSPECIFIC_CURRENT i
	POINTER g
}

PARAMETER {
	g=0
	e=		-60				:reversal potential
	postX=	0				:location x
	postY=	0				:location y
	preX=	0
	preY=	0				:location of presynaptic cell
	cellNum=-1				:the presynaptic BC
	gain=	1				:gain factor
	global_gain=1			: global gain shared for all synapses

}

INITIAL {
	g= 0.0001				:100nM, set by presynaptic release
}
ASSIGNED {
	v (millivolt)
	i (nanoamp)
	g_syn
}
 
BREAKPOINT {
	if(g < 0){
		g= 0
	}
	g_syn= g * gain * global_gain
	i = (1e-3) * g_syn * (v - e)
}
 
 
