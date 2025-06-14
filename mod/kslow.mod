
NEURON {
    
	SUFFIX kSlow
	USEION k READ ek WRITE ik
	RANGE n, h, gk, gkbar
	RANGE ninf, hinf
	RANGE ntau, htau
	GLOBAL nmin,nslope,hmin,hslope
	GLOBAL shift,tau_shift
}

PARAMETER {
	gkbar = 0.1   	(pS/um2)	: 0.12 mho/cm2
	ek=-80
	nmin=-40
	nslope=.013
	hmin=-80
	hslope=0.02	
	v 		(mV)
	
	shift=0
	tau_shift=1
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
	PI	= (pi) (1)
} 

ASSIGNED {
	ik 		(mA/cm2)
	gk		(pS/um2)
	ninf 		
	hinf	
	ntau
	htau
}
 

STATE { n h }

INITIAL {
	rates(v+shift)
	n = ninf
	h = hinf
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    gk = gkbar*(n)*h
	ik = (1e-4) * gk * (v - ek)
} 

LOCAL nexp, hexp


DERIVATIVE states {
    rates(v+shift)      
    n' =  (ninf-n)/(ntau*tau_shift)
    h' =  (hinf-h)/(htau*tau_shift)
}

PROCEDURE rates(vm) {  
:based on korngreen and sakmann 2000

	
	:Activation
	ninf = nslope*(vm-nmin)
	if(ninf<0){
		ninf=0
	}
	if(ninf>1){
		ninf=1
	}
	if(vm<-50){
		ntau=50-(-50-vm)
	}else{
		ntau=50+(-50-vm)/1.5
	}	
	if(ntau<5){ntau=5}
	
	
	:Deactivation
	hinf = 1-hslope*(vm-hmin)
	if(hinf<0){
		hinf=0
	}
	if(hinf>1){
		hinf=1
	}	
	htau=300
	if(vm<0){
		htau=-vm*10+300
	}	
}

