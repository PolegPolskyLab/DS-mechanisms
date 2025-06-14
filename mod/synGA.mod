
NEURON {
	POINT_PROCESS synGA
	:SUFFIX synGA
	NONSPECIFIC_CURRENT i_lin, i_vg
	RANGE g_linMax, g_vg, g_vgMax, g_vg, g, g_lin
	RANGE ninf 
	GLOBAL ninf80, ninf70,ninf60, ninf50, ninf40, ninf30, ninf20, ninf10
	GLOBAL ninf75, ninf65, ninf55, ninf45, ninf35, ninf25, ninf15
	GLOBAL tauA, tauB
	GLOBAL e
}

PARAMETER {
		:linear
	g_linMax= 0.1	(nS)
	g= 0		(nS)
		:voltage gated
	g_vgMax= 0.1	(nS)
	tauA= 10		(ms)
	tauB= 100		(ms)
	ninf80= 0
	ninf70= 0
	ninf60= 0
	ninf50= 0
	ninf40= 0
	ninf30= 0
	ninf20= 0
	ninf10= 0
	ninf75= 0
	ninf65= 0
	ninf55= 0
	ninf45= 0
	ninf35= 0
	ninf25= 0
	ninf15= 0

		:params
	v 				(mV)
	e= 0			(mV)
}

ASSIGNED {
	i_lin 			(nA)		
	i_vg			(nA)
	g_lin
	g_vg	
	ninf
	
}

STATE { 	
	A 		(nS)
	B 		(nS) 
}

BREAKPOINT {
		:linear
	g_lin= g * g_linMax
	
	i_lin= (1e-3) * g_lin * (v - e)

		:voltage gated
	:state_discontinuity(A, A + g_vgMax)
	:state_discontinuity(B, B + g_vgMax)
	A= A + g
	B= B + g
    SOLVE states METHOD cnexp
	if(v <= -80){
		ninf= ninf80
	}else{
		if(v <= -75){
			ninf= calc_ninf(ninf80,ninf75,-80,-75)
		}else{
			if(v <= -70){
				ninf= calc_ninf(ninf75,ninf70,-75,-70)
			}else{
				if(v <= -65){
					ninf= calc_ninf(ninf70,ninf65,-70,-65)
				}else{				
					if(v <= -60){
						ninf= calc_ninf(ninf65,ninf60,-65,-60)
					}else{
						if(v <= -55){
							ninf= calc_ninf(ninf60,ninf55,-60,-55)
						}else{
							if(v <= -50){
								ninf= calc_ninf(ninf55,ninf50,-55,-50)
							}else{
								if(v <= -45){
									ninf= calc_ninf(ninf50,ninf45,-50,-45)
								}else{
									if(v <= -40){
										ninf= calc_ninf(ninf45,ninf40,-45,-40)									
									}else{
										if(v <= -35){
											ninf= calc_ninf(ninf40,ninf35,-40,-35)									
										}else{
											if(v <= -30){
												ninf= calc_ninf(ninf35,ninf30,-35,-30)									
											}else{
												if(v <= -25){
													ninf= calc_ninf(ninf30,ninf25,-30,-25)									
												}else{
													if(v <= -20){
														ninf= calc_ninf(ninf25,ninf20,-25,-20)									
													}else{
														if(v <= -15){
															ninf= calc_ninf(ninf20,ninf15,-20,-15)									
														}else{
															if(v <= -10){
																ninf= calc_ninf(ninf15,ninf10,-15,-10)									
															}else{													
																ninf= ninf10
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	
	g_vg= (B - A) * ninf * g_vgMax
	if(g_vg < 0) {g_vg= 0}
	i_vg= (1e-3) * g_vg * (v - e)
} 

DERIVATIVE states {
	A'= -A / tauA
	B'= -B / tauB
}

FUNCTION calc_ninf(n_low,n_high,v_low,v_high) {
	calc_ninf= ((v - v_low) * n_high - (v - v_high) * n_low) / (v_high - v_low)
}