
parameter a0, b0, c0, d0; 

b0 = 2 ; 
d0 = 2 ;

scalar pbar /1/, qbar /1/ ; 

a0 = qbar + b0 * pbar ; 
c0 = qbar - d0 * pbar ; 

parameter a, b, c, d ; 

a = a0/b0; 
b = 1/b0 ; 
c = -(c0/d0) ;
d = 1/d0 ; 


positive variables 
Qd "quantity demanded"
Qs "quantity supplied"
;

variable W "social welfare" ; 

equation objfn, sup_dem ; 

objfn.. W =e= a * Qd - b * Qd * Qd / 2 
            - (c * Qs + d * Qs * Qs / 2) ; 

sup_dem.. Qs =g= Qd ; 

model welfy /all/ ; 

Qd.lo = 0.01 ; 
Qs.lo = 0.01 ; 

parameter rep ; 
solve welfy using nlp maximizing W ; 

rep('ref','Q') = qd.l ; 
rep('ref','P') = sup_dem.m ; 

a=1.2 * a ; 
solve welfy using nlp maximizing W ; 

rep('cf','Q') = qd.l ; 
rep('cf','P') = sup_dem.m ; 


execute_unload 'alldata.gdx' ; 
