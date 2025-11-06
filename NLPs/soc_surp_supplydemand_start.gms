
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

positive variables Qd, Qs ; 
variable Z "objective function" ; 

equation
eq_obj, eq_supplydemand ;

scalar tax_d /0/ ; 

eq_obj.. Z =e= a * Qd - b * Qd * Qd / 2
             - c * Qs - d * Qs * Qs / 2 - tax_d * Qd; 
* Qd * Qd equivalent to..
* power(qd, 2) ; 


eq_supplydemand.. Qd =e= Qs ; 

model surplussing /all / ;

parameter rep ; 
solve surplussing using nlp maximizing Z ; 

rep("bau","q") = qd.l ; 
rep("bau","p") = eq_supplydemand.m ; 

* increase a
a = 1.2 * a ; 
solve surplussing using nlp maximizing Z ; 
rep("cf","q") = qd.l ; 
rep("cf","p") = eq_supplydemand.m ; 

*reset out a
a = a0/b0; 
* tax that bizness
tax_d = 0.25 ; 
solve surplussing using nlp maximizing Z ; 
rep("tax","q") = qd.l ; 
rep("tax","p") = eq_supplydemand.m ; 
rep("tax","true_price") = a0 - b0 * qd.l ; 




execute_unload 'surplus.gdx' ; 


