$title Finite-Horizon Ramsey Growth Model as an NLP

SET    t   time period /0*50/;

* The model logic is different in the first and last periods of the
* model, so we define some subsets here to keep track of these
* periods.

SETS
        tfirst(t)       first period
        tlast(t)        last period
;

tfirst(t) = yes$(ord(t) eq 1);
tlast(t) = yes$(ord(t) eq card(t));

PARAMETERS
*  assumed...
        g       labor growth rate                       /0.023/
        delta   capital depreciation                    /0.04/
        i0      base year investment                    /0.3/
        c0      base year consumption                   /0.27/
        b(t)       value share of capital in C-D prod fn   
        eta     inverse intertemporal elasticity        /2/

*calibrated parameters...
        k0      initial capital stock
        l0      initial labor supply
        l(t)    labor supply
        a       prod fn scale parameter
        rho     pure time preference rate
        qref(t) steady-state output index
        beta(t) discount factor
;

b(t) = 0.65 ; 


k0 = i0 / (g+delta) ; 
l0 = (1-b("0")) * (c0 + i0) ; 
l(t) = l0 * (1+g)**(t.val) ; 
a = (c0 + i0) / ((k0 ** b("0")) * (l0**(1-b("0"))) ) ; 
*!!!
rho = b("0") * (c0 + i0) / k0 - delta ; 
qref(t) = (1 + g) ** t.val ; 
*!!! 
beta(t) = (((1+g)** eta) / (1+rho))**t.val ; 

variable w "intertemporal utility" ;

positive variable
c(t) "consumption at time t", 
K(t) "capital stock at time t",
Y(t) "production"
I(t) "investment"
; 

equations objfn, production, production_consumption, capital_stock; 

objfn..
W =e= sum(t, beta(t) * C(t)**(1-eta)) / (1-eta) ;

production(t).. y(t) =e= K(t)**gamma * L(t)**delta * M(t)**(1-gamma - delta); 

production_consumption(t).. y(t) =g= c(t) + i(t) ; 

capital_stock(t)$(not tlast(t)).. K(t+1) =e= (1-delta) * K(t) + I(t) ; 

equation kfix_first;
kfix_first(t)$tfirst(t).. K(t) =e= k0 ; 

equation end_year_investment ; 

scalar end_year /0/ ; 
end_year_investment(t)$[end_year$tlast(t)].. I(t) =e= (g+delta) * K(t) ; 

model ramsey /all/; 


k.lo(t) = 0.001 ; 
c.lo(t) = 0.001 ; 


parameter report ; 
set tfix(t) /0*25/ ; 

end_year = 0 ; 

beta(t)$tlast(t) = beta(t) * (1+rho) / (rho - g) ; 



* solve the first time around
solve ramsey using nlp maximizing W ; 
report("BAU",t,"I") = i.l(t) ; 
report("BAU",t,"C") = c.l(t) ; 

c.fx(t)$tfix(t) = c.l(t) ; 
K.fx(t)$tfix(t) = K.l(t) ; 
Y.fx(t)$tfix(t) = Y.l(t) ; 
I.fx(t)$tfix(t) = I.l(t) ; 

* now have a disruption to capital productivity
* ie a capital producitivity shock

*b(t)$(not tfix(t)) = 0.5 * b(t) ; 

M.fx(t)$(not tfix(t)) = M.l("25") - 'shock' + 'stockpile' ; 

*K.up(t)$(not tfix(t)) = 0.99 * k.l(t) ; 


solve ramsey using nlp maximizing W ; 

report("shock",t,"I") = i.l(t) ; 
report("shock",t,"C") = c.l(t) ; 
execute_unload 'ramsey_data.gdx' ; 

$exit

end_year = 1 ; 
solve ramsey using nlp maximizing W ; 
report("Naive",t,"I") = i.l(t) ; 
report("Naive",t,"C") = c.l(t) ; 

solve ramsey using nlp maximizing W ; 
report("Barr-Manne",t,"I") = i.l(t) ; 
report("Barr-Manne",t,"C") = c.l(t) ; 



execute_unload 'ramsey.gdx' ; 
