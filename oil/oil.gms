set r
/
UnitedStates
SaudiArabia
Russia
Canada
China
Iraq
Brazil
UAE
Iran
Kuwait
India
Japan
SouthKorea
Germany
ROW
/;

table ref_prod_cons(r,*)
$ondelim
$include ref_prod_cons.csv 
$offdelim
;

scalar pbar /58.64/ ; 
parameter qbar_s(r), qbar_d(r), cbar(r), epsilon(r), eta(r) ;

qbar_s(r) = ref_prod_cons(r,"production") ; 
qbar_d(r) = ref_prod_cons(r,"consumption") ;
cbar(r) = pbar ; 
epsilon(r) = 0.5 ;
eta(r) = 0.5 ; 

cbar(r) = pbar ; 

positive variables P, X(r), lambda, lambda_us ; 
equations mcc, zpc(r), russia_limit ;


scalar sw_russialimit /0/ ; 

set l(r) "leader" ; 
l(r) = no ; 

mcc.. sum(r,X(r)) =g= sum(r, qbar_d(r) * (P/pbar)**(-epsilon(r))) ; 
zpc(r)$[not l(r)].. cbar(r) * (X(r)/qbar_s(r)) ** (1/eta(r)) 
        + lambda$[sameas(r,"russia")$sw_russialimit] 
        + lambda_us$[sameas(r,"unitedstates")$sw_russialimit] 
        =g= P ; 
russia_limit$sw_russialimit.. 0.5 * qbar_s("russia") =g= X("russia") ; 

equation us_limit ; 
* us can only icnrease by 1% of reference production
us_limit$sw_russialimit.. 1.01 * qbar_s("unitedstates") =g= X("unitedstates") ; 



model oil_mcp 
/
mcc.p
zpc.x
russia_limit.lambda
us_limit.lambda_us
/; 

* setting level values of variables equal to their reference quantities/prices
p.l = pbar ; 
x.l(r) = qbar_s(r) ; 

oil_mcp.iterlim = 0 ; 

parameter rep_price, rep_quantity, rep_consumption ; 

* replicate the benchmark
solve oil_mcp using mcp ;
rep_price("bau") = p.l ; 
rep_quantity("bau",r) = x.l(r) ; 
rep_consumption("bau",r) = qbar_d(r) * (P.l/pbar) **(-epsilon(r)) ; 

* release the hounds
oil_mcp.iterlim = 1e5 ;
* limit russia
sw_russialimit = 1 ; 
* replicate the benchmark
solve oil_mcp using mcp ;
rep_price("russialimit") = p.l ; 
rep_quantity("russia_limit",r) = x.l(r) ; 
rep_consumption("russia_limit",r) = qbar_d(r) * (P.l/pbar) **(-epsilon(r)) ; 

equation eq_mpec_obj ; 
variable mpec_obj ; 

eq_mpec_obj.. mpec_obj =e= sum(r$l(r), (P - cbar(r) * (X(r)/qbar_s(r)) ** (1/eta(r))) * X(r) ) ; 

model oil_mpec
/
eq_mpec_obj
mcc.p
zpc.x
russia_limit.lambda
us_limit.lambda_us
/; 

sw_russialimit = 0 ; 

alias(r,rr) ; 

loop(rr,
l(r) = no ; 
l(rr) = yes ; 
solve oil_mpec using mpec maximizing mpec_obj ; 
x.fx(rr) = x.l(rr) ; 
) ; 

rep_price("pimax") = p.l ; 
rep_quantity("pimax",r) = x.l(r) ; 
rep_consumption("pimax",r) = qbar_d(r) * (P.l/pbar) **(-epsilon(r)) ; 


execute_unload 'oil.gdx' ; 