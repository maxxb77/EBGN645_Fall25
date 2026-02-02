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

scalar pbar /65/ ; 
parameter qbar_s(r), qbar_d(r), cbar(r), epsilon(r), eta(r) ;

qbar_s(r) = ref_prod_cons(r,"production") ; 
qbar_d(r) = ref_prod_cons(r,"consumption") ;
cbar(r) = pbar ; 
epsilon(r) = 0.5 ;
eta(r) = 0.5 ; 

set l(r) "leader" ; 
l(r) = no ; 

positive variable P, X(r) ; 
equation mcc, zpc(r) ; 
mcc.. sum(r,X(r)) =g= sum(r,qbar_d(r) * (P / pbar)**(-epsilon(r)) ) ; 

zpc(r)$(not l(r)).. cbar(r) * (X(r) / qbar_s(r)) ** (1/eta(r)) =g= P ; 

model oil_mcp
/
mcc.p,
zpc.x
/;

p.l = pbar ; 
x.l(r) = qbar_s(r) ; 

oil_mcp.iterlim = 0 ; 

parameter rep ; 
solve oil_mcp using mcp ; 
rep("mcp","P") = p.l ; 

equation eq_mpec_obj ; 
variable mpec_obj ; 
eq_mpec_obj.. mpec_obj =e= 
    sum(r$l(r),(P-cbar(r) * (X(r) / qbar_s(r)) ** (1/eta(r))) * X(r))
;

l(r) = no ; 
l("UnitedStates") = yes ; 
model mpec_oil 
/
mcc.p,
zpc.X,
eq_mpec_obj
/;

solve mpec_oil using mpec maximizing mpec_obj ; 
rep("MPEC","P") = p.l ; 

execute_unload 'oil.gdx' ; 