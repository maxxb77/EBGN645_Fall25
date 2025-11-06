$Title Simple Hourly Dispatch LP
* Maxwell Brown

option seed=7

$if not set sw_season $setglobal sw_season "day"
$if not set sw_co_reduction $setglobal sw_co_reduction 10
$if not set sw_wy_reduction $setglobal sw_wy_reduction 7.5

* enable elastic demand in the counterfactuals?
$if not set sw_elas $setglobal sw_elas 0



*=======================
*  Begin sets
*=======================

set h "hour" /h1*h24/;
set a 'temp hour' /1*24/ ; 


alias(h,hh) ; 

set h_prev(h,hh) ; 

h_prev(h,hh)$[ord(h) = ord(hh) + 1] = yes ; 

set dayhours(h) /h8*h20/; 

set k "season" /day, summer/;

set s "states"
/
  CO "Colorado",
  WY "Wyoming"
/;

*State-specific policy switches
Set SwTPS(s)  "Switch to enable or disable a TPS for state s";
Set SwCAP(s)  "Switch to enable or disable a carbon cap for state s";

*disable these for now
swtps(s) = no;
swcap(s) = no;

set f "fuels, generation technology"
        /
        bit "bituminous coal"
        dfo "distillate fuel oil"
        ng  "natural gas"
        sub "subbituminous coal"
        sun "solar"
        wat "water"
        wnd "wind"
        /;

set thermal(f) "technologies subject to ramping constraints"
/ bit, dfo, ng, sub /; 

set pc  "plant characteristics" /cap, hr, onm/;
set pid "plant id" /1*328/;
set genfeas(s,f,pid,h) "generation feasibility set, determines when a plant can generate";

*=======================
*  End sets
*=======================

*=======================
*  Begin data
*=======================

* some rudimentary data...
parameter capfac(f) "capacity availability by fuel type"
         /
         bit 0.8,
         dfo 0.65,
         ng  0.85,
         sub 0.8,
         sun 0.5,
         wat 0.55,
         wnd 0.45
         /;
* offhand guesses 

parameter cf(f,h) 'capacity factor by technology and hour' ; 

* start by broadcasting values
cf(f,h) = capfac(f) ; 

* remove those not possible
cf("sun",h)$(not dayhours(h)) = 0 ; 


parameter fcost(f) "costs by fuel ($ / MMBTU)"
         /
         bit  2.9,
         dfo  11.1,
         ng   2.97,
         sub  2.9
         /;
*estimates from EIA, taken from short term energy outlook

parameter emit(f) "lbs co2 per mmbtu for fuel burning"
         /
         bit  202,
         dfo  161,
         ng   117,
         sub  209
         /;

*Policy Parameters
*Unless otherwise defined, the Psi parameters are relative to the base case
*note that we cannot be too stringent w/o elastic demand, safety valve credits, or some capacity expansion capabilities
Parameter  red(s) "Percentage reduction of average CI or total carbon" 
                /CO 5, WY 5/;

red(s) = red(s) / 100;


parameter 
  psi_state(s), 
  psi_trade, 
  phi_state(s), 
  phi_trade ; 

* set all to zero to start
psi_state(s) = 0 ;  
psi_trade = 0 ; 
phi_state(s) = 0 ; 
phi_trade = 0 ; 


table d_in(h,k) 
$include refdem.inc
; 

parameter d(h) "demand by hour (MW)" ; 

d(h) = d_in(h,"%Sw_Season%") ; 

table plantdata(s,f,pid,pc)
$include plantdata.inc
;

set v(s,f,pid,h) "allowable combinations for technologies over all indices"; 
*populate combinations..
* - if you have a capacity factor (specific to solar)
* - if you have capacity
v(s,f,pid,h)$[cf(f,h)$plantdata(s,f,pid,"cap")] = yes ; 


*=======================
*  End data
*=======================


*=======================
*  Begin Model
*=======================

equation 
eq_cost "objective function",
eq_capcon(s,f,pid,h) "generation cannot exceed capacity",
eq_demcon(h) "generation must equal demand"
eq_carboncap(s) "(optional) total emissions cannot exceed cap",
eq_ratestandard(s) "(optional) emissions rate cannot exceed emissions standard"
;

positive variable X(s,f,pid,h) "generation in MWH"
*                  dem_b(h,b) "demand by bin"
;

variable Z "objective function value ($s)";

parameter c_v(s,f,pid) "variable cost", hr(s,f,pid) "heat rate", cap(s,f,pid), e(s,f,pid) "emissions rate" ; 

c_v(s,f,pid) = plantdata(s,f,pid,"onm") ; 
hr(s,f,pid) = plantdata(s,f,pid,"hr") ; 
cap(s,f,pid) = plantdata(s,f,pid,"cap") ; 
e(s,f,pid) = hr(s,f,pid) * emit(f) ; 


parameter psi(s), phi(s) ; 
psi(s) = 0 ; 
phi(s) = 0 ; 

scalar sw_cap /0/, sw_tps /0/ ; 

positive variable slack_state ; 
eq_cost.. Z =e= 
    sum((s,f,pid,h)$v(s,f,pid,h),X(s,f,pid,h)*(c_v(s,f,pid) + fcost(f) * hr(s,f,pid))) ; 

eq_capcon(s,f,pid,h)$v(s,f,pid,h).. cf(f,h) * cap(s,f,pid) =g= X(s,f,pid,h) ; 

eq_demcon(h).. sum((s,f,pid)$v(s,f,pid,h), X(s,f,pid,h)) =g= d(h) ; 


eq_carboncap(s)$[psi(s)$sw_cap].. 
  psi(s) =g= sum((f,pid,h)$v(s,f,pid,h),e(s,f,pid) * X(s,f,pid,h)); 

eq_ratestandard(s)$[phi(s)$sw_tps].. 
  phi(s) * sum((f,pid,h)$v(s,f,pid,h),X(s,f,pid,h)) 
           =g= sum((f,pid,h)$v(s,f,pid,h),e(s,f,pid) * X(s,f,pid,h)) ; 

model cowy /all/ ; 

parameter rep;
solve cowy using lp minimizing z ; 
rep("BAU",s,f,pid,h)$v(s,f,pid,h) = x.l(s,f,pid,h) ; 

parameter ref_emit(s) "reference emissions"; 
ref_emit(s) = sum((f,pid,h),e(s,f,pid) * X.l(s,f,pid,h)) ; 

psi(s) = (1-red(s)) * ref_emit(s) ; 
phi(s) = psi(s) / sum((f,pid,h),X.l(s,f,pid,h)) ;

sw_cap = 1 ; 
solve cowy using lp minimizing z ; 
rep("CAP",s,f,pid,h)$v(s,f,pid,h) = x.l(s,f,pid,h) ; 

sw_cap = 0 ; 
sw_tps = 1 ; 

solve cowy using lp minimizing z ; 
rep("TPS",s,f,pid,h)$v(s,f,pid,h) = x.l(s,f,pid,h) ; 

execute_unload 'cowy.gdx' ; 
