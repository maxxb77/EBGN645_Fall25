set i /hamburgers, hotdogs, frenchfries/ ; 

parameter p(i)

/
hamburgers 3, 
hotdogs 2, 
frenchfries 1
/;

scalar hbar 'hours in a week' /40/ ; 


parameter h(i)
/
hamburgers 2, 
hotdogs 1, 
frenchfries 0.5
/;

positive variable x(i) production ; 

variable z 'objective function' ; 

equation objfn, timelimit ; 

objfn.. Z =e= sum(i, p(i) * x(i)) ; 
timelimit.. hbar =g= sum(i,h(i) * x(i)) ; 

model ripley /all/ ; 

solve ripley using lp maximizing z ; 

equation eq_dual_obj, eq_zpc_x ; 
positive variable lambda_h ; 
variable dual_obj ; 

eq_dual_obj.. dual_obj =e= lambda_h * hbar ; 
eq_zpc_x(i).. lambda_h * h(i) =g= p(i) ; 

model ripley_dual /eq_dual_obj, eq_zpc_x/ ; 

solve ripley_dual minimizing dual_obj using lp ; 

model ripley_mcp 
/
timelimit.lambda_h,
eq_zpc_x.x
/;

ripley_mcp.iterlim = 0 ; 
solve ripley_mcp using mcp ; 
