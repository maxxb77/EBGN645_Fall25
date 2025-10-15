set y "year" /1*100/ ; 

scalar p  timber price / 1000 /;
scalar r  discount rate / 0.05 /;
scalar c  replanting cost / 100 /;
scalar f growth_rate / 0.1 /;  
scalar init_stock /100/ ; 

variable Z present_value;
positive variable H(y) "per-period harvest", stock(y) "standing stock";

* Equations
equation objective, stock_tracking, stock_limit ; 

set ymax(y) ; 
alias(y,yy) ; 
ymax(y)$[y.val = smax(yy,yy.val)] = yes ; 

* Model
objective.. Z =e= sum(y,(1/(1+r)**y.val) * (
        H(y) * (p - c * H(y)/2)) + ((p-c*stock(y)/20) * stock(y))$ymax(y)
        );

stock_tracking(y).. stock(y) =e= (stock(y-1) - H(y-1))*(1+f)
                                 + init_stock$(y.val=1) ; 

stock_limit(y).. stock(y) =g= H(y) ; 

stock.lo(y) = init_stock / 50; 

model managed /all/ ; 

parameter manage_out ; 

set rset "interest rates to run" /5, 7, 10/ ; 
set fset "growth rates to run"   /8, 10, 12/ ;

loop(rset,
    loop(fset,
        r = rset.val / 100 ; 
        f = fset.val / 100 ; 
        solve managed maximizing Z using nlp ; 
        manage_out(rset,fset,"stock",y) = stock.l(y) ; 
        manage_out(rset,fset,"harvest",y) = h.l(y) ; 
    )
)


execute_unload 'managed.gdx' ; 



* Model Type