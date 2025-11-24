scalar
alpha /0.75/, 
pb /2/, 
pc /4/, 
I /16/ ; 

equation eq_obj, income_limit ; 
variable U "utility", b "butter", c "cheese" ; 

eq_obj.. U =e= (b**alpha) * (c**(1-alpha)) ; 
income_limit I =g= pb * b + pc * c ; 

model butter_cheese /all/ ;

solve butter_cheese using nlp maximizing U ; 