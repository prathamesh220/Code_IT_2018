
set symbol := 1..8;

param p{symbol};
var q{symbol}>=0;
var lambda;


maximize crossentropy:
     sum {i in symbol} lambda*q[i]^(1/2)*p[i]^(1/2)+p[i]-(lambda^2/2)*p[i];

subject to Fill {i in symbol}:
    lambda*q[i]^(1/2)*p[i]^(1/2)+p[i]-(lambda^2/2)*p[i] >=1/10^6;
