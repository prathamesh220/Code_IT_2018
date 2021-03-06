param cardinality > 0  integer;
param S;
set I := {1..cardinality};

param p{i in I}=(1/i)^S/(sum {j in I}(1/j)^S);
param formatting symbolic  in I;
param formatting_2 symbolic in I;

param entropy=sum{i in I}(-p[i]*(log(p[i])/log(2)));
param entropy_2=(sum{i in I}(p[i]*(-log(p[i])/log(2))^2));
param vare= entropy_2-entropy^2;
var q{i in I}>=1/10^12, :=1/cardinality;
param pracentropy=(sum{i in I}(p[i]*(ceil(-log(p[i])/log(2)))));
param pracentropy_2=(sum{i in I}(p[i]*(ceil(-log(p[i])/log(2)))^2));

var lambda>=0;
var tilteddistribution{i in I} = (((sum {j in I}(lambda*q[j]^(1/2)*p[j]^(1/2)))+1-(lambda^2/2))/(lambda*q[i]^(1/2)*p[i]^(1/2)+p[i]-(lambda^2/2)*p[i]))^(-1) ;
var tiltedentropy=(sum {i in I}(-tilteddistribution[i]*log(tilteddistribution[i])))/log(2);
var codelength= (sum {i in I}(p[i]*(-log(tilteddistribution[i])/log(2))));
var codelength_2= (sum {i in I}(p[i]*((log(tilteddistribution[i])/log(2))^2)));
var praclengths{i in I} = ceil(log((((sum {j in I}(lambda*q[j]^(1/2)*p[j]^(1/2)))+1-(lambda^2/2))/(lambda*q[i]^(1/2)*p[i]^(1/2)+p[i]-(lambda^2/2)*p[i])))/log(2));
var praccodelength= (sum {i in I}( p[i]*praclengths[i]));
var praccodelength_2=(sum {i in I}( p[i]*(praclengths[i]^2)));
var pracavgage= praccodelength+ praccodelength_2/(2*praccodelength) - 0.5;
param pshannonage=pracentropy+pracentropy_2/(2*(pracentropy))-0.5;
param shannonage= entropy+entropy_2/(2*(entropy))-0.5;
param sum_p=sum{i in I}p[i];
var sum_td=sum{i in I}tilteddistribution[i];
var check_optimal=codelength+ codelength_2/(2*(codelength))-0.5 ;

maximize avgage:
     sum {i in I} ((lambda*q[i]^(1/2)*p[i]^(1/2)+p[i]-(lambda^2/2)*p[i])*log(((sum {j in I}(lambda*q[j]^(1/2)*p[j]^(1/2)))+1-(lambda^2/2))/(lambda*q[i]^(1/2)*p[i]^(1/2)+p[i]-(lambda^2/2)*p[i])))/log(2)-0.5;

subject to Fill {i in I}:
    lambda*q[i]^(1/2)*p[i]^(1/2)+p[i]-(lambda^2/2)*p[i] >=1/10^100;
	
subject to distribution :
    sum {i in I} q[i]=1.0;	
	
subject to lanbda_constraint:
         lambda <= 10^6;
