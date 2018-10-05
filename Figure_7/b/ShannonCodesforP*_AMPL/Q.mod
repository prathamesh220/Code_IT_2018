param cardinality > 0  integer;
param S;
set I := {1..cardinality};
set I_1:= {1..(cardinality-1)};
param p{i in I};

param formatting symbolic  in I;
param formatting_2 symbolic in I;
param mu;

param entropy=sum{i in I}(-p[i]*(log(p[i])/log(2)));
param entropy_2=sum{i in I}(p[i]*(-log(p[i])/log(2))^2);
param pracentropy=(sum{i in I}(p[i]*(ceil(-log(p[i])/log(2)))));
param pracentropy_2=(sum{i in I}(p[i]*(ceil(-log(p[i])/log(2)))^2));
param vare=entropy_2 -entropy^2;
var q{i in I}>=1/10^6, :=1/cardinality;

var lambda>=0;
var tilteddistribution{i in I} = (((sum {j in I}(lambda*mu^(1/2)*q[j]^(1/2)*p[j]^(1/2)))+1+(lambda^2*mu/2))/(lambda*mu^(1/2)*q[i]^(1/2)*p[i]^(1/2)+p[i]+(lambda^2*mu/2)*p[i]))^(-1) ;
var tiltedentropy=(sum {i in I}(-tilteddistribution[i]*log(tilteddistribution[i])))/log(2);
var codelength= (sum {i in I}(-p[i]*(log(tilteddistribution[i])/log(2))));
var codelength_2= (sum {i in I}(p[i]*((log(tilteddistribution[i])/log(2))^2)));
var praclengths{i in I} = ceil(-log(tilteddistribution[i])/log(2));
var praccodelength= (sum {i in I}( p[i]*praclengths[i]));
var praccodelength_2=(sum {i in I}( p[i]*(praclengths[i]^2)));
var pracavgage= praccodelength+ mu*praccodelength_2/(2*(1-mu*praccodelength)) ;
param pshannonage= pracentropy+mu*pracentropy_2/(2*(1-mu*pracentropy));
param shannonage= entropy+mu*entropy_2/(2*(1-mu*entropy));
param sum_p=sum{i in I}p[i];
var sum_td=sum{i in I}tilteddistribution[i];
var check_optimal=codelength+ mu*codelength_2/(2*(1-mu*codelength)) ;
maximize avgage:
     sum {i in I} ((lambda*mu^(1/2)*q[i]^(1/2)*p[i]^(1/2)+p[i]+(lambda^2*mu/2)*p[i])*log(tilteddistribution[i]^(-1)))/log(2)-(lambda^2/2);


	
subject to distribution :
    sum {i in I} q[i]=1.0;	
	
