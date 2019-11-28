
function flag=RobustNondegeneracy(G)
   G=[-1  1  0  0  0  0  0 0 ;
      0 -1  1  0  0  0  0 0; 
      0  0 -1  1  0  0  0 0 ;
      1  0  0 -1  1  0  0 -1;
      0  0  0  0 -1  1  0 0 ;
      0  0  0  0  0 -1  1 0;
      0  0  0  0  0  0 -1 1 ; 
      1 -1  0  0  0  0  0 0 ;
      0  0  0  0  1 -1  0 0];

[n,r]=size(G);
oldparam = sympref('HeavisideAtOrigin',0);

rho=heaviside(-G)'.*ones(r,n);%sym('r',[ r n]);
rho=G*rho;

V=null(G','r')';
d=size(V,1);

I=eye(n);

T=[I(:,1:n-d)'; V];

J=T*rho*inv(T);

Jr=J(1:n-d,1:n-d);

flag=n-d-rank(Jr);






