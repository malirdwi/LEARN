function flag=RobustNondegeneracy(G)
checkToolbox

[n,r]=size(G);
oldparam = sympref('HeavisideAtOrigin',0);

rho0=heaviside(-G)'.*ones(r,n);%sym('r',[ r n]);
rho=G*rho0;

V=null(G','r')';
d=size(V,1);

f=0;
while(f==0)
I=rand(n);

T=[I(:,1:n-d)'; V];
if(rank(T)==n)
    f=1;
end
end

J=T*rho*inv(T);

Jr=J(1:n-d,1:n-d);

flag=n-d-rank(Jr);

flag=1-flag;







