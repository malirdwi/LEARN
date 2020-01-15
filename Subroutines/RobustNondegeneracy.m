function flag=RobustNondegeneracy(varargin)
checkToolbox

if nargin==1
G=varargin{1};
oldparam = sympref('HeavisideAtOrigin',0);
 A=heaviside(-G);
B=heaviside(G);
elseif nargin==2
    A=varargin{1};
    B=varargin{2};
    G=B-A;
else
        error('Error: Wrong Number of Arguments');
        return;
end

[n,r]=size(G);
oldparam = sympref('HeavisideAtOrigin',0);

rho0=A'.*ones(r,n);%sym('r',[ r n]);
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







