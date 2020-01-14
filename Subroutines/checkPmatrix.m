function flag=checkPmatrix(varargin)
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


flag=1;
[n,r]=size(G); 

rho=A'.*sym('r',[ r n]); 


J=-G*rho;

 
 
n=size(J,1);
for j=1:n
    
    Comb=coeffP(1:n,j);
    m=size(Comb,1);
    
    for k=1:m
        WAQ=expand(det(J(Comb(k,:),Comb(k,:))));
        if(WAQ~=sym(0))
           % coeffs(WAQ)
              if(min(double(coeffs(WAQ)))<0)
                flag=0;
                return;
            end
        end
    end
    
end
