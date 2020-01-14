function [x]=IsAS1(varargin)

if nargin==1
G=varargin{1};
elseif nargin==2
    A=varargin{1};
    B=varargin{2};
    G=B-A;
else
        error('Error: Wrong Number of Arguments');
        return;
    end
 

cvx_begin quiet
variable x(size(G,2))
G*x==0
minimize(sum(x))
x>=ones(size(G,2),1)
cvx_end

if(min(v)>0)
else
v=0;
end
