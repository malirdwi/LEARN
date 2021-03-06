function [x]=IsConservative(varargin)
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
variable x(size(G,1))
G'*x==0
minimize(sum(x))
x>=ones(size(G,1),1)
cvx_end

if(min(x)>0)
else
x=[];
end
