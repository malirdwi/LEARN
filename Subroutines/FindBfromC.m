function [B]=FindBfromC(varargin)
if nargin==2
G=varargin{1};
C=varargin{2};
elseif nargin==2
    A=varargin{1};
    B=varargin{2};
    G=B-A;
    C=varargin{3};
else
        error('Error: Wrong Number of Arguments');
        return;
end

B=[];
for k=1:size(C,1)
cvx_begin quiet
variable b(1,size(G,1))
b*G==C(k,:);
%b>=ones(1,size(G,1));
minimize(sum(abs(b)))
cvx_end

B=[B;b];
end
B=round(1e4*B)/1e4;
 