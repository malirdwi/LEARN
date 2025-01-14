function L2=FindLambda(varargin)
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


  [n r]=size(G);
 oldparam = sympref('HeavisideAtOrigin',0);
 A=heaviside(-G);
B=heaviside(G);


dR=A';
 
 [jj ii]=find(dR);


 [m r ]=size(C);
 s=length(ii);

Q=zeros(r,r,s);
Jq=zeros(n,n,s);
  L2=zeros(m,m,s);

 for l=1:s

     I=eye(r);
        Q(:,:,l)=I(:,jj(l))*G(ii(l),:);
             I2=eye(n);
        Jq(:,:,l)=G(:,jj(l))*I2(ii(l),:);


 TT=C*I(:,jj(l))*G(ii(l),:);
 if(norm(TT)~=0)
     cvx_begin quiet
            variables L11(m,m)  %symmetric


dL11=(diag(L11));
oL11=L11-diag(dL11);


                 sum(abs(oL11'))<=-(dL11') %linf
              minimize(sum((dL11')+sum(abs(oL11')))) %linfobjective


    TT== (L11)*C;


 cvx_end
 else

 end
   L2(:,:,l)=L11;
  if(sum(cvx_status)~=sum('Solved')&&sum(cvx_status)~=sum('Inaccurate/Solved'))
  error('The provided function is not a GLF for the provided network. double check')
  break;
   L2(:,:,l)=NaN*L11;
 end
 end

r=sym('rho',[s,1]);
assume(r>=0);;
 J=zeros(m,m);


  L2=round(1e4*L2)/1e4;
  fprintf(1,'The Jacobian G*dR/dx satisfies  C*G*dR/dx= \sum_{ell=1}^s dR_{j_l}/dx_{i_l} Lambda_l x')
for j=1:s
fprintf(1,'The matrix Lambda_%d multiplying the partial derivative dR%d/dx%d is \n',j,jj(j),ii(j))
L2(:,:,j)
J=J+r(j)*(L2(:,:,j));
end


  fprintf(1,'The weighted infinity norm of is upper bounded by the max of the following vector')

vpa(sum( (diag(diag(J))+abs(J-diag(diag(J))))')',3)

 x=input('enter any number of to save the tensor containing the Lambda matrices');
