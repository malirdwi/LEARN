 function status=checkRLFV2(G,C)% 
checkToolbox
AS1=IsAS1(G);
if(AS1==0)
    fprintf('This subroutine checks RLFs for networks that have a positive vector in the kernel of the stoichiometri matrix. The given network does not satsify that.');
    return;
end

C=RemoveRedundant(C);
C=[C;-C];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [n r]=size(G);
  
 
C=C';
 
[n r]=size(G);
G1=(sign(G+1e-9)-1)/2;

 v=null(G); v=v/max(v);
 
 Ct=C(1:r-1,:);
 m=size(C);
 m=m(2);
status=1;
 
 for j=1:m/2
     
 
  C1=  C * Tlyap(j,m,r) ;
  

%       cvx_begin  %%check nonnegativity
% variable V(m-1,1)
% minimize(0);
% subject to
% C1*V<=C(:,j);
% V>=1e-7*ones(m-1,1);
% 
% cvx_end 
% 
% if(1-strncmpi(cvx_status,'Solved',6))
%     break;
% end
 

for l=1:r
 
    cvx_begin quiet
variable L(m-1,n)
minimize(0);
subject to
C1*L<=G'*diag(C(l,j)*G1(:,l));
L>=0;

cvx_end
if(1-strncmpi(cvx_status,'Solved',6))
    if(1-strncmpi(cvx_status,'Inaccurate/Solved',17))
    break;
end
end
 
if(1-strncmpi(cvx_status,'Solved',6))
    if(1-strncmpi(cvx_status,'Inaccurate/Solved',17))
    break;
end
end
 
 end
 end
  
% clc
  disp('==========================================================');
if(1-strncmpi(cvx_status,'Solved',6));
    status=0;
    disp('NOT VALID LYAPUNOV FUNCTION');
  %  fprintf(1, ' Check with c%u, R%u, with input x%u \n',j,l,W(w));
else
    disp('VALID LYAPUNOV FUNCTION');

end

 
