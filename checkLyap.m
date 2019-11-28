function status=checkRLF(G,C)% 




[n r]=size(G);
 
 
 m=length(C);
 
status=1;
 
 for j=1:m
 C1=[];
 for i=1:m
     if(i~=j)
     C1=[C1; C(j,:)-C(i,:)];
     end
 end
  C1=C1';
  
 Sp=[]; Sm=[];
 
 for i=1:n
    
     [x,fval,flag]=linprog(zeros(m-1,1),[],[],C1,G(i,:)',zeros(m-1,1));

if(flag==1)
        Sp=[Sp i];
end

[x,fval,flag]=linprog(zeros(m-1,1),[],[],-C1,G(i,:)',zeros(m-1,1));

if(flag==1)
      Sm=[Sm i];
end
    

 end
 

for l=1:r
    W=find(G(:,l)<0);
    
    if(C(j,l)>1e-8)
       
        for w=1:length(W)
            t3=find(Sm==W(w));
            if(length(t3)==0)
                status=0;
                break;
            end
        end
        
    else if(C(j,l)<-1e-8)
          for w=1:length(W)
            t3=find(Sp==W(w));
            if(length(t3)==0)
                status=0;
                break;
            end
        end
    end
 end
 if(status==0)
     break;
 end
end
 
 if(status==0)
     break;
 end
 
 end
 
  
 %clc
  disp('==========================================================');
if(status==0);
    disp('NOT VALID LYAPUNOV FUNCTION');
    fprintf(1, ' Check with c%u, R%u, with input x%u \n',j,l,W(w));
else
    disp('VALID LYAPUNOV FUNCTION');
end

 
% clear x
% cvx_begin
% variable x(m,1)
% minimize(0);
% subject to
% C1*x<=g;
% x>=0;
% 
% cvx_end 
% x