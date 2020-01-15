function status=checkRLF(G,C)% 
checkToolbox
file=1;

fprintf(file,'--------------------------------\n');
fprintf(file,'Welcome to LEARN v1.0, Jan 2020\n');
fprintf(file,'Developed by M. Ali Al-Radhawi malirdwi@{northeastern.edu,mit.edu,gmail.com}\n\n');
fprintf(file,'This subroutine tries to verify if a given function is a Robust Lyapunov Function for a given reaction network.\n');
fprintf(file,'Rarely, this subroutine faces numerical problems due to the LP solver linprog. \n');
fprintf(file,'If you witness an unexpected output, please use the alternative subroutine checkRLF2 which uses cvx.\n');
fprintf(file,'--------------------------------\n');


AS1=IsAS1(G);
if(AS1==0)
    fprintf('This subroutine checks RLFs for networks that have a positive vector in the kernel of the stoichiometri matrix. The given network does not satsify that.');
    return;
end

C=RemoveRedundant(C);
C=[C;-C];

[n r]=size(G);
 
 
 m=length(C);
 
status=1;
  options = optimoptions('linprog','Display','off','Algorithm','interior-point-legacy');
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
    
     [x,fval,flag]=linprog(zeros(m-1,1),[],[],C1,G(i,:)',zeros(m-1,1),[],options);

if(flag==1)
        Sp=[Sp i];
end

[x,fval,flag]=linprog(zeros(m-1,1),[],[],-C1,G(i,:)',zeros(m-1,1),[],options);

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
    fprintf(1,'NOT VALID LYAPUNOV FUNCTION');
    fprintf(1, ' Check with c%u, R%u, with input x%u \n',j,l,W(w));
else
    fprintf(1,'SUCCESS!! VALID RLF');
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
