function status=checkRLF(G,C)% 
checkToolbox
file=1;

 status=0;


AS1=IsAS1(G);
if(AS1==0)
   % fprintf('This subroutine checks RLFs for networks that have a positive vector in the kernel of the stoichiometri matrix. The given network does not satsify that.');
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
 
 
 
