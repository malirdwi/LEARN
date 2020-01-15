function [Q,cvx,H2,Xi]=ConstructL1(G,H2)

    checkToolbox

  cvx=1;

  H=[G;H2  ];


     


V=null(G,'r');
 


[n r]=size(G);
[ne r]=size(H);
S=[0:2^ne-1]';
S= 2*(double(dec2bin(S))-48-0.5);
 options = optimoptions('linprog','Display','off');

S1=[];
j=1;jj=1;
for i=1:2^ne   %%% REMOVING INCONSISTENT SYSTEM OF INEQUALITES
   
 [x,y,flag]=linprog(zeros(1,r),-diag(S(i,:))*H,-1e-5*ones(ne,1),[],[],zeros(r,1),[],[],options);
 if(flag>0)
     S1(j,:)=S(i,:);
    j=j+1;
 else
     S2(jj,:)=S(i,:);
     jj=jj+1;
 end
end
%S1=S; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"!!!!!!!!!!!!!!!!%%%%%%%%%%
[m ne]=size(S1);
 


C=[];
r0=[];
for j=1:m  %%% EXCLUDING RATES WITH MULTIPLE CONFLICTING INPUTS
   
    for l=1:r
        In=abs(min([zeros(n,1) G(:,l)]')).*S1(j,1:n);
        if(max(In)*min(In)<0)  %be careful from zero!
            C(j,l)=0;
            if(length(find(r0==l))==0)
                r0=[r0 l];
            end
        elseif(max(In)>0)
            C(j,l)=-1;
            if(length(find(r0==l))==0)
                r0=[r0 l];
            end
        elseif(min(In)<0)
            C(j,l)=1;
            if(length(find(r0==l))==0)
                r0=[r0 l];
            end
        else
            C(j,l)=NaN; %%reaction with no inputs
            
            
        end
    end
    
        
end
 
 
s=length(find(C>0));
q=sym('q',[s 1]);
cvx_solver sedumi

  cvx_begin quiet
variable Q(m/2,r)
variable L(1,ne)

minimize(sum(sum(abs(Q))))
 
 
subject to
 
 
Q(:,r0).*C(1:m/2,r0)>=0;  %%% to have descreasing function
% 
for j=1:m/2
    for l=1:r
        if(C(j,l)==0)
            Q(j,l)==0;   %% zeros what have to be zeroed
        end
    end
end 

 Q*V==0;  
 L>=0;
 
Qf=Q;


Lf=[L;L];
for j=1:m/2
    

       sum(L(1,:))>=2;

    Q(j,:)==L(1,:)*diag(S1(j,:))*H; 
  
          
end




cvx_end 


if sum(cvx_status)==sum('Failed')
    Q=[];
    Xi=[];
    return;
end
if sum(cvx_status)==sum('Infeasible')
    Q=[];
    Xi=[];
    return;
end
Xi=L;
Q=full([Q])/min(abs(Q(find(abs(Q)>1e-5))));
Q=Q(1:m/2,:);
 



