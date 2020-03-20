function [Q,cvx,H2,Xi]=ConstructLP(G,varargin)
checkToolbox
if nargin>1
H2=varargin{1};
else
    H2=[];
end
if nargin>2
w=varargin{2};
else
    w=1;
end

cvx=0;
if nargin>3
    cvx=varargin{3};
end
    
if w==1
    [Q,cvx,H2,Xi]=ConstructL1(G,H2);
    return;
end

  H=[G;H2  ];


     


V=null(G,'r');
 


[n r]=size(G);
[ne r]=size(H);
S=[0:2^ne-1]';
S= 2*(double(dec2bin(S))-48-0.5);
 options = optimoptions('linprog','Display','off');
 warning off;
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

 [Z,Z1,Zx,Z1v]=RegNbhd(S1,H);
 cvx_solver sedumi

  cvx_begin  quiet
variable Q(m/2,r)
variable L(m/2,ne)
variable L2(m,ne,m)

minimize(sum(sum(abs(Q))))
 
 
subject to
Q(:,r0).*C(1:m/2,r0)>=0;  %%% to have descreasing function
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
if cvx==1
    L2>=0
end
Lf=[L;L];
for j=1:m/2
       sum(L(j,:))>=2;
if w==1
    Q(j,:)==L(1,:)*diag(S1(j,:))*H; 
else
   Q(j,:)== L(j,:)*diag(S1(j,:))*H; 
end   
         Qf(m-j+1,:)=-Q(j,:);  
end
for j=1:m
    for j2=j+1:m
        bb=find(abs(sign(S1(j,:)-S1(j2,:)))>0);
          if(Z(j,j2))
         (Qf(j,:)-Qf(j2,:)) == L2(j2,bb,j)*diag(S1(j,bb))*H(bb,:);  
        end       
         
    end
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

Q=full([Q])/min(abs(Q(find(abs(Q)>1e-5))));
Q=Q(1:m/2,:);
 
if length(find(L2<-1e-8))>0
    cvx=0;
else
    cvx=1;
end
Xi=L;
