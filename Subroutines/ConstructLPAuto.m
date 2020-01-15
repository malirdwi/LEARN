function [Q,cvx,H2,Xi]=ConstructLPAuto(A,B,varargin)
 if nargin>2
H2=varargin{2};
else
    H2=[];
end
if nargin>3
w=varargin{3};
else
    w=0;
end

cvx=0;
if nargin>4
    cvx=varargin{4};
end


  
    
 
    G=B-A;
    
    v=IsAS1(G);
if(min(v)>0)
    AS1=1;
else 
    AS1=0;
    Q=[];
    return;
end

 
  H=[G;H2  ];


     

 V=null(G,'r');
 


[n r]=size(G);
[ne r]=size(H);
S=[0:2^ne-1]';
S= 2*(double(dec2bin(S))-48-0.5);

S1=[];
j=1;jj=1;
for i=1:2^ne   %%% REMOVING INCONSISTENT SYSTEM OF INEQUALITES
    options = optimoptions('linprog','Display','off');
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
 
 

%%%%%%%%% REFINING THE PARTITIONS



C=[];
r0=[];
for j=1:m  %%% EXCLUDING RATES WITH MULTIPLE CONFLICTING INPUTS
   
    for l=1:r
        In=max([zeros(n,1) A(:,l)]').*S1(j,1:n);
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


 

Z2=Z1;
  cvx_begin quiet
variable Q(m,r)
variable L(m,ne)
variable L0(m,r)
variable L2(m,ne,m)

    minimize(sum(sum(Q(1,1))))
 
 
subject to
 
uy=[ 1:m/2     ];% 4 5 6  8
uy=[uy m+1-flip(uy)];
Q(uy,r0).*C(uy,r0)>=0;  %%% to have descreasing function
% 
for j=uy
    for l=1:r
        if(C(j,l)==0)
            Q(j,l)==0;   %% zeros what have to be zeroed
        end
    end
end 
 

 Q*V==0;   % Kernel of the function must be subset of the ketnel of G
 L>=0;
 
L0==0; %%%%%%%%%%%!!!!!!!!!!!!!!!!
 
for j=1:m
    
 
       sum(L(j,:))+sum(L0(j,:))>=1; 
 
if w==1
    Q(j,:)==L(1,:)*diag(S1(j,:))*H; 
else
   Q(j,:)== L(j,:)*diag(S1(j,:))*H; 
end
 
     if(j<=m/2)
         Q(j,:)==-Q(m-j+1,:);  
     end
     %%%%%% continu9ety  covexity
     ex=[       ];
   ex=[ex wrev(m+1-ex)];
     Z2(:,ex)=0*Z2(:,ex);
     Z2(ex,:)=0*Z2(ex,:);
     
 
   
   
   
    for j2=j+1:m
        
     
    
        bb=find(abs(sign(S1(j,:)-S1(j2,:)))>0);
 
            P=null(H(bb,:));
             if(length(P)>0)
             (Q(j,:)-Q(j2,:))*P==0;
             end
 
            
    end
end
    
if cvx==1
    L2>=0
end


cvx_end 
 


 
if sum(cvx_status)==sum('Failed')
    Q=[];
    return;
end
if sum(cvx_status)==sum('Infeasible')
    Q=[];
    return;
end

Q=round(Q,5);
Q=full([Q])/min(abs(Q(find(Q))));
null(Q,'r');

if length(find(L2<-1e-8))>0
    cvx=0;
else
    cvx=1;
end
Xi=L;
