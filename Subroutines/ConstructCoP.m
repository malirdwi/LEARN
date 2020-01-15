function [P,H2]=ConstructCoP(G,varargin)
checkToolbox
if nargin>1
H2=varargin{1};
else
    H2=[];
end
 warning off;

 
    


dR=(-G+abs(G))'/2 ; dR=sign(dR);

[jj,ii]=find(dR);
IJ=[jj,ii];
ml=size(IJ,1);
 E=zeros(size(dR,1),size(dR,2),ml);
for k=1:ml
  E(IJ(k,1),IJ(k,2),k)=1;
end




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
 


 % 
[Z,Z1,Zx,Z1v]=RegNbhd(S1,H);
 


 V=null(G,'r');


tol=0.5;
cvx_precision=[(2.22e-16)^0.5,(2.22e-16)^0.5,(2.22e-16)^0.375];
cvx_begin sdp quiet
cvx_quiet true
 variables Y(n,r,m/2) A(ne,ne,m*ml/2) B(ne,ne,m*ml/2) C(ne,ne,m/2) D(ne,ne,m/2)   LLQ(m,r,m)
 
 T=0;
 for j=1:m/2
     T=T+trace(Y(:,:,j)'*G);
 end
%  maximize(log_det(Y(:,:,1)'*G))
 minimize(T)
 
%[1:r]*(Y(:,:,1)'*G)*[1:r]'
%sum(sum(sum(D)))>tol;

for j=1:m/2
    
  
    
   (Y(:,:,j)'*G)== (Y(:,:,j)'*G)';
 
 
  
   C(:,:,j)== C(:,:,j)';
   D(:,:,j)== D(:,:,j)'; 
   
   
   C(:,:,j)>=0;
   

   
   
   (Y(:,:,j)'*G)*V==0;
   [1:r]*(Y(:,:,j)'*G)*[1:r]'>tol;
   
 
   
   
 diag(reshape(Y(:,:,j).*dR',[n r]))>=0;
   
   
   
      
   [(Y(:,:,j)'*G) ]>= [(diag(S1(j,:))*H)'*(C(:,:,j)+D(:,:,j))*(diag(S1(j,:))*H)  ];
   
   

for k=1:ml
      A(:,:,(j-1)*ml+k)== A(:,:,(j-1)*ml+k)';
      B(:,:,(j-1)*ml+k)== B(:,:,(j-1)*ml+k)';
      
      
      A(:,:,(j-1)*ml+k)>=0;
      diag(reshape(B(:,:,(j-1)*ml+k),[(ne)^2,1]))>=0;
   
   
   [(E(:,:,k)*G)'*(Y(:,:,j)'*G)+ (Y(:,:,j)'*G)*E(:,:,k)*G  ] +  [(diag(S1(j,:))*H)'*(A(:,:,(j-1)*ml+k)+B(:,:,(j-1)*ml+k))*(diag(S1(j,:))*H)  ] <=0 ;
  %  (diag(S1(1,:))*H)'*( (E(:,:,k)*G)'*P1+ P1*E(:,:,k)*G )*(diag(S1(1,:))*H) <=0 ;
end
   

end

 
P(:,:,1)=(Y(:,:,1)'*G);
 

for j=1:m
% diag(diag(LLQ(:,:,j)))>=0;

    for j2=j+1:m
        
      
     
    
        
         if(Zx(j,j2)>01)
             bb=find(abs(sign(S1(j,:)-S1(j2,:)))>0);
               if(j>m/2)
                     P1=(Y(:,:,m+1-j)'*G);
               else
                   P1=(Y(:,:,j)'*G);
               end
               
               if(j2>m/2)
                     P2=(Y(:,:,m+1-j2)'*G);
               else
                   P2=(Y(:,:,j2)'*G);
               end
       
            P1-P2 == H(bb,:)'*LLQ(j2,:,j) + LLQ(j2,:,j)'*H(bb,:); 
 
             end
        
            
    end






   
   
end
  


cvx_end
 

if(cvx_status(1)=='S' | cvx_status(3)=='a' | cvx_status(1)=='F');
    for j=1:m/2 
    P(:,:,j)=(Y(:,:,j)'*G);
   T=T+norm(P(:,:,j));
    end
     
   
%eig(P)
end
