function [flag,deadlock]=checkSiphons(varargin)
if nargin==1
G=varargin{1};
oldparam = sympref('HeavisideAtOrigin',0);
 A=heaviside(-G);
B=heaviside(G);
elseif nargin==2
    A=varargin{1};
    B=varargin{2};
    G=B-A;
else
        error('Error: Wrong Number of Arguments');
        return;
end
    

 flag=0;
 deadlock=0;

 
[n,r]=size(G);
M=[sign(A) sign(B)  eye(n)];
for j=1:r
   y1=find(M(:,j+r)-M(:,j)>0);      
   y2=find(M(:,j)>0);              
   for c1=1:length(y1)     
      for c2=1:length(y2)  
                           
                                                     R = sign(M(y1(c1),:)+M(y2(c2),:));  

                          [mr,mc]=size(M);
                          F=0;
                          for i=1:mr,
                             if (R == M(i,:)),
                                F =1;            
                                break
                             end
                          end
                          if (F ==0)             
                             M = [M ; R];   
                          end
                       end
   end
   y1=sort(y1);
   for c=length(y1):-1:1
      if min(M(y1(c),1:r) - M(y1(c),r+1:2*r)) < 0,  
         M(y1(c),:) = [];
      end
   end
end
[mr,mc]=size(M);
 
   S=M(:,2*r+1:2*r+n);
   
   for i=1:size(S,1)
       
           s=find(S(i,:)>0);
           Gr=G(s,:);
           cvx_begin quiet
variable x(size(Gr,1))
Gr'*x==0
minimize(sum(x))
sum(x)>=1;
x>=0;
cvx_end

if(min(x)>=0)
else
    flag=1;
    if(deadlockcheck(S(i,:),A)==1)
        deadlock=1;
        return;
    end
        
        
 
end
 


   end
   
    

    function f=deadlockcheck(S,A)
       
        [n,r]=size(A);
        O=[1:r];
        I=find(S);
        for i=I
            O(find(A(i,:)))=0*O(find(A(i,:)));
        end
        if sum(O)==0
            f=1;
        else
            f=0;
        end
        
        
           
       
        









   

