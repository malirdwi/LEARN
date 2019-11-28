function [C]=ConstructIterate(varargin)
if nargin==1
G=varargin{1};
auto=0;
elseif nargin==2
    A=varargin{1};
    B=varargin{2};
    G=B-A;
    auto=1;
else
        error('Error: Wrong Number of Arguments');
        return;
end

 
v=IsAS1(G);
if(min(v)>0)
    AS1=1;
else 
    AS1=0;
end

if AS1==0 && auto==0
    C=NoPositiveSS(G);
elseif AS1==0 && auto==1
    C=[]; %!!
elseif auto==1
    [C]=ConstructIterate(A,B);
else
    
C=[G;-G];

[n r]=size(G);
m=1;
kk=1;
G1=G.*((sign(G+1e-9)-1)/2);
G2=G.*(sign(G-1e-9)+1)/2;
 
  Gs=G1*G2';  
  llq=3;
while(m<500*n)
    if(length(C(:,1))<m)
      %  fprintf('SUCCESS at %u th iteration \n',m-1);
        C=RemoveRedundant(C);
        break;
    end
   Sin=find( abs(C(m,:)) > 0 ); 
 
   ww=C(m,:);
   for i=1:length(Sin)
       
       in=find(G(:,Sin(i))<0);
       for j=1:length(in)
           red=1;
                    sign(C(m,Sin(i)))*G(in(j),:);

         c1=C(m,:)+sign(C(m,Sin(i)))*G(in(j),:);
 
        for mm=1:length(C(:,1)) %check if redundant
          if(max(abs(c1-C(mm,:)))==0)
              red=0;
              break;
          end
        end
      %(max(abs(c1))>0)&&
        if((max(abs(c1))>0)&&(red>0))
         
           % c1
           MMM(llq)=m;
          llq=llq+1;
             C=[C;c1];  %% 
       
        end
       end
   end
     m=m+1;
 
  
end

if(length(C(:,1))>=m)
       C=[];
   %    fprintf('FAILURE WITH %u iterations \n',m-1);
       
end
if(rank(C)-rank(G)~=0)
C=[];
end

end
 
