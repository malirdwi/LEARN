function [C]=ConstructIterateA(A,B)
if nargin==1
error('Error: Wrong Number of Arguments');
end

 
G=B-A;  
C=G;
C=[C;-C];
 
 
[n r]=size(G);
m=1;
kk=1;
 
 
 
 
  llq=3;
while(m<500*n)
    if(length(C(:,1))<m)
        return;
    end
   Sin=find( abs(C(m,:)) > 0 ); 
   %disp('----')
   ww=C(m,:);
   for i=1:length(Sin)
       
       in=find(A(:,Sin(i))>0);
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
        if((max(abs(c1))>0)&&(red>0))
         
            
           MMM(llq)=m;
          llq=llq+1;
             C=[C;c1];  %% 
       
        end
       end
   end
     m=m+1;
   %kk=kk+1;
  
end

if(length(C(:,1))>=m)
       C=[];
   %    fprintf('FAILURE WITH %u iterations \n',m-1);
       
end
if(rank(C)-rank(G)~=0)
C=[];
end

end
 
