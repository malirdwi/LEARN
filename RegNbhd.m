function[Z,Z1,Zx,Z1v]=RegionNeighbors(S1,G);

[m r]=size(S1);
[n r]=size(G);
 
for j=1:m
    for j2=1:m
        Z(j,j2)=length(find(abs(sign(S1(j,:)-S1(j2,:)))>0));
        
    end
    
end

Z1=sparse(0.5*((-Z+2)+abs(-Z+2)-4*eye(m)));
Zx=sparse(0*Z1);
[x y]=find(Z1);
for j=1:length(find(Z1)')
          Zx(x(j),y(j))=find(abs(sign(S1(x(j),:)-S1(y(j),:)))>0) ;
    
end
     




Z1v=sparse(Z1);
dR=sign(0.5*(-G+abs(-G)))';

for j=1:m
    for q=find(Z1(j,:)) 
        
        ss=[];
        for l=1:r
            for i=1:n
                if(i~=Zx(j,q))
                     ss=[ss G(Zx(j,q),l)*dR(l,i)*S1(j,i)];
                end
               
            end
        end
        
        if(max(ss)*min(ss)>=0)
            s1=max(ss(find(ss)));
               s2=sign(S1(j,Zx(j,q)));
               if(s1*s2>0)
                Z1v(j,q)=0;
               end
        end
        
       
    end
    
end