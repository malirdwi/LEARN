function C=ConstructGraphical(G)

[flag,G2,I]=checkMnetwork(G);
C=[];
if flag==0
    return;
end


[n,r]=size(G);
[n,r2]=size(G2);


for i=1:n
    if(length(find(G2(i,:)<0))>1)
        return;
    end
end

for j=r2-(r-r2)+1:r2
    In=find(G2(:,j)>0);
    for i=In
        
    if(length(find(G2(i,:)>0))>1)
        return;
    end
    
    end
end

for i=1:size(G2,2)
    for j=i+1:size(G2,2)
        
        c=zeros(1,r);
        
        c(I(i))=1;
        c(I(j))=-1;
        
        b1=FindRev(G,I(i));
        b2=FindRev(G,I(j));
        if(b1~=0)
            c(b1)=-1;
        end
        if(b2~=0)
            c(b2)=1;
        end
        
        C=[C;c];
    end
end






        




        