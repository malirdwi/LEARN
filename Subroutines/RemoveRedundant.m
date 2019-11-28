function C1=RemoveRedundant(C)

[m,r]=size(C);

I=[1:m];

i=1;

for i=1:m
 if(C(i,1)<Inf)
    for j=i+1:m
        if( C(i,:)+C(j,:)==0)
            C(j,1)=Inf;
        end
    end
 end
end

C1=[];
for i=1:m
 if(C(i,1)<Inf)
    C1=[C1;C(i,:)];
    end
 end
end
 
