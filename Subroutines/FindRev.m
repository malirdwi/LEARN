function b=FindRev(G,a)

[n,r]=size(G);

for i=1:r
    if(G(:,a)+G(:,i)==0)
        b=i;
        return;
    end
end

b=0;
