function H2=AddRowtoH(G)

H=G;
H2=[];

[n,r]=size(G);


for i=1:n
    
    for j=i+1:n
          h=G(i,:)-G(j,:);
        
        for k=1:size(H,1)
            
          h2=H(k,:);
          a1=find(h);
          a2=find(h2);
          flag=1;
          if(length(a1)==length(a2))
             if(norm(a1-a2)==0)
                 zz=h(a1)./h2(a2);
                 zz=zz/zz(1);
                 if norm(zz-ones(1,length(zz)))==0
                     flag=0;
                     break;
                 end
             end
          end
        end
          if flag==1
              H=[H;h];
              H2=[H2;h];
          end
        
    end
end
          
