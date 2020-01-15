function  T=Tlyap(j,m,r)
T1=-eye(m-1);
  
T=[ T1(1:j-1,:) ; ones(1,m-1); T1(j:m-1,:)];
 
