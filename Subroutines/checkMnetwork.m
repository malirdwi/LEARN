function [flag,G2,I]=checkMnetwork(varargin)
if nargin==1
G=varargin{1};
elseif nargin==2
    A=varargin{1};
    B=varargin{2};
    G=B-A;
    if(max(max(abs(A.*B)))>0)%non-autocatalytic
        flag=0;
    end
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


if AS1==0
    flag=0;
    return;
end
[n,r]=size(G);

rev=[];
irrev=[];
tg=ones(1,r);
for j=1:r
    if(tg(j)==1)
    for k=j+1:r
        if( G(:,k)+G(:,j)==0)
            rev=[rev; j k];
            tg(j)=0;
            tg(k)=0;
            break;
        end
    end
    end
    
    if(tg(j)==1)
        irrev=[irrev;j];
    end
end
 
G2=[G(:,irrev)];

flag=0;
if(size(rev,1)==0)
    if(size(null(G2),2)~=1)
    flag=0;
    I=irrev;
    return;
end
    flag=1;
    I=irrev;
    return;
else
for j=0:2^size(rev,1)-1
    jb=1+dec2bin(j,size(rev,1))-48;
    I=[];
    for k=1:size(rev,1)
   
        G2=[G2 G(:,rev(k,jb(k)))];
        I=[I;rev(k,jb(k))];
    end
    
    v=IsAS1(G2);
    if min(v)>0
          flag=1;
        break;
    end
         
    G2=[G(:,irrev)];
    
end
end
if(flag==0)
    G2=[]; 
    return;
end
 

I=[irrev; I];

 
