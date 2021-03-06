function flag=checkSiphonCondition(varargin)

if nargin==1
    auto=0;
G=varargin{1};
elseif nargin==2
    A=varargin{1};
    B=varargin{2};
    G=B-A;
     if(max(max(abs(A.*B)))>0)%non-autocatalytic
        auto=1;
     else
         auto=0;
     end
    
else
        error('Error: Wrong Number of Arguments');
        return;
end

file=1;
if auto==0
[siphons,deadlock]=checkSiphons(G);
else
[siphons,deadlock]=checkSiphons(A,B);
end

Mnetwork=checkMnetwork(G);
v=IsConservative(G);


if(min(v)>0)
    conservative=1; 

else
    conservative=0; 
end

if siphons==0
   fprintf(file,'The critical siphon necessary condition is satisfied.\n') ;
elseif deadlock==1
     fprintf(file,'There is a critical deadlock. A PWL RLF cannot exist.\n') ;
     elseif conservative==1 && Mnetwork==1
     fprintf(file,'This is a conservative M-network with a critical siphon. A PWL RLF cannot exist.\n') ;
     elseif conservative==1
     fprintf(file,'This is a conservative network with critical siphon. If there exists a positive non-degenerate steady state then a PWL RLF cannot exist.\n'); 
else
     fprintf(file,'There is a critical siphon. The necessary condition test is inconclusive but a PWL RLF is unlikely to exist.\n') ;
end
