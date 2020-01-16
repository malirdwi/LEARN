function LEARN(varargin)

if ~exist('ConstructGraphical','file')
 if ~exist('coeffP','file')
        addpath('./Subroutines');
    end
 end
    
checkToolbox
file=1;
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
    
fprintf(file,'--------------------------------\n');
fprintf(file,'Welcome to LEARN v1.0, Jan 2020\n');
fprintf(file,'Developed by M. Ali Al-Radhawi malirdwi@{northeastern.edu,mit.edu,gmail.com}\n\n');
fprintf(file,'LEARN tries to construct a Robust Lyapunov Function for a given reaction network.\n');
fprintf(file,'--------------------------------\n');
[n,r]=size(G);


fprintf(file, 'The network has %d species and %d reactions.\n',n,r);

fprintf(file, 'The stoichiometric space is %d-dimensional.\n',rank(G));

v=IsAS1(G);

if(min(v)>0)
    AS1=1;
    fprintf(file, 'The network has a positive vector in the kernel of the stoichiometry matrix, i.e. it has the potential for positive steady states.\n');

else
    AS1=0;
    fprintf(file, 'The network does not have the potential for positive steady states, i.e. there is no positive vector in the kernel of the stoichiometry matrix\n');
end


v=IsConservative(G);


if(min(v)>0)
    conservative=1;
    fprintf(file, 'The network is conservative. \n');

else
    conservative=0;
    fprintf(file, 'The network is not conservative.\n');
end


if auto==0
[siphons,deadlock]=checkSiphons(G);
elseif auto==1
    [siphons,deadlock]=checkSiphons(A,B);
end

if(siphons==1)
    fprintf(file, 'The network has critical siphons. It is NOT structurally persistent. \n');
else
    fprintf(file, 'The network has no critical siphons. It is structurally persistent. \n');
end





fprintf(file,'--------------------------------\n');

fprintf(file,'LEARN will check some necessary conditions\n');
Mnetwork=checkMnetwork(G);

if auto==0
Mnetwork=checkMnetwork(G);
elseif auto==1
    Mnetwork=0;
end
siph1=1;
fprintf(file,'Necessary Condition # 1 ....\n');
if siphons==0
   fprintf(file,'The critical siphon necessary condition is satisfied.\n') ;
elseif deadlock==1
     fprintf(file,'There is a critical deadlock. A PWL RLF cannot exist.\n') ;
     siph1=0;
     elseif conservative==1 && Mnetwork==1
     fprintf(file,'This is a conservative M-network with a critical siphon. A PWL RLF cannot exist.\n') ;
     siph1=0;
     elseif conservative==1
     fprintf(file,'This is a conservative network with critical siphon. If there exists a positive non-degenerate steady state then a PWL RLF cannot exist.\n');
else
     fprintf(file,'There is a critical siphon. The necessary condition test is inconclusive but a PWL RLF is unlikely to exist.\n') ;
end

fprintf(file,'Necessary Condition # 2 ....\n');
if auto==0
[n1,Q]=SignPatternCheck(G);
elseif auto==1
   [n1,Q]=SignPatternCheck(A,B);
end 

if(n1==1)
    fprintf(file,'The sign pattern necessary condition is satisfied.\n');
else
        fprintf(file,'The sign pattern necessary condition is violated. A PWL RLF does not exist\n');
   %     return;
end

fprintf(file,'Necessary Condition # 3 ....\n');

if auto==0
Pmatrix=checkPmatrix(G);
elseif auto==1
   Pmatrix=checkPmatrix(A,B);
end 
 

if(Pmatrix==1)
    fprintf(file,'The P matrix necessary condition is satisfied. \n');
else
        fprintf(file,'The P matrix necessary condition is violated. A PWL RLF does not exist\n');
     %   return;
end

necessary=siph1*Pmatrix*n1;
fprintf(file,'--------------------------------\n');
fprintf(file,'LEARN will search for a PWL RLF\n');

fprintf(file,'Method # 1: Graphical Method .. \n') ;

if(Mnetwork==1)
   fprintf(file,'This is an M-network. The graphical criteria will be checked \n\n');
   C=ConstructGraphical(G);
   PrintResult(C,G,conservative,siphons,file);
else
    fprintf(file,'This is not an M-network. Method # 1 is not applicable. \n\n');
end
fprintf(file,'--------------------------------\n');
fprintf(file,'Method # 2: Iterative Method .. \n');


if auto==0
C=ConstructIterate(G);
elseif auto==1
   C=ConstructIterate(A,B);
end 
    
    if(AS1==1)
    PrintResult(C,G,conservative,siphons,file);
    else
           PrintResult2(C,G,conservative,siphons,file);
    end
    
    
  

fprintf(file,'--------------------------------\n');
fprintf(file,'Method # 3: Linear Programming Method .. \n');
if(AS1==1)
fprintf(file,'The partition matrix H is set to the default choice H=the stoichiometry matrix .. \n') ;   
    if auto==0
   [C,cvx,H2,Xi]=ConstructLP(G);
    elseif auto==1
[C,cvx,H2]=ConstructLPAuto(A,B,[],1,0);
    end
      if cvx==0
       cvx=checkRLF_quiet(G,round(C,6));
   end
   if cvx==1 || length(C)==0
       PrintResult(round(C,5),G,conservative,siphons,file);
       if(length(C)>0)
       fprintf(file,'Please note that this function is a Sum-of-Currents RLF which can alternatively be written \n') ;   
       fprintf(file,'  as V(x)= sum_i xi_i |dot x_i|, where xi=[xi_1 .... xi_n]=   \n') ;
       PrintMatrix(Xi/min(abs(Xi(find(abs(Xi)>1e-5)))),file);
       end
   else
    fprintf(file,'SUCCESS!! A continuous PWL RLF has been found .. \n') ;   
    fprintf(file,'The set of steady states is Lyapunov stable .. \n')  ;
        fprintf(file,'The coefficient matrix is given as follows .. \n') ;  
        PrintMatrix(round(C,5),file);
        if(length(H2)>0)
         fprintf(file,'The parititon matrix is given as follows .. \n') ;
         PrintMatrix(H2,file);
        end
   end

else
    fprintf(file,'This method for constructing a PWL RLF has failed.\n')
end





       


fprintf(file,'THE END.\n')
 
