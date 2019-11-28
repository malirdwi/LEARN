function PrintResult2(C,G,conservative,siphons,file)
 if(length(C)==0)
         fprintf(file,'This method for constructing a PWL RLF has failed.\n')
    else
       fprintf(file,'Success!! A PWL RLF has been found.\n')
        
       fprintf(file,'The following is always a Lyapunov function for any monotone kinetics: V(x)=|| C*R(x) ||_infty, \n where C is given as follows:\n')
       
       PrintMatrix(C,file)          
       fprintf(file,'\n')
       
       f3=RobustNondegeneracy(G);
       if(f3==0)
            
           if conservative==1 && siphons==0
               fprintf(file,'The robust non-degeneracy test is passed. \n Since the network is conservative and with no critical sihpons then the following holds: \n There exists a unique postive globally asymptotically stable steady state  in each stoichiometric class.\n')
           else
              fprintf(file,'The robust non-degeneracy test is passed. \nAny positive steady state is a unique globally asymptotically stable relative to its stoichiometric class.\n')
           end 
           
           
       else
           fprintf(file,'The set of steady states is Lyapunov stable.\n')
       end
 end
 
