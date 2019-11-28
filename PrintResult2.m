function PrintResult2(C,G,conservative,siphons,file)
 if(length(C)==0)
         fprintf(file,'This method for constructing a PWL RLF has failed.\n')
    else
       fprintf(file,'Success!! A PWL RLF has been found.\n')
        
       fprintf(file,'The following is always a Lyapunov function for any monotone kinetics: V(x)=max_k c_k^T R, \n where C=[c_1^T,...,c_m^T]^T is given as follows:\n')
       
       PrintMatrix(C,file)          
       fprintf(file,'\n')
       
               fprintf(file,'The set of steady states is Lyapunov stable. Check the graphical LaSalle algorithm for global stability.\n')
        
 end
 