function PrintMatrix(C,file)

       for ii = 1:size(C,1)
           for jj=1:size(C,2)
               if(C(ii,jj)>=0)
     fprintf(file,' %g\t',C(ii,jj));
               else
                fprintf(file,'%g\t',C(ii,jj));   
               end
           end
    fprintf(file,'\n');
end
    
