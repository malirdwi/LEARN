function coeff = coeffP(vec,k)

 

p = length(vec);

 
 
if p == k
   coeff = vec(:).';
elseif p == k + 1
   tt = vec(:).';
   coeff   = tt(ones(p,1),:);
   coeff(1:p+1:p*p) = [];
   coeff = reshape(coeff,p,p-1);
elseif k == 1
   coeff = vec(:);
elseif   (k > 3 || p-k < 4) && p < 17
   rr = 2.^(p);
   nc = rr;

   for ct = 1:p
      ss = (0:1);
      nc = nc/2;
      nrp = rr./(2*nc);
      ss = ss(ones(1,nrp),:);
      ss = ss(:);
      ss = ss(:,ones(1,nc));
      x(:,p-ct+1) = ss(:);
   end

   index = x(sum(x,2) == k,:);
   nr = size(index,1);
   [rr,~] = find(index');
   coeff = reshape(vec(rr),k,nr).';
else 
   X = zeros(1,0); % want 1 row even for k=0
 
      vec = vec.';
  
   if k < p && k > 1
      for index = 1:p-k+1
         Y = coeff3(vec(index+1:p),k-1);
         X = [X; [vec(ones(size(Y,1),1),index) Y]];
      end
   end
   coeff = X;
end
