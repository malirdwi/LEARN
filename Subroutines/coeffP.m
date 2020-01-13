function coeff = coeffP(vec,k)

 

[m, n] = size(vec);

 

if n == 1
   p = m;
   ff = 1;
else
   ff = 0;
end

if n == k
   coeff = vec(:).';
elseif n == k + 1
   tt = vec(:).';
   coeff   = tt(ones(n,1),:);
   coeff(1:n+1:n*n) = [];
   coeff = reshape(coeff,n,n-1);
elseif k == 1
   coeff = vec(:);
elseif   (k > 3 || n-k < 4) && n < 17
   rr = 2.^(n);
   nc = rr;

   for ct = 1:n
      ss = (0:1);
      nc = nc/2;
      nrp = rr./(2*nc);
      ss = ss(ones(1,nrp),:);
      ss = ss(:);
      ss = ss(:,ones(1,nc));
      x(:,n-ct+1) = ss(:);
   end

   index = x(sum(x,2) == k,:);
   nr = size(index,1);
   [rr,~] = find(index');
   coeff = reshape(vec(rr),k,nr).';
else 
   X = zeros(1,0); % want 1 row even for k=0
   if n == 1,
      vec = vec.';
   end
   if k < n && k > 1
      for index = 1:n-k+1
         Y = coeffP(vec(index+1:n),k-1);
         X = [X; [vec(ones(size(Y,1),1),index) Y]];
      end
   end
   coeff = X;
end
