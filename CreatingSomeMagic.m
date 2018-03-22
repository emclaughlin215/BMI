N = size(trial,1);
K = size(trial,2);
magic = zeros(K,1);

for k=1:K
   for n=1:N
       magic(k,1) = magic(k,1)+norm(trial(n,k).handPos(1:2,end)-trial(n,k).handPos(1:2,1))/N;
   end
end