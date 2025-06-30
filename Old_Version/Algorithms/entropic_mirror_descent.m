
function f =  entropic_mirror_descent(J,s,T, alpha)
   f_1 = zeros(s,1)/s;
   f_t = repmat(f_1,1,T);
   f_t(:,1) = 1/s;
   eig_val_holder = zeros(1, T);
   for t = 1:T
       [V, d] = eigs(J(diag(f_t(:, t))), 1, 'largestreal');
       eig_val_holder(t) = d;
       u = V(1:length(V)/2);
       v = V(length(V)/2+1:end);
       theta = -alpha*(u.^2 + v.^2);
       beta = sqrt(2*log(s)/(T*norm(theta,'inf')^2));
       h = f_t(:,t).*exp(-beta.*theta);
       f_t(:,t+1) = h/sum(h);
   end
   
   [~, min_idx] = min(eig_val_holder);
   f = f_t(:,min_idx);
end