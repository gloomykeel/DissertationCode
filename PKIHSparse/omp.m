function xr = omp(n, Phi, y, s, sigma);

% initialization
res = y;
active_set = [];

i = 0;
while (i<s) && (norm(res) > sigma)

    % maximal correlation
    [val, idx] = max(abs(Phi'*res));

    % update the active set
    active_set = [idx active_set];

    % new residual
    res = y-Phi(:,active_set)*pinv(Phi(:,active_set))*y;
    i = i+1;
    
end

xr = zeros(n,1);
%
xr_active_set = pinv(Phi(:,active_set))*y;
xr(active_set) = xr_active_set;
