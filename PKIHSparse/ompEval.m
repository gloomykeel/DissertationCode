function [ sucess ] = ompEval( N,k,r,d,q,type,Psi,trial_num,o )
%WHT1D Summary of this function goes here
%   Detailed explanation goes here
sucess = 0;
for i=1:trial_num
    % create a sparse signal in Psi domain
    alp = [randn(k,1); zeros(N-k,1)];
    p = randperm(N);
    alp = alp(p);
    x = Psi*alp;
    
    p = ceil(N/(q*d));
    pvec = randperm(N);
    switch(type)
        case 1
            qvec = randperm(N);
            H = kron(eye(p*q),hadamard(d)/sqrt(d)); % base matrix
            A = H(pvec(1:r),qvec); % srmg
        case 2
            H = kron(eye(p*q),hadamard(d)/sqrt(d)); % base matrix
            A = H(pvec(1:r),1:N)*sparse(1:N,1:N,2*round(rand(N,1))-1); % srml full
        case 3
            kih = kron(eye(p),kron(hadamard(d)/sqrt(d),sparse(eye(q))));
            rv = reshape(1:(ceil(N/o)*o),o,[])';
            A = kih(pvec(1:r),rv(rv<=N)); % pkih
        case 4
            A = orth(randn(N));
            A = A(pvec(1:r),:);
    end
    
    % observation
    y = A*x;
    
    % reconstruction using OMP alg.
    xr = omp(N, A*Psi, y, k, 0);
    
    % if snr >50, it is considered as perfect reconstruction
    if Qsnr(alp,xr)>50
        sucess = sucess+1;
    end
end

end

