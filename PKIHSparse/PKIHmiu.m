for j=[128 256 512 1024]
N=j;
d=4;
q=5;
o=3;
psi = dctmtx(N);
p = ceil(N/(q*d));
mius = zeros(100,4);
for i=1:100
    pvec = randperm(N);
    pvec=pvec(1:N);
    qvec = randperm(N);
    H = kron(eye(p*q),hadamard(d)/sqrt(d)); % base matrix
    A = H(pvec,qvec); % srmg
    b=A*psi;
    mius(i,1) = max(max(b));
    
    H = kron(eye(p*q),hadamard(d)/sqrt(d)); % base matrix
    A = H(pvec,1:N)*sparse(1:N,1:N,2*round(rand(N,1))-1); % srml full
    b=A*psi;
    mius(i,2) = max(max(b));
    
    kih = kron(eye(p),kron(hadamard(d)/sqrt(d),sparse(eye(q))));
    rv = reshape(1:(ceil(N/o)*o),o,[])';
    A = kih(pvec,rv(rv<=N)); % pkih
    b=A*psi;
    mius(i,3) = max(max(b));
    
    A = orth(randn(N));
    b=A*psi;
    mius(i,4) = max(max(b));
end
mm = mean(mius);
hold on
plot(j,mm(1),'o')
plot(j,mm(2),'s')
plot(j,mm(3),'*')
plot(j,mm(4),'^')
end