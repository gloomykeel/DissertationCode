function [ Hand,HandT ] = getMtxHandle( mtxName,N,K,blkSize,m,n )
%GETMTXHANDLE 根据不同的输入参数，返回不同的测量矩阵句柄
%   此处显示详细说明

mtxName = upper(mtxName);
switch(mtxName)
    case 'SRMG'
        % partial block Hadamard with random permutation;
        mtxB = hadamard(blkSize)/sqrt(blkSize);
        p = randperm(N);
        q = randperm(N);
        mtx = kron(eye(ceil(N./blkSize)),mtxB);
        mtx = mtx(q(1:K),p);
    case 'SRML'
        % partial Hadamard with random sign reversal;
        mtx = hadamard(blkSize)/sqrt(blkSize);
        p = 2*round(rand(N,1))-1;
        mtx = mtx*diag(p);
        q = randperm(N);
        mtx = mtx(q(1:K),:);
    case 'KRON'
        I = eye(ceil(N./blkSize));
        H = hadamard(blkSize);
        p = randperm(ceil(N./blkSize));
        Km = round(sqrt(K*ratio));
        Kn = round(K/Km);
        I = I(p(1:Km),:);
        mtx = kron(I,zeros(Kn,blkSize));
        for i=1:Km
            tmp = randperm(blkSize);
            tmp = [tmp,tmp];
            num = find(tmp==1);
            tmp = tmp(num:(num+Kn-1));
            tmpH = H(tmp,:);
            tMtx = kron(I(i,:),tmpH);
            mtx(((i-1)*Kn+1):i*Kn,:) = tMtx;
        end
    case 'TOEP'
        tM = 40;
        tN = 64;
        blkM = K/tM;
        blkN = N/tN;
        T = cell(blkM,blkN);
        for i=1:(blkM*blkN)
            tVec = zeros(tN,1);
            tp = randperm((tN-tM+1));
            %     tp = randperm(tN);
            tVec(tp(1)) = (rand<0.5)*2-1;
            tmp = zeros(tM,1); tmp(1) = tVec(1);
            T{i} = circshift(toeplitz(tmp,tVec),ceil(rand*tN),2);
            %     T{i} = toeplitz(tmp,tVec);
        end
        T = reshape(T,blkM,[]);
        mtx = cell2mat(T);
    case 'PKUH'
        Km = round(sqrt(K*m/n));
        Kn = round(sqrt(K*n/m));
        mtxB = hadamard(blkSize)/sqrt(blkSize);
        mtx1 = kron(eye(ceil(m./blkSize)),mtxB);
        mtx2 = eye(n);
        p = randperm(m);
        q = randperm(m);
        tmp = [q,q];
        num = find(tmp==1);
        tmp = tmp(num:(num+Km-1));
        mtx1 = mtx1(tmp,p);
        p = randperm(n);
        q = randperm(n);
        mtx2 = mtx2(q(1:Kn),p);
        mtx = kron(mtx2,mtx1);
    case 'KGAU'
        Km = round(sqrt(K*m/n));
        Kn = round(sqrt(K*n/m));
        mtx1 = orth(randn(m));
        mtx2 = orth(randn(n));
        mtx1 = mtx1(1:Km,:);
        mtx2 = mtx2(1:Kn,:);
        mtx = kron(mtx1,mtx2);
    case 'KSRMG'
        mtxB = hadamard(blkSize)/sqrt(blkSize);
        Km = round(sqrt(K*m/n));
        Kn = round(sqrt(K*n/m));
        mtx1 = kron(eye(ceil(m./blkSize)),mtxB);
        mtx2 = kron(eye(ceil(n./blkSize)),mtxB);
        p = randperm(m);
        q = randperm(m);
        tmp = [q,q];
        num = find(tmp==1);
        tmp = tmp(num:(num+Km-1));
        mtx1 = mtx1(tmp,p);
        p = randperm(n);
        q = randperm(n);
        tmp = [q,q];
        num = find(tmp==1);
        tmp = tmp(num:(num+Kn-1));
        mtx2 = mtx2(tmp,p);
        mtx = kron(mtx1,mtx2);
    case 'RW'
        TTL = 150;
        seqs = zeros(TTL,K);
        mtx = zeros(K,N);
        for i=1:K
            numMat = reshape(1:N,m,n);
            NoSeq = ceil(rand*N);
            [ NoSeq ] = rwalk( numMat,0,NoSeq,TTL,1 );
            seqs(:,i) = NoSeq;
            mtx(i,unique(NoSeq)) = 1;
        end
    otherwise
        mtx = orth(randn(N));
        mtx = mtx(1:K,:);
end

Hand = @(z) mtx*z;
mtxT = mtx';
HandT = @(z) mtxT*z;

end

