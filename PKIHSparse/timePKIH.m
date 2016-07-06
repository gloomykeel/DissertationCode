clear all
close all
res=[];
iter=1;
for N = 500:500:5000%[5e2,1e3,2e3,5e3]%,1e4];
    tASRM = [];
    tASRML = [];
    tAPKUH = [];
    tAPKUH4 = [];
    for j=1:10        
        d = 4;
        q = 13;
        p = ceil(N/(q*d));
        eN = p*d*q;
        tic
        for i=1:100
            Mtx = kron(sparse(eye(p*q)),hadamard(d));
            r = randperm(eN);
            Mtx(:,r) = Mtx;
        end
        tSRM = toc;
        
        tic
        for i=1:100
            Mtx = kron(sparse(eye(p*q)),hadamard(d));
            r= sparse(1:eN,1:eN,2*round(rand(eN,1))-1);
            Mtx = Mtx*r;
        end
        tSRML = toc;
        
        tic
        for i=1:100
            mx = kron(eye(p),kron(hadamard(d),sparse(eye(q)))); % kron(I2,kron(H,circshift(I1,3)))
        end
        tPKUH = toc;
        
        q = 44; p = ceil(N/(q*d));
        eN = p*d*q;
        tic
        for i=1:100
            mx = kron(eye(p),kron(hadamard(d),sparse(eye(q)))); % kron(I2,kron(H,circshift(I1,3)))
        end
        tPKUH4 = toc;
        
        tASRM = [tASRM tSRM];
        tASRML = [tASRML tSRML];
        tAPKUH = [tAPKUH tPKUH];
        tAPKUH4 = [tAPKUH4 tPKUH4];
        [tASRM;tASRML;tAPKUH;tAPKUH4]
    end
    tmp.tASRM=tASRM;
    tmp.meanSRM = mean(tASRM');
    tmp.tASRML=tASRML;
    tmp.meanSRML = mean(tASRML');
    tmp.tAPKUH=tAPKUH;
    tmp.meanPKUH = mean(tAPKUH');
    tmp.tAPKUH4=tAPKUH4;
    tmp.meanPKUH4 = mean(tAPKUH4');
    res=[res,tmp];
    iter=iter+1;
end
meanSRM=[res.meanSRM];
meanSRML=[res.meanSRML];
meanPKUH=[res.meanPKUH];
meanPKUH4=[res.meanPKUH4];
figure; plot([meanSRM' meanSRML' meanPKUH' meanPKUH4'])

