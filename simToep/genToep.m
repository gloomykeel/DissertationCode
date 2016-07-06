N = 960;
M = 240;
tM = 16;
tN = 16;
blkM = M/tM;
blkN = N/tN;
typeStr = 'ZH';
value = 1;
stats = [];
for k = 60:10:100;
    trials = 100;
    iter = 1;
    stat = zeros(trials,3);
    for iter = 1:trials
        T = cell(blkM,blkN);
        for i=1:(blkM*blkN)
            rVec = zeros(tN,1);
            cVec = zeros(tM,1);
            trp = randperm(tN);
            tcp = randperm(tM);
            %             if value==1; value=1; else value=-1; end
            switch(upper(typeStr))
                case 'DC'
                    if(trp(1)==1); trp(1)=trp(2); end
                    rVec(trp(1)) = value;
                    cVec(trp(1)) = -value;
                    cVec(1) = rVec(1);
                    T{i} = toeplitz(cVec,rVec);
                case 'DCSJ'
                    rVec(trp(1)) = value;
                    cVec(tcp(1)) = -value;
                    cVec(1) = rVec(1);
                    T{i} = toeplitz(cVec,rVec);
                case 'DCZH'
                    rVec(trp(1)) = value;
                    cVec(tcp(1)) = -value;
                    cVec(1) = rVec(1);
                    T{i} = circshift(toeplitz(cVec,rVec),ceil(rand*tN),2);
                case 'ZH'
                    rVec(trp(1)) = (rand(1,1)<0.5)*2-1;
                    cVec(tcp(1)) = (rand(1,1)<0.5)*2-1;
                    cVec(1) = rVec(1);
                    T{i} = circshift(toeplitz(cVec,rVec),ceil(rand*tN),2);
                case 'GAU'
                    T{i} = zeros(tM,tN);
            end
        end
        if(strcmpi(typeStr,'GAU'))
            mtx = orth(randn(N));
            mtx = mtx(1:M,:);
        else
            mtx = cell2mat(T);
        end
        
        T = toeplitz(0:-1:-blkM+1,0:blkN-1);
        CT = num2cell(T);
        for kl=-blkM+1:blkN-1
            rVec = zeros(tN,1);
            cVec = zeros(tM,1);
            trp = randperm(tN);
            tcp = randperm(tM);
            rVec(trp(1)) = value;
            cVec(tcp(1)) = -value;
            cVec(1) = rVec(1);
            ind = find(T==kl);
            to=randn(tM,tN);toeplitz(cVec,rVec);
            for kj=1:length(ind)
                CT{ind(kj)}=to;
            end
        end
        mtx = cell2mat(CT);
        
        x = zeros(N,1);
        s = randperm(N);
        x(s(1:k)) = 1;
        y = mtx*x;
        
        tau = norm(y,'fro')/sqrt(k)/8;%tau = 0.001*max(abs(PhiT(y)));
        [alp,alp_debias,objective,times,debias_start,mses,taus]= ...
            GPSR_BB(y,mtx,tau,...
            'Debias',1,...
            'Monotone',1,...
            'Initialization',0,...
            'StopCriterion',3,...
            'ToleranceA',1e-3,...
            'ToleranceD',1e-6,...
            'Verbose',0,...
            'MaxiterA',200);
        
        if isempty(alp_debias)
            MSE = 1000;
            PSNR = 0;
            NORM = 1000;
        else
            [MSE,PSNR,NORM] = statistics(x,alp_debias);
        end
        stat(iter,:) = [MSE,PSNR,NORM];
        fprintf('iter %d, MSE %e, PSNR %f, NORM %f\n',iter,MSE,PSNR,NORM);
    end
    stats = [stats,stat];
end
cVec = stats(:,2:3:end);
figure(1);hold on;plot(sum(cVec>30)/trials,'-*')
