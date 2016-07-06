strm = RandStream.getGlobalStream;
if(~strm.Seed)
    %     RandStream.setGlobalStream(RandStream('mt19937ar','Seed','shuffle'));
    RandStream.setGlobalStream(RandStream('mt19937ar','Seed',sum(100*clock)));
    strm = RandStream.getGlobalStream;
end

N = 960;
M = 240;
blkSize=40;[16,24,40,48,60];
for bd=1:length(blkSize)
tM = blkSize(bd); tN = tM;
blkM = M/tM; blkN = N/tN;
TStr = 'T'; % 'T' 'C' 'R' 'G'
blkStr = 'DC'; % 'S' 'SR' 'SC' 'CR' 'TG' 'C' 'T' 'TR' 'TC' 'DC'
recStr = 'CVX'; % GPSR
Sbg = 60;
Stp = 5;
trials = 100;
fileName = sprintf('%sblk%s_N%dM%dtM%dtN%dS%dD%d_%d_%d',...
    TStr,blkStr,N,M,tM,tN,Sbg,Stp,strm.Seed,mod(uint32(sum(clock())),19931));

stats = [];
for k = Sbg:Stp:100;
    T = [];
    stat = zeros(trials,3);
    switch(upper(TStr))
        case 'T'
            T = toeplitz(0:-1:-blkM+1,0:blkN-1);
        case 'C'
            T = toeplitz([0,fliplr(blkN-blkM+1:blkN-1)],0:blkN-1);
        case 'R'
            T = zeros(blkM,blkN);
        case 'G'
            T = [];
    end
    
    for iter = 1:trials        
        BT = num2cell(T);
        switch(upper(TStr))
            case {'C','T'}
                for kl=-blkM+1:blkN-1
                    To = genToepblk( tM,tN,blkStr );
                    ind = find(T==kl);
                    for kj=1:length(ind)
                        BT{ind(kj)}=To;
                    end
                end
            case 'R'
                for kl=1:(blkM*blkN)
                    BT{kl} = genToepblk( tM,tN,blkStr );
                end
            case 'G'
                BT = [];
                mtx = orth(randn(N));
                mtx = mtx(1:M,:);
        end        
        if(iscell(BT)), mtx = cell2mat(BT); end
        
        x = zeros(N,1);
        s = randperm(N);
        x(s(1:k)) = randn(k,1);
        y = mtx*x;
        
        tic;
        if(strcmpi(recStr,'CVX'))
            cvx_begin quiet
            variable alp_debias(N);
            minimize( norm(alp_debias,1) );
            subject to
            norm(mtx*alp_debias - y, 2) <= 1e-3;
            cvx_end
        elseif(strcmpi(recStr,'GPSR'))
            tau = norm(y,'fro')/sqrt(k)/8;%tau = 0.001*max(abs(PhiT(y)));
            [alp,alp_debias,objective,times,debias_start,mses,taus]= ...
                GPSR_BB(y,mtx*dctm,tau,...
                'Debias',1,...
                'Monotone',1,...
                'Initialization',0,...
                'StopCriterion',3,...
                'ToleranceA',1e-3,...
                'ToleranceD',1e-6,...
                'Verbose',0,...
                'MaxiterA',200);
        end
        tmpTime = toc;
        
        if isempty(alp_debias)
            MSE = 1000;
            PSNR = 0;
        else
            [MSE,PSNR,~] = statistics(x,alp_debias);
        end
        NORM = tmpTime;
        stat(iter,:) = [MSE,PSNR,NORM];
        fprintf('k %d, iter %d, MSE %e, PSNR %f, Time %f\n',k ,iter,MSE,PSNR,NORM);
    end
    stats = cat(2,stats,stat);
    save(fileName,'stats');
    if(sum(stat(:,2)>30)<1); break; end
end
s = stats(:,2:3:end);
figure(1);hold on;plot(Sbg:Stp:(Sbg-Stp+Stp*length(s(1,:))),sum(s>30)/trials,'-*')
figure(2);hold on;plot(Sbg:Stp:(Sbg-Stp+Stp*length(s(1,:))),sum(s>30)/trials,'-*')
end
