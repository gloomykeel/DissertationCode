
strm = RandStream.getGlobalStream;
if(~strm.Seed)
%     RandStream.setGlobalStream(RandStream('mt19937ar','Seed','shuffle'));
    RandStream.setGlobalStream(RandStream('mt19937ar','Seed',sum(100*clock)));
    strm = RandStream.getGlobalStream;
end

% set options
m = 32; n = 64;
N = m*n;
K = round(N*0.25);
blkSize = 4;
mtxName = 'TOEP';
mtxName = 'GAU';
mtxName = 'SRMG';
mtxName = 'PKUH';
mtxName = 'KSRMG';
mtxName = 'KGAU';
mtxNames = {'PKUH';'KSRMG';'KGAU';};
trials = 200;
threshDB = 30;
for Nind=1:length(mtxNames)
    mtxName = mtxNames{Nind};
% if (strcmp(mtxName,'PKUH')); fctr = 2; else fctr = 1;end
for fctr=1:2
fileName = sprintf('%s8_%d_%d_%d',mtxName,fctr,strm.Seed,mod(uint32(sum(clock())),19931));
Z = 0; stats = [];
for k = 60:5:200
% for k = 120
    stat = zeros(trials,3);
    parfor iter =1:trials
        s = zeros(N,1);
        x = randperm(N);
        s(x(1:k)) = 2*round(rand(k,1))-1;
        
        % generate sensing matrix (handle)
        [Phi,PhiT] = getMtxHandle(mtxName,N,K,blkSize,m*fctr,n/fctr);
        % get sparse basis
        [PB,PBt,B,Bt] = getBasisHandle('DCT',m,n,Phi,PhiT);
        
        x = B(s);
        % do sensing
        y = Phi(x);
        % do recover
        % Reconstruction. Using GPSR_BB modules
        tau = norm(y,'fro')/sqrt(N)/8;%tau = 0.001*max(abs(PhiT(y)));
        % [alp,alp_debias0,objective,times,debias_start,mses]= ...
        %      GPSR_B(y,B,tau,...
        %      'AT', Bt,'Debias',1,'Initialization',0,'StopCriterion',1,......
        %      'ToleranceA',1e-4,'MaxiterA',200);
%         cvx_begin
%         variable alp_debias(N);
%         minimize( norm(alp_debias,1) );
%         subject to
%         norm(PB(alp_debias) - y, 2) <= 1e-3;
%         cvx_end
        [alp,alp_debias,objective,times,debias_start,mses,taus]= ...
            GPSR_BB(y,PB,tau,'AT', PBt,...
            'Debias',1,...
            'Monotone',1,...
            'Initialization',0,...
            'StopCriterion',3,...
            'ToleranceA',1e-3,...
            'ToleranceD',1e-6,...
            'Verbose',0,...
            'MaxiterA',200);
        if isempty(alp_debias)
            MSE = 10;
            PSNR = 0;
            NORM = 10;
        else
        % Transform from the Basis domain to the spatial domain
        xr = Bt(alp_debias');
        % statistics
        [MSE,PSNR,NORM] = statistics(x,xr);
        end
        stat(iter,:) = [MSE,PSNR,NORM];
        fprintf('iter %d-%d, MSE %e, PSNR %f, NORM %f\n',k,iter,MSE,PSNR,NORM);
    end
    stats = [stats,stat];
    save(fileName,'stats');
    if sum(stat(:,2)>threshDB)<1
        Z = Z+1;
        if Z>1; break; end
    else
        Z = 0;
    end
end
end
end
tmp = stats(:,2:3:end);
plot(sum(tmp>threshDB),'-*')