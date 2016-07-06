
s = RandStream.getGlobalStream;
if(~s.Seed)
%     RandStream.setGlobalStream(RandStream('mt19937ar','Seed','shuffle'));
    RandStream.setGlobalStream(RandStream('mt19937ar','Seed',sum(100*clock)));
    s = RandStream.getGlobalStream;
end

msg.Seed = s.Seed;
fileName = sprintf('MsgRes%d_%s.mat',s.Seed,datestr(datetime,'HH_MM_SS_FFF'));
% parameters
N = 512;
msg.N = N; % signal length
k = 30;
msg.k = k;% sparsity
% sparsifying transform
Psi = dctmtx(N);
msg.blkSz=4;
trial_num=200;
msg.trial_num = trial_num;
msg.start_samp_num = 70;

msg.srmg = []; msg.srml = []; msg.pkih = []; msg.gau = [];
for i=1:100
    srmgs = []; srmls = []; pkihs = []; gaus = [];
for j =1:10
    samp_num = msg.start_samp_num+(j-1)*10;
    message=sprintf('Iter %d Sampling No.=%d',i,samp_num);
    disp(message);
    srmg = ompEval(N,k,samp_num,msg.blkSz,1,1,Psi,trial_num)/trial_num;
    srmgs = [srmgs;srmg];
    srml = ompEval(N,k,samp_num,msg.blkSz,1,2,Psi,trial_num)/trial_num;
    srmls = [srmls;srml];
    pkih = ompEval(N,k,samp_num,msg.blkSz,5,3,Psi,trial_num,3)/trial_num;
    pkihs = [pkihs;pkih];
    gau = ompEval(N,k,samp_num,msg.blkSz,1,4,Psi,trial_num)/trial_num;
    gaus = [gaus;gau];
    disp([srmgs,srmls,pkihs,gaus]);
end
    msg.srmg = [msg.srmg,srmgs];
    msg.srml = [msg.srml,srmls];
    msg.pkih = [msg.pkih,pkihs];
    msg.gau = [msg.gau,gaus];
    save(fileName,'msg');
end
    
    