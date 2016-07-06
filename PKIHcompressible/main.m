
path(path,genpath('./wavelet_routines'))
s = RandStream.getGlobalStream;
if(~s.Seed)
    %     RandStream.setGlobalStream(RandStream('mt19937ar','Seed','shuffle'));
    RandStream.setGlobalStream(RandStream('mt19937ar','Seed',sum(100*clock)));
    s = RandStream.getGlobalStream;
end

imgs = {'lena256.bmp';'camera.png';'peppers256.png';'boat256.bmp'};%;'dscgray.jpg'};
%imgs = {'lena256.bmp'};
samRate=0.1:0.05:0.5;
snrs = {};
for i = 1:length(imgs)
    % Read the input image;
    disp(imgs{i})
    img = double(imread(imgs{i}));
    psnrs=[];
    trials = 10;
    wpft = zeros(length(samRate),trials);srmg=wpft;srml=wpft;srml5=wpft;pkih=wpft;
    fileName = sprintf('MsgRes%d_%s.mat',s.Seed,datestr(clock,'HH_MM_SS_FFF'));
    for j=1:trials
        disp(j);tic;
        % Benchmark: Wavelet domain measurement through partial FFT
        wpft(:,j) = pkihGPSR(img, samRate, 'WPFFT', 0);
        psnrs = cat(3,psnrs,wpft);
        %     [~,psnr_bf]=pfft_cswt2d(imgs{i}, samRate);
        srmg(:,j) = pkihGPSR(img, samRate, 'BWHT', 1, 32);
        psnrs = cat(3,psnrs,srmg);
        srml(:,j) = pkihGPSR(img, samRate, 'BWHT', 2, 32);
        psnrs = cat(3,psnrs,srml);
        srml5(:,j) = pkihGPSR(img, samRate, 'BWHT', 2, 512);
        psnrs = cat(3,psnrs,srml5);
        pkih(:,j) = pkihGPSR(img, samRate, 'BWHT', 3, 32, 39, 37);
        psnrs = cat(3,psnrs,pkih);
        save(fileName,'psnrs');
        disp([wpft(:,j),srmg(:,j),srml(:,j),srml5(:,j),pkih(:,j)]);toc;
    end
    snrs = cat(1,snrs,psnrs);
    fileName = sprintf('MsgSnr%d_%s.mat',s.Seed,datestr(clock,'HH_MM_SS_FFF'));
    save(fileName,'snrs');
end

