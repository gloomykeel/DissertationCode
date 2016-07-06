% This function implements fast spatial-domain compressive sampling of 2D
% images through structurally random matrices (SRM). Details of SRM can be
% found in the following paper:
function [psnr_val]=pkihGPSR(x, samRate, trans_mode, rand_type, blkSz, q, o)

% Test the input parameters:
if nargin<4,
    error('Wrong Number of Input Parameters');
end

% keep an original copy of the input signal
x0 = x;
[m,n] = size(x);
% Total number of pixels;
N = m*n;
K = round(N.*samRate);
x = x(:);

% Sparsifying transform: 9-7 Wavelet
[h0,h1,f0,f1] = filter9_7();
L = floor(log2(m))-3;               % Level of decomposition

% Initialize the output parameters;
psnr_val=[];
% Main loop to test the reconstruction performance for each K(i);
for i =1:length(K)
    Ki=round(K(i));
    % Define selected samples
    svec=randperm(N);
    svec = svec(1:Ki)';
    % Define the random vector for the input samples:
    switch(rand_type)
        case 0
            % Get the wavelet transform coefficients:
            x = dwt2d(x0,h0,h1,L);
            x = x(:);
            % Define selected samples
            svec = randperm(round(N/2)-1)+1;
            svec = svec(1:round(Ki/2))';
            % Define Sampling Operator;
            B = @(z) fft1d_f(z, svec, 1:N);
            % Define the transpose of the Sampling Operator;
            Bt = @(z) fft1d_t(z, N, svec, 1:N);
            Phi = B;
        case 1
            rvec = randperm(N); % Random permutation;
            % Define Sampling Operator;
            Phi = @(z) blk_f1d(z, svec, rvec, trans_mode, blkSz);
            % Define the transpose of the Sampling Operator;
            Phi_T = @(z) blk_t1d(z, N, svec, rvec,trans_mode, blkSz);
        case 2
            rvec = 2*round(rand(N,1))-1; % Random filipping the sign;
            % Define Sampling Operator;
            Phi = @(z) blk_f1d(z, svec, rvec, trans_mode, blkSz);
            % Define the transpose of the Sampling Operator;
            Phi_T = @(z) blk_t1d(z, N, svec, rvec,trans_mode, blkSz);
            
        case 3
            rvec = reshape(1:ceil(N/(q*blkSz))*(q*blkSz),q,[])';
            rvec = reshape(rvec,blkSz,[]);
            rvec = rvec(:);
            % ql 11.12
            eN=ceil(length(rvec)/o)*o;
            rv = reshape(1:eN,o,[])';
            rv = rv(rv<=length(rvec));
            rvec = rvec(rv);
            %~ql
            % Define Sampling Operator;
            Phi = @(z) blk_frr(z, svec, rvec, trans_mode, blkSz);
            % Define the transpose of the Sampling Operator;
            Phi_T = @(z) blk_trr(z, N, svec, rvec,trans_mode, blkSz);
        case 4
            eN = ceil(N/(q*blkSz))*(q*blkSz);
            rvec = reshape(zeros(eN,1),blkSz,[]);
            mvec=double(mseq(dec2binvec(8219),dec2binvec(8219),lcm(eN*3,blkSz)));
            mvec = reshape(mvec,[],blkSz);
            ind = find(mvec);
            mvec(ind(1:eN))=1:eN;
            for l=1:blkSz
                tmp = mvec(:,l);
                tmp = tmp(tmp>0);
                rvec(l,1:length(tmp))=tmp;
            end
    end
    
    if(rand_type)
        % Define the equivalent sampling operator for the wavelet coefficients;
        B = @(z) A_idwt1d(Phi,z,f0,f1,L,m,n);
        Bt = @(z) At_dwt1d(Phi_T,z,h0,h1,L,m,n);
    end
    
    % getting measurements
    y = Phi(x);
    
    % Reconstruction. Using GPSR_BB modules
    tau = norm(y,'fro')/sqrt(Ki)/9;
    [alp,alp_debias,objective,times,debias_start,mses]= ...
        GPSR_BB(y,B,tau,...
        'AT', Bt,'Debias',1,'Initialization',0,...
        'StopCriterion',1,'ToleranceA',0.0001,'ToleranceD',0.00005,'Verbose',1);
    
    % Transform from the WT domain to the spatial domain
    alp_debias = reshape(alp_debias,m,n);
    xr = idwt2d(alp_debias,f0,f1,L);
    
    psnr_val  = [psnr_val,Qsnr(x0,xr)];
end


