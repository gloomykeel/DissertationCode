function [ PB,PBt,B,Bt ] = getBasisHandle( name,m,n,Phi,Phi_T )
%GETBASISHANDLE 得到变换基底的句柄
%   此处显示详细说明

name = upper(name);

switch(name)
    case 'DCT'
        mtxm = dctmtx(m);
        mtxn = dctmtx(n);
        mtx = kron(mtxn,mtxm);
        B = @(z) mtx*z;
        PB = @(z) Phi(mtx*z);
        mtx = mtx';
        Bt = @(z) z*mtx;
        PBt = @(z) mtx*Phi_T(z);
    case 'EYE'
        B = @(z) z;
        Bt = @(z) z;
        PB = @(z) Phi(z);
        PBt = @(z) Phi_T(z);
    case 'DWT'
        % Sparsifying transform: 9-7 Wavelet
        [h0,h1,f0,f1] = filter9_7();
        L = floor(log2(m))-3;               % Level of decomposition
        
        B = @(z) dwt2d(reshape(z,m,n),f0,f1,L);
        Bt = @(z) idwt2d(reshape(z,m,n),f0,f1,L);
        % Define the equivalent sampling operator for the wavelet coefficients;
        PB = @(z) A_idwt1d(Phi,z,f0,f1,L,m,n);
        PBt = @(z) At_dwt1d(Phi_T,z,h0,h1,L,m,n);
    otherwise
end

end

