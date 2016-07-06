function [ MSE,PSNR,NORM ] = statistics( x,xr )
%STATISTICS �����źż��ָ��źŵ�ͳ��ֵ
%   �˴���ʾ��ϸ˵��
x = double(x(:));
xr = double(xr(:));
NORMX = norm(x);
NORM = norm(x-xr);
MSE = NORM*NORM/(NORMX*NORMX);
PSNR = 20*log10(NORMX/NORM);

end

