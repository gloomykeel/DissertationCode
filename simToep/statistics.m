function [ MSE,PSNR,NORM ] = statistics( x,xr )
%STATISTICS 计算信号及恢复信号的统计值
%   此处显示详细说明
x = double(x(:));
xr = double(xr(:));
NORMX = norm(x);
NORM = norm(x-xr);
MSE = NORM*NORM/(NORMX*NORMX);
PSNR = 20*log10(NORMX/NORM);

end

