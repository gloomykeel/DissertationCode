% psnr calc of two images 
% motified by ql  
function y = psnr(im1,im2)
[m,n] = size(im1);
mse = norm(double(im1(:)) - double(im2(:)));
mse = (mse*mse)/(m*n);
if mse > 0
    y = 10*log10(255^2/mse);
else
    disp('infinite psnr');
end