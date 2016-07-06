function y = Qsnr(x,xr);
x = x(:);
xr = xr(:);
y = 20*log10(norm(x)/norm(x-xr));
