% 计算概率，理论完善前没有特殊用处
s=40;
m=240/s;
n=960/s;
sp=100;
d=ceil((1/2/sp*s+1).^2)
pr = 0;
for k=d:m
    nk = nchoosek(m,k);
    pr = pr + nk.*s.^-k;
end

n*(n-1)/s*pr
d=1/10;
k=10;
n*(n-1)/s*nchoosek(m,ceil(m*d./k))*s.^-ceil(1+m*d./k)

% b=60;
% m=30;
% for i=1:m
%     nk = nchoosek(m,i);
%     pr = nk*b.^-i;
%     hold on
%     plot(i,pr,'*')
% end

fid=fopen('rnd','wb');
for i=1:1024*1024
w=fwrite(rand(1024)>0.5,'uint8')
end
fclose(fid);