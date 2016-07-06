a={};b={};
load('snrs212934.mat')
a=snrs212934;
load('snrs210755.mat')
b=snrs210755;
a=cat(2,a,b);
[m,n] = size(a);
for i=1:(m*n);a{i}=mean(a{i},3); end
c=zeros(9,5,m);for i=1:m;for j=1:13; c(:,:,i)=c(:,:,i)+a{i,j}; end;end;c=c./j;
x=0.1:0.05:0.5;
for i=1:4
figure(i);
plot(x,c(:,1,i),'-',...
    x,c(:,2,i),'-o',...
    x,c(:,3,i),'-s',...
    x,c(:,4,i),'-v',...
    x,c(:,5,i),'-*')
legend('WPFFT','SRMG','SRML','SRML512','PKIH(39,37)','Location','northwest')
xlabel('Sampling Rate')
ylabel('PSNR(dB)')
end
