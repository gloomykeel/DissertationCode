function [ To ] = genToepblk( tM,tN,typeStr )
%GENTOEPBLK 按照输入生成对应的托普利兹块矩阵
%   此处显示详细说明
To=[];
value = 1;
rVec = zeros(tN,1); cVec = zeros(tM,1);
trp = randperm(tN); tcp = randperm(tM);
%             if value==1; value=1; else value=-1; end
if(trp(1)==1); trp(1)=trp(2); end
switch(upper(typeStr))
    case 'S'
        rVec(trp(1)) = value;
        cVec(trp(1)) = -value;
        cVec(1) = rVec(1);
        To = toeplitz(cVec,rVec);
    case 'SR'
        rVec(trp(1)) = value;
        cVec(tcp(1)) = -value;
        cVec(1) = rVec(1);
        To = toeplitz(cVec,rVec);
    case 'SC'
        rVec(trp(1)) = value;
        cVec(tcp(1)) = -value;
        cVec(1) = rVec(1);
        To = circshift(toeplitz(cVec,rVec),ceil(rand*tN),2);
    case 'CR'
        rVec(trp(1)) = (rand(1,1)<0.5)*2-1;
        cVec(tcp(1)) = (rand(1,1)<0.5)*2-1;
        cVec(1) = rVec(1);
        To = circshift(toeplitz(cVec,rVec),ceil(rand*tN),2);
    case 'TG'
        To = orth(randn(max(tM,tN)));
        To = To(1:tM,1:tN);
    case 'C'
        rVec(trp(1)) = value;
        rVec(trp(2)) = -value;
        cVec = flipud(rVec);
        cVec = [rVec(1);cVec(2:end)];
        To = toeplitz(cVec,rVec);
    case 'T'
        rVec(trp(1)) = value;
        cVec(tN-trp(1)+2) = -value;
        To = toeplitz(cVec,rVec);
    case 'TR'
        rVec(trp(1)) = (rand(1,1)<0.5)*2-1;
        cVec(tN-trp(1)+2) = (rand(1,1)<0.5)*2-1;
        To = toeplitz(cVec,rVec);
    case 'TC'
        rVec(trp(1)) = value;
        cVec(tN-trp(1)+2) = -value;
        To = circshift(toeplitz(cVec,rVec),ceil(rand*tN),2);
    case 'DC'
        cVec = ones(tM,1);
        cVec(2:2:end)=-1;
        cVec = diag(cVec);
        To = circshift(cVec,ceil(rand*tN),2);
end

end