function [ naLen,pathLen ] = pathStat(  )
%PATHSTAT 此处显示有关此函数的摘要
%   此处显示详细说明
n = 8;
N = n*n;
blk = 4;
M = N/2;
A = reshape(randperm(N),[],blk);
pathLen=0;
naLen=0;
for i=1:M
    numMat = reshape(1:N,n,n);
    ind = A(ceil(rand*N/blk),:);
    curi = mod(ind,n);
    curj = ceil(ind/n);
    tmp = sqrt(curi.*curi+curj.*curj);
    iij = find(tmp==min(tmp));
    pathLen = pathLen + curi(iij(1))+curj(iij(1));
    naLen = naLen + curi(iij(1))+curj(iij(1)) + 1;
    curi = mod(ind(1),n);if(curi==0);curi=n;end;curj = ceil(ind(1)/n);
    walkSeq = [];
    while ~isempty(ind)
        walkSeq = [walkSeq;numMat(curi,curj)];
        dist = inf;
        ind(ind == numMat(curi,curj)) = [];
        for j=1:length(ind)
            rows = mod(ind(j),n);if(rows==0);rows=n;end;cols = ceil(ind(j)/n);
%             [rows,cols] = find(numMat==ind(j)); % find position
            tmp = sqrt((curi-rows).^2+(curj-cols).^2);
            if tmp < dist
                dist = tmp;
                curi = rows;
                curj = cols;
            end
        end
    end
    for j=1:(length(walkSeq)-1)
%         [curi,curj] = find(numMat==walkSeq(j));
%         [rows,cols] = find(numMat==walkSeq(j+1));
        curi = mod(walkSeq(j),n);if(curi==0);curi=n;end;curj = ceil(walkSeq(j)/n);
        rows = mod(walkSeq(j+1),n);if(rows==0);rows=n;end;cols = ceil(walkSeq(j+1)/n);
        pathLen = pathLen + abs(curi-rows)+abs(curj-cols);
    end
end

end

