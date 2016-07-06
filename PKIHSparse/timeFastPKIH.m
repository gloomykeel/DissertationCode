clear all
close all
res=[];
iter=1;
for N = [5e3,1e4,2e4,5e4,1e5,2e5,5e5,1e6]%];
    tASRM = [];
    tAPKUH = [];
    tASRML = [];
    tASRML5 = [];
    d=32;
    q=37;
    for j=1:10
        p=ceil(N/(q*d));
        eN = p*d*q;
        x = rand(eN,1);        
        tic
        for i=1:100
            H = hadamard(d)/sqrt(d);
            r = randperm(eN);
            x(r) = x;
            y = H*reshape(x,d,[]);
        end
        tSRM = toc;
                
        vec = reshape(1:eN,q,[])';
        vec = vec(:);
        
        tic
        for i=1:100
            H = hadamard(d)/sqrt(d);
%             vec = reshape(1:eN,q,[])';
%             x(vec(:))=x;
            x(vec) = x;
            y = H*reshape(x,d,[]);
        end
        tPKUH = toc;
        
        tic
        for i=1:100
            H = hadamard(d)/sqrt(d);
            r = 2*round(rand(eN,1))-1;
            x = x.*r;
            y = H*reshape(x,d,[]);
        end
        tSRML = toc;
        
        d=512;
        eN = ceil(eN/(d))*(d);
        x = rand(eN,1);
        tic
        for i=1:100
            H = hadamard(d)/sqrt(d);
            r = 2*round(rand(eN,1))-1;
            x = x.*r;
            y = H*reshape(x,d,[]);
        end
        tSRML5 = toc;
        
        tASRM = [tASRM tSRM];
        tASRML = [tASRML tSRML];
        tASRML5 = [tASRML5 tSRML5];
        tAPKUH = [tAPKUH tPKUH];
        [tASRM;tASRML;tASRML5;tAPKUH] 
    end
    tmp.tASRM=tASRM;
    tmp.meanSRM = mean(tASRM');
    tmp.tASRML=tASRML;
    tmp.meanSRML = mean(tASRML');
    tmp.tASRML5=tASRML5;
    tmp.meanSRML5 = mean(tASRML5');
    tmp.tAPKUH=tAPKUH;
    tmp.meanPKUH = mean(tAPKUH');
    res=[res,tmp];
    iter=iter+1;
end
meanSRM=[res.meanSRM];
meanSRML=[res.meanSRML];
meanSRML5=[res.meanSRML5];
meanPKUH=[res.meanPKUH];
figure; plot([meanSRM' meanSRML' meanSRML5' meanPKUH'])

