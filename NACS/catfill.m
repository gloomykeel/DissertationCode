function [ y ] = catfill( a,b )
%CATFILL �����������������������Ӻϣ�������
%   �˴���ʾ��ϸ˵��
[~,na] = size(a);
[mb,nb] = size(b);
x = zeros(mb,na);
x(1:mb,1:nb) = b;
y = [a;x];
end

