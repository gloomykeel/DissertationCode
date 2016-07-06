function [ y ] = catfill( a,b )
%CATFILL 根据输入参数，将矩阵纵向接合，横向补零
%   此处显示详细说明
[~,na] = size(a);
[mb,nb] = size(b);
x = zeros(mb,na);
x(1:mb,1:nb) = b;
y = [a;x];
end

