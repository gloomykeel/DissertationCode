function [ xc,yc ] = getFigData( fnam )
%GETFIGDATA 此处显示有关此函数的摘要
%   此处显示详细说明
if(isnumeric(fnam))
    figure(fnam);
else
    open(fnam); % fnam是文件名
end
lh = findall(gca,'type','line'); % 从当前图(gca)中取出曲线的handle
xc = get(lh,'xdata'); % 取出x轴数据，注意，这个x和y是以cell的数据结构保存的
yc = get(lh,'ydata'); % 取出y轴数据
% for i=1:17;xc{i}=(40+(xc{i}-1)*5);end
% for i=1:17;set(lh(i),'xdata',xc{i});end
end

