function [ xc,yc ] = getFigData( fnam )
%GETFIGDATA �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
if(isnumeric(fnam))
    figure(fnam);
else
    open(fnam); % fnam���ļ���
end
lh = findall(gca,'type','line'); % �ӵ�ǰͼ(gca)��ȡ�����ߵ�handle
xc = get(lh,'xdata'); % ȡ��x�����ݣ�ע�⣬���x��y����cell�����ݽṹ�����
yc = get(lh,'ydata'); % ȡ��y������
% for i=1:17;xc{i}=(40+(xc{i}-1)*5);end
% for i=1:17;set(lh(i),'xdata',xc{i});end
end

