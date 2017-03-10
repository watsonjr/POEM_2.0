%
% Sweep over benthic production
%
baserun
result0 = result;

factor = linspace(0.9,0,10);
param.tEnd = 50;

for i = 1:length(factor)
    param.K(3:5) = 5*factor(i);
    result = poem(param, result);
    B(i,:) = result.y(end,:);
end
%%
ix = cell([3,1]);
FF = 1;
LP = 2;
LD = 3;

ix{FF} = 6:7;
ix{LP} = 8:10;
ix{LD} = 11:13;
col = cell([3,1]);
col{FF} = 'b-';
col{LP} = 'm-';
col{LD} = 'g-';

clf
for i = 1:3
    plot(factor, B(:,ix{i}(end)), col{i})
    hold on
end
