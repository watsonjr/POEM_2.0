%
% Sweep over benthic production
%
baserun
result0 = result;

param.F = [0.03 0.3 0 0.03 0.3 0 0.03 0.3]';

factor = linspace(0, 2, 10);
param.tEnd = 50;

for i = 2:length(factor)
    param.r(3:5) = 10*factor(i);
    param.K(3:5) =  1*factor(i);
    result(i) = poem(param, result(i-1));
    B(i,:) = result(i).y(end,:);
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
