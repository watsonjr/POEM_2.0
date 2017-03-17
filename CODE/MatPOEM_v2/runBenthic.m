%
% Sweep over benthic production to test carrying capacity
%
baserun
result0 = result;

factor = linspace(2.0,0.1,20);
param.tEnd = 50;

%%
B = NaN*ones(length(factor),13);
for i = 1:length(factor)
    param.K(3:5) = factor(i);
    result = poem(param, result);
    results(i) = result;
    B(i,:) = result.y(end,:);
end
%%

save(['/Volumes/GFDL/CSV/Matlab_new_size/',simname,'/bentrun.mat'],'results','factor','B')

%Plot
ix = cell([3,1]);
FF = 1;
LP = 2;
LD = 3;

ix{FF} = 6:7;
ix{LP} = 8:10;
ix{LD} = 11:13;
col = cell([3,1]);
col{FF} = 'r-';
col{LP} = 'b-';
col{LD} = 'k-';

clf
%subplot(2,2,1)
for i = 1:3
    plot(factor, log10(B(:,ix{i}(end))), col{i},'LineWidth',2)
    hold on
end
legend('MF','LP','LD')
xlabel('Benthic K')
ylabel('log10 Biomass')
print('-dpng',['/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/bentrun_' simname])


