%Fishing rates
F = linspace(0, 1, 10);

%Lenth of model run (years)
param.tEnd = 50;
Bf = NaN*ones(length(F),8);
Y = NaN*ones(length(F),1);
for i = 1:length(F)
    %Fish F
    param.F = [0.1*F(i) F(i)  0 0 0               0 0 0]';
    result = poem(param, result);
    fishing(i) = result;
    Y(i) = sum(result.Y(end,:));
    Bf(i,:) = result.B(end,:);    
end

%% Save
dirname
save(['/Volumes/GFDL/CSV/Matlab_new_size/',simname,'/fishing_' harv '.mat'],'fishing','Y','Bf','F')

%Plot yield vs. Fishing rate
clf
plot(F,Y,'bo-')
xlabel('Fishing rate (yr^-^1)')
ylabel([harv ' Yield'])
title(['Repro effic = ' num2str(param.RE)])
print('-dpng',['/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/fishing_' harv '_' simname])

    