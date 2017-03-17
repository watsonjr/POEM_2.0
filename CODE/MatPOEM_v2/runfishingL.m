%Run baserun to spinup biomasses
baserun
simname0 = simname;

%Fishing rates
F = linspace(0, 1, 10);

%Lenth of model run (years)
param.tEnd = 50;
B = NaN*ones(length(F),8);
Y = NaN*ones(length(F),1);
for i = 1:length(F)
    %Fish only large fishes
    param.F = param.dt * [0 0   0 0.1*F(i) F(i)     0 0.1*F(i) F(i)]';
    result = poem(param, result);
    results(i) = result;
    Y(i) = sum(result.Y(end,:));
    B(i,:) = result.B(end,:);    
end

%% Save
dirname
save(['/Volumes/GFDL/CSV/Matlab_new_size/',simname0,'/fishing_' harv '.csv'],'results','Y','B','F')

%Plot yield vs. Fishing rate
clf
plot(F,Y,'o-')
xlabel('Fishing rate (yr^-^1)')
ylabel([harv ' Yield'])
print('-dpng',['/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/fishing_' harv '_' simname0])

    