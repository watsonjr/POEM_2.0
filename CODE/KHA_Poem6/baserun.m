%
% run all 3 types
%

%Set parameters
baseparameters

%Fishing rate
param.F = 0.0*[0 0 0 0 0 0 0.1 1]';
%Length of run (years)
param.tEnd = 200;
%Simname & directory
dirname
%Results
result = poem(param);
save(['/Volumes/GFDL/CSV/Matlab_new_size/',simname,'/baserun.mat'],'result')
%Plot
plotPoem(param, result)
print('-dpng',['/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/baserun_' simname])

