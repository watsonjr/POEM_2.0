%% POEM Make file

clear all
close all

%%%%!! EXPERIMENTS
testoneloc = false;
testlocs = false;
testonelocP = false;
oneloc_fishing = false;
oneloc_hind_pristine = false;
oneloc_fore_pristine = false;
spinup_pristine = true;
climatol = false;
pre_industrial = false;
historic_pristine = false;
historic_fished = false;
forecast_pristine = false;
forecast_fished = false;

tic
if testoneloc
    Testoneloc()
end
if testlocs
    Testlocs()
end
if testonelocP
    TestonelocP()
end
if spinup_pristine
    Spinup_pristine()
end
if climatol
    Climatol_pristine()
end
if pre_industrial
    Pre_industrial()
end
if historic_pristine
    Historic_pristine()
end
if historic_fished
    Historic_fished()
end
if forecast_pristine
    Forecast_pristine()
end
if forecast_fished
    Forecast_fished()
end
toc
