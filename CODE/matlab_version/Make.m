%% POEM Make file

clear all
close all

%%%%!! EXPERIMENTS
testoneloc = false;
testlocs = true;
testonelocP = false;
oneloc_fishing = false;
oneloc_hind_pristine = false;
oneloc_fore_pristine = false;
spinup_pristine = false;
pre_industrial = false;
historic_pristine = false;
historic_fished = false;
forecast_pristine = false;
forecast_fished = false;

if testoneloc
    tic
    Testoneloc()
    toc
end
if testlocs
    tic
    Testlocs()
    toc
end
if testonelocP
    tic
    TestonelocP()
    toc
end
if spinup_pristine
    tic
    Spinup_pristine()
    toc
end
if pre_industrial
    tic
    Pre_industrial()
    toc
end
if historic_pristine
    tic
    Historic_pristine()
    toc
end
if historic_fished
    tic
    Historic_fished()
    toc
end
if forecast_pristine
    tic
    Forecast_pristine()
    toc
end
if forecast_fished
    tic
    Forecast_fished()
    toc
end
