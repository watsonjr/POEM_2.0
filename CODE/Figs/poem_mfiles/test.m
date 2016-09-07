%Consolodate spinup sim output into tables

clear all
close all

datap = '/Volumes/GFDL/CSV/';
d = dir(datap);
d2={d.name};
d2=d2';
d3=d2(4:163);