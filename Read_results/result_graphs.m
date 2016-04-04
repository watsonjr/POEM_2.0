clear all
%read in datafile from POEM2 COBALTS
%PISC_data = csvread('/volumes/bigDisk/POEM2/CODE/Data_hindcast/CSV2/Data3_PISC.csv');
%plot(PISC_data(:, 10))

%read in datafile from POEM2 TP
subplot(2, 2, 1)
PISC_data = csvread('/volumes/bigDisk/POEM2/CODE/Data_hindcast/CSV2/DataTP_PISC.csv');
plot(PISC_data(:, 1:10))
title('piscivores')
ylabel('g WW m-2');
xlabel('time (days)')
legend('bin1', 'bin2', 'bin3', 'bin4', 'bin5', 'bin6', 'bin7', 'bin8', 'bin9', 'bin10');

subplot(2, 2, 2)
PLAN_data = csvread('/volumes/bigDisk/POEM2/CODE/Data_hindcast/CSV2/DataTP_PLAN.csv');
plot(PLAN_data(:, 1:15))
title('planktivores')
ylabel('g WW m-2');
xlabel('time (days)')
legend('bin1', 'bin2', 'bin3', 'bin4', 'bin5', 'bin6', 'bin7', 'bin8', 'bin9', 'bin10', 'bin11', 'bin12', 'bin13', 'bin14', 'bin15');

subplot(2, 2, 3)
DETR_data = csvread('/volumes/bigDisk/POEM2/CODE/Data_hindcast/CSV2/DataTP_DETR.csv');
plot(DETR_data(:, 1:7))
title('benthos')
ylabel('g WW m-2');
xlabel('time (days)')
legend('bin1', 'bin2', 'bin3', 'bin4', 'bin5', 'bin6', 'bin7');

subplot(2, 2, 4)
W_data = csvread('/volumes/bigDisk/POEM2/CODE/Data_hindcast/CSV2/DataTP_W.csv');
plot(W_data)
title('detritus')
ylabel('g WW m-2');
xlabel('time (days)')
legend('bin1');