%! Create a directory for output
tfcrit = num2str(int64(100*param.fc));
tcfn = num2str(param.h);
tefn = num2str(param.gamma);
td = num2str(1000+int64(100*param.D));
tj = num2str(1000+int64(100*param.J));
tsm = num2str(1000+int64(100*param.Sm));
%tbe = num2str(100+int64(100*param.BE));
tmort = num2str(1000+int64(100*mean(param.nmrt)));
tcc = num2str(1000+int64(100*mean(param.K(3:5))));
tre = num2str(100000+int64(round(100000*param.RE)));
if(exist('frate','var'))
    tfish = num2str(1000+int64(100*frate));
end
if (param.F(2) > 0)
    if (param.F(5) > 0 && param.F(8) > 0)
        harv='All';
    else
        harv='F';
    end
else
    if (param.F(5) > 0 && param.F(8) > 0)
        harv = 'L';
    elseif (param.F(5) > 0)
        harv = 'P';
    elseif (param.F(8) > 0)
        harv = 'D';
    else
        harv=[];
    end
end
simname = ['Enc',tefn,'_cmax-met',tcfn,'_fc',tfcrit,'_nmort',tmort(2:end),'_D',td(2:end),'_J',tj(2:end),'_Sm',tsm(2:end),'_RE',tre(2:end),'_CC',tcc(2:end)];
if (~isdir(['/Volumes/GFDL/CSV/Matlab_new_size/',simname]))
    mkdir(['/Volumes/GFDL/CSV/Matlab_new_size/',simname])
end
if (~isdir(['/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/',simname]))
    mkdir(['/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/',simname])
end