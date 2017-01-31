%%% Biomass recruiting to size-class (g m-2 d-1)
function rec = sub_rec_larv(X,bio,rfrac)
  %X: repro rate
  %rfrac: repro efficiency ~ sex ratio * egg survival
  
  %! Global constant RE
  rec = rfrac .* X .* bio;

end