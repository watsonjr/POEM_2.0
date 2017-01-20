%%% Update biomass
function [bio_out, pred] = sub_update_be(bio_in,BE,det,CC,con,bio)
    %bio_in = benthic biomass
    %con = biomass specific consumption rate by MD and LD
    %bio = biomass of MD and LD
    
    %! Peterson & Wroblewski: daily & uses dry weight
    %nmort = exp(0.063*(temp-15.0)) * 5.26e-3 * (wgt/9.0)^(-0.25);
    % Use Wei et al mean (1.3956 gWW) or median (0.7656) wgt
%     wgt = 0.7656;
%     omort = exp(0.063*(Tb-15.0)) * 5.26e-3 * (wgt/9.0)^(-0.25);
%     eaten = con.*bio;
%     pred = sum(eaten);
%     bio_out = bio_in - sum(eaten) - bio_in*omort;
    
    %! Quadratic mortality from carrying capacity
  r = BE*det;
  eaten = con.*bio;
  pred = sum(eaten);
  % Chemostat
  bio_out = bio_in + r * (r*CC - bio_in) - pred;
  % Logistic
  %bio_out = bio_in + r * bio_in * (1 - bio_in/CC) - pred;
end
