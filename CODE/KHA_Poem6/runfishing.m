baserun

%Fishing rates
F = linspace(0, 1, 10);
%Length of run (years)
param.tEnd = 150;
for i = 1:length(F)
    param.F = [0 0.1*F(i) F(i)]';
    %param.F = [0 F(i)  0 0 F(i)  0 0 F(i)]';
    result = poem(param, result);
    Y(i) = sum(result.Y(end,:));
    B(i,:) = result.B(end,:);
end
%Plot yield vs. fishing rate
clf
plot(F,Y,'o-')

