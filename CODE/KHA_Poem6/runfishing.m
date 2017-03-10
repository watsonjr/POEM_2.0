baserun

F = linspace(0, 1, 10);

param.tEnd = 50;
for i = 1:length(F)
    param.F = [0.1*F(i) F(i)  0 0 0     0 0 0   ]';
    %param.F = [0 F(i)  0 0 F(i)  0 0 F(i)]';
    result = poem(param, result);
    Y(i) = sum(result.Y(end,:));
    B(i,:) = result.B(end,:);    
end

clf
plot(F,Y,'o-')

