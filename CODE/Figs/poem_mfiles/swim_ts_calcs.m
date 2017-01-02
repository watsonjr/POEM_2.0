wgt = 2.5; 	
T = 30.0;
w = ((3.9*wgt.^0.13 * exp(0.149*T)) /100);
dw = w * 60 * 60 *24

wgt2 = 2.5*1e3; 	
w2 = ((3.9*wgt2.^0.13 * exp(0.149*T)) /100);
dw2 = w2 * 60 * 60 *24

min(dyte(:))
mean(dyte(:))
median(dyte(:))

min(dxtn(:))
mean(dxtn(:))
median(dxtn(:))

min(sqrt(dat(:)))
mean(sqrt(dat(:)))
median(sqrt(dat(:)))

dist = min(sqrt(dat(:)));
dist / (w * 60 * 60)
dist / (w2 * 60 * 60)