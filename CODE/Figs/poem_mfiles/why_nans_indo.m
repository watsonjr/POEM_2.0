%%
all_mg = NaN*ones(size(geolon_t));
all_mg(grid(:,1)) = all_bio(:,15);

log10_all_mg = log10(all_mg);
log10_all_bio = log10(all_bio(:,15));

figure(50)
axesm ('mollweid','MapLatLimit',[-90 90],'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(all_mg)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
hcb = colorbar('h');
ylim(hcb,[-1 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('1951-2000 log10 mean biomass All Fishes (g m^-^2)')
stamp(cfile)
%print('-dpng',[ppath 'Hist_pristine_All_v2.png'])

%%
on = find(grid(:,2)>90 & grid(:,2)<120);
at = find(grid(:,3)>0 & grid(:,3)<60);
ll = intersect(on,at);

test1 = all_bio(ll,15);
test2 = log10_all_bio(ll);
test3 = b_mean5000(ll);

%% 1800-1850
time = 1:size(SP.bio,2);
mo=(time-1)/12;
mo=mo+1760;
yr50=find(mo>=1800 & mo<1850);

sp_smean=mean(SP.bio(:,yr50),2);
neg=find(sp_smean<=0);
nn=find(isnan(sp_smean));
test=grid(nn,:);
id=find(test(:,2)<-180);
test(id,2)=test(id,2)+360;

ld_smean=mean(LD.bio(:,yr50),2);
neg2=find(ld_smean<=0);
nn2=find(isnan(ld_smean));
test2=grid(nn2,:);
id=find(test2(:,2)<-180);
test2(id,2)=test2(id,2)+360;





