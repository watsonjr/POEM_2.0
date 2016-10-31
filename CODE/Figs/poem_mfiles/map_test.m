figure(1)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(AllD2)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%%
% figure(2)
% axesm ('Stereographic','lat',45,'lon',-157.5,'rad',90)
% surfm(geolat_t,geolon_t,real(log10(AllD2)))
% colormap('jet')              
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%%
% figure(3)
% axesm ('Orthographic','lat',45,'lon',-157.5,'rad',90)
% surfm(geolat_t,geolon_t,real(log10(AllD2)))
% colormap('jet')              
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%%
% figure(4)
% axesm ('Azimuthal Equal-area','lat',45,'lon',-157.5,'rad',90,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(geolat_t,geolon_t,real(log10(AllD2)))
% colormap('jet')              
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%%
% figure(5)
% axesm ('Azimuthal Equidistant','lat',45,'lon',-157.5,'rad',90,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(geolat_t,geolon_t,real(log10(AllD2)))
% colormap('jet')              
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%% just weird looking
% figure(6)
% axesm ('Gnomonic','MapLatLimit',latlim,'MapLonLimit',[-255 -60],'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(geolat_t,geolon_t,real(log10(AllD2)))
% colormap('jet')              
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%%
% figure(7)
% axesm ('Albers Equal-Area Conic','lat',latlim,'lon',[-255 -60],'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(geolat_t,geolon_t,real(log10(AllD2)))
% colormap('jet')              
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%%
% figure(8)
% axesm ('Lambert Conformal Conic','lon',[-255 -60],'lat',latlim)
% surfm(geolat_t,geolon_t,real(log10(AllD2)))
% colormap('jet')              
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%%
figure(9)
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(AllD2)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%%
% figure(10)
% axesm ('Miller Cylindrical','MapLatLimit',latlim,'MapLonLimit',[-255 -60],'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(geolat_t,geolon_t,real(log10(AllD2)))
% colormap('jet')              
% load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%%
% figure(11)
% axesm ('Equidistant Cylindrical','MapLatLimit',latlim,'MapLonLimit',[-255 -60],'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(geolat_t,geolon_t,real(log10(AllD2)))
% colormap('jet')              
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%%
% figure(12)
% axesm ('Oblique Mercator','MapLatLimit',latlim,'MapLonLimit',[-255 -60],'aspect',.8)
% surfm(geolat_t,geolon_t,real(log10(AllD2)))
% colormap('jet')              
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%%
% figure(13)
% axesm ('Transverse Mercator','MapLatLimit',latlim,'MapLonLimit',[-255 -60],'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(geolat_t,geolon_t,real(log10(AllD2)))
% colormap('jet')              
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%%
% figure(14)
% for l=6:7,
%  m_proj('sinusoidal','long',Slongs(l,:),'lat',Slats(l,:));
%  m_grid('fontsize',6,'xticklabels',[],'xtick',[-180:30:360],...
%         'ytick',[-80:20:80],'yticklabels',[],'linest','-','color',[.9 .9 .9]);
%  m_coast('patch','g');
% end;
% surfm(geolat_t,geolon_t,real(log10(AllD2)))
% colormap('jet')              
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%%
% figure(15)
% axesm ('Gall-Peters','MapLatLimit',latlim,'MapLonLimit',[-255 -60],'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(geolat_t,geolon_t,real(log10(AllD2)))
% colormap('jet')              
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%%
% figure(16)
% m_proj('Hammer-Aitoff','MapLatLimit',latlim,'MapLonLimit',[-255 -60],'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(geolat_t,geolon_t,real(log10(AllD2)))
% colormap('jet')              
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%%
figure(17)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(AllD2)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%% just weird looking
% figure(18)
% axesm ('UTM','MapLatLimit',latlim,'MapLonLimit',[-255 -60],'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(geolat_t,geolon_t,real(log10(AllD2)))
% colormap('jet')              
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);