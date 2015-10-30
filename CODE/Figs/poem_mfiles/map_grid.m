% Grid data for mapping

clear all
close all

dpath = '/Users/Colleen/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';

grid = csvread([dpath 'grid_csv.csv']);
%fix lon shift
id=find(grid(:,2)<-180);
grid(id,2)=grid(id,2)+360;

%%
x=-180:180;
y=-90:90;
[X,Y]=meshgrid(x,y);
Z=griddata(grid(:,2),grid(:,3),grid(:,4),X,Y);

% Interpolate model output to map
%Check that it worked with depth

figure(1)
m_proj('miller','lat',82);
m_pcolor(X,Y,Z); hold on;
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;

%% Make higher res map
x=-180:0.25:180;
y=-90:0.25:90;
[X,Y]=meshgrid(x,y);
Z=griddata(grid(:,2),grid(:,3),grid(:,4),X,Y);

figure(2)
m_proj('miller','lat',82);
m_pcolor(X,Y,Z); hold on;
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid('linestyle','none');