%this will compute some stuff from the daymet swe dataset
% dave hill
% jun 2022

clear all
close all
clc
addpath ktaub

pathname='/nfs/attic/dfh/data/daymet/orders/data/Daymet_Daily_V4/data';
yearrange=1980:2021;
numyears=length(yearrange);

nx=7814;ny=8075; %size of grid for daymet data...

%check to see if snowcovered days has already been calculated. If not, then
%calculate it and save file.
if ~exist('snowcovereddays.mat')
    for j=1:numyears
        yearrange(j)
        filename=fullfile(pathname,['daymet_v4_daily_na_swe_' num2str(yearrange(j)) '.nc']);
        %loop over a water year to compute no. of snow covered days
        for k=1:365
            k;
            %read in daily swe.
            swe=ncread(filename,'swe',[1 1 k],[Inf Inf 1]);
            if k==1
                counter=zeros(size(swe));
            end
            I=find(swe>0);
            swe(I)=1;swe(~I)=0;
            counter=counter+swe;
        end
            numdays(:,:,j)=counter; %num snow covered days for that year
    end
    %note this array is 7814 x 8075 x 42 (years)
    m=matfile('snowcovereddays.mat','Writable',true);
    m.numdays=numdays;
else
    disp('File already exists.')
    m=matfile('snowcovereddays.mat');
    numdays=m.numdays;
end
%next, move on to compute trend analysis....Only do this for cells that do
%not contain NaNs for any years. This will mask out ocean, lakes, and so
%on.

if ~exist('pvals.mat') & ~exist('slopes.mat')
    
    %initialize grids with nans
    slope=NaN*ones(nx,ny);
    pvalue=NaN*ones(nx,ny);
    
    %Find locations of valid cells.
    tmp=mean(numdays,3);
    [I,J]=find(isnan(tmp)==0);  %locations of valid cells
    %begin loop
    for j=1:length(I)
        tmp=numdays(I(j),J(j),:);
        inputdata=[yearrange' tmp(:)];
        [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] ...
            = ktaub(inputdata, .1, 0);
        slope(I(j),J(j))=sen;
        pvalue(I(j),J(j))=sig;
        if mod(j,1000)==0
            disp([num2str(j/length(I)) ' done'])
        end
    end
    n=matfile('pvals.mat','Writable',true);
    n.pvalue=pvalue;
    o=matfile('slopes.mat','Writable',true);
    o.slope=slope;
    else
    disp('Files already exist.')
    load pvals.mat
    load slopes.mat
end

%let's now make some plots...
%get lon / lat from one grid
filename=fullfile(pathname,['daymet_v4_daily_na_swe_' num2str(yearrange(1)) '.nc']);
Lon=ncread(filename,'lon');
Lat=ncread(filename,'lat');

%plot number of snow covered days...
meandays=nanmean(numdays,3);
clear numdays
figure(1)
p1=geoshow(Lat,Lon,meandays,'DisplayType','surface')
xlabel('Lon');ylabel('Lat');title('Mean annual snow covered days 1980-2021')
colorbar
caxis([0 320])
axis([-180 -50 10 85])
hold on
states = shaperead('usastatehi','UseGeoCoords',true);
p2=geoshow([states.Lat],[states.Lon],'Color','black','LineWidth',0.5);
p2.ZData(1:length(p2.YData))=400;

figure(2)
I=find(meandays==0);
maskedmeandays=meandays;
maskedmeandays(I)=NaN;
p1=geoshow(Lat,Lon,maskedmeandays,'DisplayType','surface')
xlabel('Lon');ylabel('Lat');title('Mean annual snow covered days 1980-2021 (zeros masked out)')
colorbar
caxis([0 320])
axis([-180 -50 10 85])
hold on
states = shaperead('usastatehi','UseGeoCoords',true);
p2=geoshow([states.Lat],[states.Lon],'Color','black','LineWidth',0.5);
p2.ZData(1:length(p2.YData))=400;

%plot the trend (days per year) in the number of snow covered days.

figure(3)
p1=geoshow(Lat,Lon,slope,'DisplayType','surface')
xlabel('Lon');ylabel('Lat');title('Rate of change (days / yr) in snow covered days')
colorbar
caxis([-1 1])
axis([-180 -50 10 85])
hold on
states = shaperead('usastatehi','UseGeoCoords',true);
p2=geoshow([states.Lat],[states.Lon],'Color','black','LineWidth',0.5);
p2.ZData(1:length(p2.YData))=400;

xxx

figure(4)
sigslope=slope;
I=find(pvalue>0.1);
sigslope(I)=NaN;
p1=geoshow(Lat,Lon,sigslope,'DisplayType','surface')
xlabel('Lon');ylabel('Lat');title('Rate of change (days / yr) in snow covered days (p < 0.05)')
colorbar
caxis([-1 1])
axis([-125 -65 25 50])
hold on
states = shaperead('usastatehi','UseGeoCoords',true);
p2=geoshow([states.Lat],[states.Lon],'Color','black','LineWidth',0.5);
p2.ZData(1:length(p2.YData))=400;


%look at climate velocity results.
%compute gradients
%[dx,dy]=gradient(meandays,lat,lon);
[aspect,slope2,gradN,gradE]=gradientm(Lat,Lon,meandays);

%figure(5);
%p1=geoshow(Lat,Lon,slope2,'DisplayType','surface')

slopesmooth=imgaussfilt(slope,'Filtersize',19);
figure
pcolor(slopesmooth');shading flat.

%compute vector components...
vx=slope./gradE;
vy=slope./gradN;
speed=sqrt(vx.^2+vy.^2);
%figure(5)
%p1=geoshow(Lat,Lon,speed,'DisplayType','surface');


