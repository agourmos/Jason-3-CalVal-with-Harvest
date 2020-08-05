%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adam Gourmos
% Harvest Cal/Val for Professor Leben
% Corrected Data
clear all;
clc;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SECTION A: Define all variables
latHarv = 34.4688755;
lonHarv = -120.68267131;
invNum = 0;
count = 0;
countVar = 1;
countBub = 1;
countRad = 1;
countLat = 1;
countCol = 1;
numBins = 21;
cdfSla1pt = [];
cdfSla4pt = [];
cdfSla9pt = [];
cdfSla16pt = [];
cdfTimeHold = [];
cdfTimeFinal = [];
cdfSlaFinal = [];
cdfSlaExtrapFinal = [];
noaaTime = [];
noaaN1_1 = [];
idxX = [];
binVec = [33.57, 33.61, 33.66, 33.71, 33.76, 33.81, 33.855, 33.905,... 
    33.95, 34, 34.045, 34.09, 34.14, 34.19, 34.235, 34.28, 34.33, 34.38, 34.43,...
    34.47, 34.52 34.6]; % Built after plotting all lon/lat values
% Building to seperate major values by bin. Seperated by latitude values
latAvg = zeros(1,21);
slaAvg = zeros(1,21);
binSla = NaN*ones(21,131); 
binLat = NaN*ones(21,131); 
binLon = NaN*ones(21,131); 
binTime = NaN*ones(21,131); 
binDry = NaN*ones(21,131);
cdfSlaExtrap = NaN*ones(1,131);

%% SECTION B: Pull the data from the .mat file and edit them
matfile = 'lidaroverflights.mat';
data = load(matfile);
time = data.lidardataoverflights.time;
noaaTimeComp = time;
time.Format = 'yyyy-MM-dd'; %Converts time to a comparable format for the netCDF data
noaaTimeComp.Format =  'yyyy-MM-dd HH:mm';
lid = data.lidardataoverflights.l_mean*100;
bub = data.lidardataoverflights.bub*100;
rad = data.lidardataoverflights.rad*100;

% This section finds the index of where the NaN's are in rad and bubbler
% data in order to add them back in for graphing and rkae them out for
% further analysis
bubNaN = find(isnan(bub));
bubNotNaN = find(~isnan(bub));
radNaN = find(isnan(rad));
radNotNaN = find(~isnan(rad));

% This for statement gets rid of all 'NaN's from the .mat/.csv file
bub = bub(~isnan(bub));
rad = rad(~isnan(rad));

% Pulling mean values for each list and applying them to correct for any
% mean values amongst data sets
lidMean = mean(lid);
bubMean = mean(bub);
radMean = mean(rad);

lid = (lid-lidMean)*-1;
bub = bub-bubMean;
rad = rad-radMean;

lidMean = mean(lid);
bubMean = mean(bub);
radMean = mean(rad);

% This next section finds values that are outside of n standard deviations
% and removes for bub, rad, and lid data
n = 3;

% Lidar Outliars
not_ind_lid = find((lid < (n * std(lid) + lidMean)) & ...
    (lid >( -n *std(lid) +  lidMean)));
lid = lid(not_ind_lid);

% Bubbler Outliars
not_ind_bub = find(~((bub < (n * std(bub) + bubMean)) & ...
    (bub >( -n *std(bub) +  bubMean))));
bub(not_ind_bub) = NaN;
bubHold = bub;

% Radar Outliars
not_ind_rad = find(~((rad < (n * std(rad) + radMean)) & ...
    (rad >( -n *std(rad) +  radMean))));
rad(not_ind_rad) = NaN;
radHold = rad;

% Redefines bub and rad arrays
bub = [];
rad = [];

% This next part insert NaN into the rad and bub vectors
% For bubbler data
for i = 1:(length(bubNaN)+length(bubNotNaN))
    if any(i == bubNotNaN)
        bub = [bub bubHold(countBub)];
        countBub = countBub +1;
    else
        bub = [bub NaN];
    end
end
% For radar data
for i = 1:(length(radNaN)+length(radNotNaN))
    if any(i == radNotNaN)
        rad = [rad radHold(countRad)];
        countRad = countRad +1;
    else
        rad = [rad NaN];
    end
end

%% SECTION C: This section pulls all data from the netCDF files pulled from RADs
cd radsData % Opens the directory to where all the netCDF files are
cdfData = dir; %Names all files in the directory
cdfData(1:2) = []; %Removes . and .. from list in directory
countCDF = 1; %Used for indexing

for i = 1:length(cdfData)
    name = cdfData(i).name;
    if length(name) == 17
        cycles(i) = str2double(name(12:15));
    else
        cycles(i) = str2double(name(12:14));
    end
    vecLen = length(ncread(name,'time'));
    cdfTime(countCDF:(countCDF+vecLen-1)) = ncread(name,'time');
    cdfSla(countCDF:(countCDF+vecLen-1)) = ncread(name,'sla');
    cdfLat(countCDF:(countCDF+vecLen-1)) = ncread(name,'lat');
    cdfLon(countCDF:(countCDF+vecLen-1)) = ncread(name,'lon');
    cdfDac(countCDF:(countCDF+vecLen-1)) = ncread(name,'dac');
    cdfDry(countCDF:(countCDF+vecLen-1)) = ncread(name,'dry_tropo_ecmwf');
    cdfIono(countCDF:(countCDF+vecLen-1)) = ncread(name,'iono_alt');
    cdfMss(countCDF:(countCDF+vecLen-1)) = ncread(name,'mss_dtu15');
    cdfRef(countCDF:(countCDF+vecLen-1)) = ncread(name,'ref_frame_offset');
    cdfSsb(countCDF:(countCDF+vecLen-1)) = ncread(name,'ssb_cls');
    cdfLoad(countCDF:(countCDF+vecLen-1)) = ncread(name,'tide_load_fes14');
    cdfOcean(countCDF:(countCDF+vecLen-1)) = ncread(name,'tide_ocean_fes14');
    cdfPole(countCDF:(countCDF+vecLen-1)) = ncread(name,'tide_pole');
    cdfSolid(countCDF:(countCDF+vecLen-1)) = ncread(name,'tide_solid');
    cdfWet(countCDF:(countCDF+vecLen-1)) = ncread(name,'wet_tropo_rad');
    countCDF = countCDF + vecLen;
end   
cd ../ %Returns to the original directory
% The time correction here finds the difference between POSIX and UNIX and
% adds it to cdfTime in order to align with datetime natural fucntionality
% of working with POSIX
timeDiff = 86400 * (datenum('01-Jan-1985') - datenum('01-Jan-1970'));
cdfTime = datetime(cdfTime+timeDiff,'ConvertFrom','posixtime'); 
cdfTime.Format = 'yyyy-MM-dd'; %Converts time to a comparable format with

%% SUB-SECTION C: This section pulls all data from the Diagtool files and 
% preps them for comaprison between multiple datasets
cd NoaaData % Opens the directory to where all the netCDF files are
noaaData = dir; %Names all files in the directory
noaaData(1:2) = []; %Removes . and .. from list in directory
noaaCount = 1;

% This for statement pulls data from all the files an fixes them
for i = 1:length(noaaData)
    vars = readtable(noaaData(i).name);
    noaaTime = [noaaTime datetime(table2cell(vars(:,1)),'InputFormat','yyyy/MM/dd HH:mm')'];
    noaaN1_1 = [noaaN1_1 table2array(vars(:,6))'];
end
cd ../ %Returns to the original directory
noaaTime.Format = 'yyyy-MM-dd HH:mm';
noaaTime = noaaTime';
noaaN1_1 = noaaN1_1'*100;
noaaMean = nanmean(noaaN1_1);
noaaN1_1 = noaaN1_1 - noaaMean;
noaaN1Final = zeros(length(noaaTimeComp),1);

% This next sectiosn compares dates AND time to the .csv time points
for i = 1:length(noaaTime)
   if seconds(noaaTime(i)-noaaTimeComp(noaaCount)) > seconds(0)
       valHold = noaaN1_1(i);

       while isnan(valHold)
           valHold = noaaN1_1(i+valCount);
           valCount = valCount+1;
       end
       noaaN1Final(noaaCount) = valHold;
       noaaCount = noaaCount + 1;
       valCount = 1;
   else
   end
   if noaaCount == length(noaaTimeComp)+1
       break
   end
end
noaaTime.Format =  'yyyy-MM-dd';

%% SECTION D: This section corrects the netCDF data by applying corrections
% through reverse-RPN, smoothing data, date comparison, and pulling closest 
% data to harvest

% Correction through reverse-RPN
cdfSla = (cdfSla+cdfDac+cdfSolid+cdfLoad+cdfRef+cdfOcean)*100; 

% Begin smoothing data by averaging as well as seperating each cycle by
% through datetime analysis. Also finds the latitude and ssh and puts them
% into a matrix and seperates them by bin
for i = 1:length(cdfTime)
    if i == (length(cdfTime)) % When the last cycle is being analyzed
        distances = minDist(cdfLat((i-count):i),cdfLon((i-count):i));
        [val, idx] = sort(distances);
        slaHeld = cdfSla((i-count):i);
        latHeld = cdfLat((i-count):i);
        cdfSla1pt = [cdfSla1pt slaHeld(idx(1))];
        cdfSla4pt = [cdfSla4pt mean(slaHeld(idx(1:4)))];
        cdfSla9pt = [cdfSla9pt mean(slaHeld(idx(1:9)))];
        cdfSla16pt = [cdfSla16pt mean(slaHeld(idx(1:16)))];
        cdfTimeHold = [cdfTimeHold cdfTime(i)];
        latHeld = cdfLat((i-count):i);
        lonHeld = cdfLon((i-count):i);
        dryHeld = cdfDry((i-count):(i));
        for j = 1:(length(binVec)-1)
            if countLat <= length(slaHeld)
                if binVec(j) < latHeld(countLat) && binVec(j+1) > latHeld(countLat)
% Compares the latitude at each binned value to the actual latitudes and 
% seperates them by bin    
                    binSla(j,countCol) = slaHeld(countLat);
                    binLat(j,countCol) = latHeld(countLat);
                    binLon(j,countCol) = lonHeld(countLat);
                    binDry(j,countCol) = dryHeld(countLat);
                    countLat = countLat +1;
                else
                    continue
                end
            else
                continue
            end
        end
    elseif datestr(cdfTime(i),'yyyy-mm-dd') == datestr(cdfTime(i+1),...
            'yyyy-mm-dd')
    % This condition occurs when the date is the same as the next value
    % (same cycle)
        count = count + 1;
    else % This condition occurs when the cycle changes.
        distances = minDist(cdfLat((i-count):(i+1)),cdfLon((i-count):(i+1)));
        [val, idx] = sort(distances);
        slaHeld = cdfSla((i-count):(i));
        cdfSla1pt = [cdfSla1pt mean(slaHeld(idx(1)))];
        cdfSla4pt = [cdfSla4pt mean(slaHeld(idx(1:4)))];
        % These try statements check to see if there are less then 9 or 16
        % values in a single cycle. If there is, it inserts a NaN.
        try 
            cdfSla9pt = [cdfSla9pt mean(slaHeld(idx(1:9)))];
        catch 
            cdfSla9pt = [cdfSla9pt NaN];
        end
        try
            cdfSla16pt = [cdfSla16pt mean(slaHeld(idx(1:16)))];
        catch 
            cdfSla16pt = [cdfSla16pt NaN];
        end
        cdfTimeHold = [cdfTimeHold cdfTime(i)];
        latHeld = cdfLat((i-count):(i));
        lonHeld = cdfLon((i-count):(i));
        dryHeld = cdfDry((i-count):(i));
        for j = 1:(length(binVec)-1)
            if countLat <= length(slaHeld) 
                if binVec(j) < latHeld(countLat) && binVec(j+1) > latHeld(countLat)
% Compares the latitude at each binned value to the actual latitudes and 
% seperates them by bin                    
                    binSla(j,countCol) = slaHeld(countLat);
                    binLat(j,countCol) = latHeld(countLat);
                    binLon(j,countCol) = lonHeld(countLat);
                    binDry(j,countCol) = dryHeld(countLat);
                    countLat = countLat +1;
                else
                    continue
                end
            else
                continue
            end
        end
        countLat = 1;
        count = 0;
        countCol = countCol + 1;
    end
end
% Sorting each cycle with each point
[cycles, idx] = sort(cycles);
cdfSla1pt = cdfSla1pt(idx);
cdfSla4pt = cdfSla4pt(idx);
cdfSla9pt = cdfSla9pt(idx);
cdfSla16pt = cdfSla16pt(idx);
cdfTimeHold = cdfTimeHold(idx);
binSla = binSla(:,idx);
binLat = binLat(:,idx);
binLon = binLon(:,idx);

% This for statement averages the latitude at each bin 
[r, c] = size(binLat);
for i = 1:r
    holdVec = binLat(i,:);
    holdVec2 = binSla(i,:);
    holdIdx = (~isnan(holdVec));
    latAvg(i) = mean(holdVec(holdIdx));
    slaAvg(i) = mean(holdVec2(holdIdx));
end

% This for statement finds the last 5 values in each cycle and uses them to
% create a linear fit, then find the theroretical value at harvest. If
% there are less than 2 points, then the value is set to NaN.
nn = 3;
binSlaComp = binSla;
binLatComp = binLat;
for i = 1:length(cdfSlaExtrap)
    holdLat = binLat(:,i);
    extrapLat = holdLat(~isnan(holdLat));
    holdSla = binSla(:,i);    
    extrapSla = holdSla(~isnan(holdSla));
    extrapMean = mean(extrapSla);
    not_ind_sla = find((extrapSla < (nn * std(extrapSla) + extrapMean)) & ...
     (extrapSla >( -nn *std(extrapSla) +  extrapMean)));
    holdSla2 = extrapSla(not_ind_sla);
    holdLat2 = extrapLat(not_ind_sla);
    if length(holdSla) < 2
        continue
    else
        p = polyfit(holdLat2,holdSla2,1);
        cdfSlaExtrap(i) = (p(1)*latHarv)+p(2);
    end
    nanIdx = find(isnan(holdLat));
    for j = 1:length(nanIdx)
        latHoldOn = latAvg(nanIdx(j));
        binLatComp(nanIdx(j),i) = latHoldOn;
        binSlaComp(nanIdx(j),i) = (p(1)*latHoldOn)+p(2); 
    end
end

% This next for statement compares the dates from the netCDF file to the
% .csv data in order to compare the 2 datasets
for i = 1:length(cdfTimeHold)
    for j = 1:length(time)
        if datestr(cdfTimeHold(i),'yyyy-mm-dd') == datestr(time(j),'yyyy-mm-dd')
            cdfTimeFinal = [cdfTimeFinal cdfTimeHold(i)];
            cdfSlaFinal = [cdfSlaFinal cdfSla1pt(i)];
            cdfSlaExtrapFinal = [cdfSlaExtrapFinal cdfSlaExtrap(i)];
        end
    end
end

% Find the root mean square of the RADs data
rmsSub = cdfSlaFinal - lid';
rms = sqrt(mean(rmsSub.^2));

% Find the root mean square of the NOAA data
rmsNoaaSub = cdfSlaFinal - noaaN1Final';
rmsNoaa = sqrt(mean(rmsNoaaSub.^2));

% Repeat the rms analysis for 4pt mean vector
rmsSub4pt = cdfSla1pt - cdfSla4pt;
rms4pt = sqrt(mean(rmsSub4pt.^2));

% Repeat the rms analysis for 9pt mean vector
rmsSub9pt = cdfSla1pt - cdfSla9pt;
rms9pt = sqrt(nanmean(rmsSub9pt.^2));

% Repeat the rms analysis for 16pt mean vector
rmsSub16pt = cdfSla1pt - cdfSla16pt;
rms16pt = sqrt(nanmean(rmsSub16pt.^2));

% Find the fit for the RMS data
xFill = datenum(cdfTimeFinal);
xxFill = linspace(min(xFill),max(xFill),100);
dd = datetime(datestr(xxFill));
dd.Format = 'yyyy-MM-dd';
fitted_y = polyval(polyfit(xFill,rmsSub,1),xxFill);
fitted_y_noaa = polyval(polyfit(xFill,rmsNoaaSub,1),xxFill);
cdfTimeFit = cdfTimeFinal;

% Finds the dates that are missing from the track and inserts NaN 
cdfTimeVar = cdfTimeHold(find(cdfTimeHold == cdfTimeFinal(1)):...
    find(cdfTimeHold == cdfTimeFinal(end))); % Finds orginal dates at 
% beggining and end of track
cdfTimeHold2 = NaT(1,length(cdfTimeVar));
cdfSlaHold2 = zeros(1,length(cdfTimeVar));
timeHold2 = NaT(1,length(cdfTimeVar));
lid2 = zeros(1,length(cdfTimeVar));
noaaTimeHold = NaT(1,length(cdfTimeVar));
noaaN1FinalHold = zeros(1,length(cdfTimeVar));
cdfSlaExtrapFinalHold = zeros(1,length(cdfTimeVar));
for i = 1:length(cdfTimeVar)
    if isempty(find(cdfTimeVar(i) == cdfTimeFinal)) 
        cdfTimeHold2(i) = cdfTimeVar(i);
        timeHold2(i) = cdfTimeVar(i);
        noaaTimeHold(i) = cdfTimeVar(i);
        cdfSlaHold2(i) = NaN;
        lid2(i) = NaN;
        noaaN1FinalHold(i) = NaN;
        cdfSlaExtrapFinalHold(i) = NaN;
        idxX = [idxX i];
    else
        cdfTimeHold2(i) = cdfTimeVar(i);
        cdfSlaHold2(i) = cdfSlaFinal(countVar);
        timeHold2(i) = cdfTimeVar(i);
        lid2(i) = lid(countVar);
        noaaTimeHold(i) = cdfTimeVar(i);
        noaaN1FinalHold(i) = noaaN1Final(countVar);
        cdfSlaExtrapFinalHold(i) = cdfSlaExtrapFinal(countVar);
        countVar = countVar +1;
    end
end
Xs = zeros(1,length(idxX));

% Redefines final variables
cdfTimeFinal = cdfTimeHold2;
cdfSlaFinal = cdfSlaHold2;
time = timeHold2;
lid = lid2;
noaaTimeFinal = noaaTimeHold;
noaaN1Final = noaaN1FinalHold;
cdfSlaExtrapFinal = cdfSlaExtrapFinalHold;

%% SECTION E: Begin plotting all types of data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1 plots flightpath as well as all data together
fig1 = figure(1);
% Subplot 1, figure(1)
ax1 = subplot(2,2,1);
hold on
plot(cdfTimeFinal,cdfSlaFinal,'r-*');
plot(time, lid,'b-*');
plot(time(idxX),Xs,'kx')
%scatter(time, bub,300,'.')
%scatter(time, rad,300,'.')
hold off
ylabel('SSH (cm)')
xtickangle(45)
grid on
ylim([-1.5,1.5]*100)
title('(a) Jason-3 Data Versus Lidar Data')
legend(ax1,{'Jason-3/RADs','lidar'},'Location','southwest')

% Subplot 2, figure(1)
ax2 = subplot(2,2,3);
hold on
scatter(cdfTimeFit, rmsSub,300,'.')
plot(dd,fitted_y,'LineWidth',2)
legend(ax2,{'Subtracted Data','Trendline'},'Location','southwest')
ylabel('Jason-3 Data - Lidar Data (cm)')
title('(b) Difference betewen Jason-3 and Lidar Data')
xtickangle(45)
grid on
hold off

% Subplot 3, figure(1)
ax3 = subplot(2,2,[2 4]);
hold on
grid on
worldmap([33.5 34.7],[-121.5 -120]);
%worldmap([30 40],[-125 -115]);
load coastlines
graph1 = geoshow(cdfLat,cdfLon,...
    'displaytype','point',...
    'Marker','.',...
    'MarkerSize',1,...
    'MarkerEdgeColor', 'b');
graph2 = geoshow(latHarv,lonHarv,...
    'displaytype','point',...
    'Marker','.',...
    'MarkerSize',20,...
    'MarkerEdgeColor', 'r');
idk = 34.6;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%geoshow([idk,idk],[-180,0])
title('(c) Jason-3 Overflight Path')
[latcells, loncells] = polysplit(coastlat, coastlon);
numel(latcells);
plotm(coastlat, coastlon,'k');
hold off
legend([graph1, graph2],{'Flightpath','Harvest'},'Location','southwest')

% Save figure
set(fig1, 'Position',  [0,0, 2000, 1200])% Change size of image
saveas(fig1,'HarvestCalVal.png')
close Figure 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2 plots the difference between different number averages
% Subplot 1, figure(2)
fig2 = figure(2);
ax1 = subplot(3,4,1);
hold on
graph1 = plot(cycles,cdfSla1pt,'k','LineWidth',1.2);
graph2 = plot(cycles,cdfSla4pt,'color','r','LineWidth',1.2);
scatter(cycles,cdfSla1pt,'*k')
scatter(cycles,cdfSla4pt,'*r')
ylabel('SSH (cm)')
title('Smotthing Comparison for n = 1,4,9,16')
xlim([30 161])
grid on
legend([graph1, graph2],{'1pt Average','4pt Average'},'Location','southwest')
hold off

% Subplot 2, figure(2)
ax2 = subplot(3,4,5);
hold on
graph3 = plot(cycles,cdfSla1pt,'k','LineWidth',1.2);
graph4 = plot(cycles,cdfSla9pt,'color','g','LineWidth',1.2);
scatter(cycles,cdfSla1pt,'*k')
scatter(cycles,cdfSla9pt,'*g')
ylabel('SSH (cm)')
xlim([30 161])
grid on
legend([graph3,graph4],{'1pt Average','9pt Average'},'Location','southwest')
hold off

% Subplot 3, figure(2)
ax3 = subplot(3,4,9);
hold on
graph5 = plot(cycles,cdfSla1pt,'k','LineWidth',1.2);
graph6 = plot(cycles,cdfSla16pt,'color','b','LineWidth',1.2);
scatter(cycles,cdfSla1pt,'*k')
scatter(cycles,cdfSla16pt,'*b')
ylabel('SSH (cm)')
xlabel('Cycle Number')
xlim([30 161])
grid on
legend([graph5,graph6],{'1pt Average','16pt Average'},'Location','southwest')
hold off

% Subplot 4, figure(2)
ax4 = subplot(3,4,[2 3 4 6 7 8 10 11 12]);
hold on
plot(cycles,cdfSla4pt-cdfSla1pt,'color','r','LineWidth',1.2)
plot(cycles,cdfSla9pt-cdfSla1pt,'color','g','LineWidth',1.2)
plot(cycles,cdfSla16pt-cdfSla1pt,'color','b','LineWidth',1.2)
legend({'1pt v 4pt Average','1pt v 9pt Average','1pt v 16pt Average'},...
    'Location','southwest')
title('Difference Between Neareest Point to Harvest to n-pt Averages')
xlim([30 161])
grid on
hold off

% Save figure 2
sgtitle('Harvest Cal/Val Analysis')
set(fig2, 'Position',  [0, 0, 2000, 1200])% Change size of image
saveas(fig2,'HarvestCalValAnalysis.png')
close Figure 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3
% 1pt v 4pt
fig3 = figure(3);
hold on
graph1 = plot(cycles,cdfSla1pt,'k','LineWidth',1.2);
graph2 = plot(cycles,cdfSla4pt,'color','r','LineWidth',1.2);
graph3 = plot(cycles,cdfSla4pt-cdfSla1pt,'--','color',...
    [0.6350, 0.0780, 0.1840],'LineWidth',1.2);
scatter(cycles,cdfSla1pt,'*k')
scatter(cycles,cdfSla4pt,'*r')
ylabel('SSH (cm)')
xlabel('Cycle Number')
title('Nearest Point to Harvest Versus the Nearest 4pt Average')
xlim([30 161])
grid on
legend([graph1, graph2, graph3],{'1pt Average','4pt Average','4pt-1pt'},...
    'Location','southwest')
hold off
set(fig3, 'Position',  [0, 0, 2000, 1200])% Change size of image
saveas(fig3,'CalVal1ptv4pt.png')
close Figure 3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4
% 1pt v 9pt
fig4 = figure(4);
hold on
graph1 = plot(cycles,cdfSla1pt,'k','LineWidth',1.2);
graph2 = plot(cycles,cdfSla9pt,'color','g','LineWidth',1.2);
graph3 = plot(cycles,cdfSla9pt-cdfSla1pt,'--','color',...
    [0 0.5 0],'LineWidth',1.2);
scatter(cycles,cdfSla1pt,'*k')
scatter(cycles,cdfSla9pt,'*g')
ylabel('SSH (cm)')
xlabel('Cycle Number')
title('Nearest Point to Harvest Versus the Nearest 9pt Average')
xlim([30 161])
grid on
legend([graph1, graph2, graph3],{'1pt Average','9pt Average','9pt-1pt'},...
    'Location','southwest')
hold off
set(fig4, 'Position',  [0, 0, 2000, 1200])% Change size of image
saveas(fig4,'CalVal1ptv9pt.png')
close Figure 4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5
% 1pt v 16pt
fig5 = figure(5);
hold on
graph1 = plot(cycles,cdfSla1pt,'k','LineWidth',1.2);
graph2 = plot(cycles,cdfSla16pt,'color','b','LineWidth',1.2);
graph3 = plot(cycles,cdfSla16pt-cdfSla1pt,'--','color',...
    [0 0 128/255],'LineWidth',1.2);
scatter(cycles,cdfSla1pt,'*k')
scatter(cycles,cdfSla16pt,'*b')
ylabel('SSH (cm)')
xlabel('Cycle Number')
title('Nearest Point to Harvest Versus the Nearest 16pt Average')
xlim([30 161])
grid on
legend([graph1, graph2, graph3],{'1pt Average','16pt Average','16pt-1pt'},...
    'Location','southwest')
hold off
set(fig5, 'Position',  [0, 0, 2000, 1200])% Change size of image
saveas(fig5,'CalVal1ptv16pt.png')
close Figure 5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6
% 1pt v 4pt, 9pt, 16pt
fig6 = figure(6);
hold on
plot(cycles,cdfSla4pt-cdfSla1pt,'color','r','LineWidth',1.2)
plot(cycles,cdfSla9pt-cdfSla1pt,'color','g','LineWidth',1.2)
plot(cycles,cdfSla16pt-cdfSla1pt,'color','b','LineWidth',1.2)
legend({'1pt v 4pt Average','1pt v 9pt Average','1pt v 16pt Average'},...
    'Location','southwest')
title('Difference Between Neareest Point to Harvest to n-pt Averages')
xlim([30 161])
ylabel('SSH (cm)')
xlabel('Cycle Number')
grid on
hold off
set(fig6, 'Position',  [0, 0, 2000, 1200])% Change size of image
saveas(fig6,'CalValAnalysisAll.png')
close Figure 6

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 7
% Plots rms value for each mean point value
fig7 = figure(7);
hold on
plot([4, 9, 16],[rms4pt, rms9pt, rms16pt],'r-*');
%plot(polyval(polyfit(sqrt([4, 9, 16]),[rms4pt, rms9pt, rms16pt],1),...
%    linspace(0,16,100)))
xlim([2, 18])
xlabel('Mean of n Closest Points to Harvest')
ylabel('Root Mean Square Value')
title('Measurement of Precision of Mean of n Closest Points to Harvest')
grid on
hold off
set(fig7, 'Position',  [0, 0, 1200, 1200])% Change size of image
saveas(fig7,'PrecisionOfAnalysis.png')
close Figure 7

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 8 plots flightpath as well as all data together
fig8 = figure(8);
% Subplot 1, figure(8)
ax1 = subplot(2,2,1);
hold on
plot(cdfTimeFinal,cdfSlaFinal,'g-*');
plot(time, noaaN1Final,'b-*');
plot(time(idxX),Xs,'kx')
%scatter(time, bub,300,'.')
%scatter(time, rad,300,'.')
hold off
ylabel('SSH (cm)')
xtickangle(45)
grid on
ylim([-1.5,1.5]*100)
title('(a) Jason-3 Data Versus NOAA Radar Data')
legend(ax1,{'Jason-3/RADs','NOAA Radar'},'Location','southwest')

% Subplot 2, figure(8)
ax2 = subplot(2,2,3);
hold on
scatter(cdfTimeFit, rmsNoaaSub,300,'.')
plot(dd,fitted_y_noaa,'LineWidth',2)
legend(ax2,{'Subtracted Data','Trendline'},'Location','southwest')
ylabel('Jason-3 Data - Lidar Data (cm)')
title('(b) Difference betewen Jason-3 and Lidar Data')
xtickangle(45)
grid on
hold off

% Subplot 3, figure(1)
ax3 = subplot(2,2,[2 4]);
hold on
grid on
worldmap([33.5 34.7],[-121.5 -120]);
%worldmap([30 40],[-125 -115]);
load coastlines
graph1 = geoshow(cdfLat,cdfLon,...
    'displaytype','point',...
    'Marker','.',...
    'MarkerSize',1,...
    'MarkerEdgeColor', 'b');
graph2 = geoshow(latHarv,lonHarv,...
    'displaytype','point',...
    'Marker','.',...
    'MarkerSize',20,...
    'MarkerEdgeColor', 'g');
idk = 33.57;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%geoshow([idk,idk],[-180,0])
title('(c) Jason-3 Overflight Path')
[latcells, loncells] = polysplit(coastlat, coastlon);
numel(latcells);
plotm(coastlat, coastlon,'k');
hold off
legend([graph1, graph2],{'Flightpath','Harvest'},'Location','southwest')

% Save figure
set(fig8, 'Position',  [0,0, 2000, 1200])% Change size of image
saveas(fig8,'HarvestCalVal_NOAA.png')
%close Figure 8

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 9
fig9 = figure(9);
hold on
plot(cdfTimeFinal,cdfSlaFinal,'g-*');
plot(time, noaaN1Final,'b-*');
plot(time, lid,'r-*');
plot(time(idxX),Xs,'kx')
hold off
ylabel('SSH (cm)')
xlabel('Date')
xtickangle(45)
grid on
ylim([-1.5,1.5]*100)
title({'(a) Jason-3 Data Versus NOAA Radar Data','Versus Lidar at Harvest Data'})
legend({'Jason-3/RADs','NOAA Radar', 'Lidar at Harvest'},'Location','southwest')

% Save figure
set(fig9, 'Position',  [0,0, 2000, 1200])% Change size of image
saveas(fig9,'HarvestCalValAllData.png')
close Figure 9

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 10 This section plots cdfSLA for the 1pt value vs the extrapolated value

fig10 = figure(10);
hold on
plot(cdfTimeFinal,cdfSlaFinal,'r-*');
%plot(time, lid,'c-*');
plot(cdfTimeFinal, cdfSlaExtrapFinal,'b-*');
plot(time(idxX),Xs,'kx')
%scatter(time, bub,300,'.')
%scatter(time, rad,300,'.')
hold off
ylabel('SSH (cm)')
xtickangle(45)
grid on
ylim([-1.5,1.5]*100)
title('(a) Jason-3 Data Versus Extrapolated Jason-3 Data')
legend({'Jason-3/RADs','Extrapolated Jason-3 Data'},'Location','southwest')

% Save figure
set(fig10, 'Position',  [0,0, 2000, 1200])% Change size of image
saveas(fig10,'HarvestCalVaExtrapProof.png')
close Figure 10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 11 pltos the map of the overlfight data only

fig11 = figure(11);
worldmap([33.5 34.7],[-121.5 -120]);
load coastlines
graph1 = geoshow(cdfLat,cdfLon,...
    'displaytype','point',...
    'Marker','.',...
    'MarkerSize',1,...
    'MarkerEdgeColor', 'b');
graph2 = geoshow(latHarv,lonHarv,...
    'displaytype','point',...
    'Marker','.',...
    'MarkerSize',20,...
    'MarkerEdgeColor', 'r');
for i = 1:length(binVec)
    hold on
    geoshow([binVec(i),binVec(i)],[-180,0])
end
grid on
title('(c) Jason-3 Overflight Path')
xlabel('Longitude')
ylabel('Latitude')
[latcells, loncells] = polysplit(coastlat, coastlon);
numel(latcells);
plotm(coastlat, coastlon,'k');

% Save figure
set(fig11, 'Position',  [0,0, 2000, 1200])% Change size of image
saveas(fig11,'Jason3OverflightRoute.png')
close Figure 11
