% sn_blueMOT_analysis.m
% Optimized analysis for Sn blueMOT dynamics with robust file naming from original pos strings

clearvars; close all; clc;

%% User settings
simTrapDynamics = true;       % simulate trajectories with photon scatter if true
moleculeName    = 'Sn';        % molecule identifier
bfieldGrad      = 15;          % radial B-field gradient (G/cm)
dataTimestamp   = '20250813_1035';
laserWaistMM    = 7;           % beam waist radius (mm)
useRandPhase    = true;        % include random laser phases

%% Data folder setup
baseFolder = 'SnData/SnCBMOTBGradGPerCM';
subName    = sprintf('%dWaistMM%dDate%s', bfieldGrad, laserWaistMM, dataTimestamp);
if useRandPhase
    subName = [subName '_WithRandPhase'];
end
% concatenate without extra slash to match existing on-disk structure
dataFolder = strcat(baseFolder, subName);

%% Verify data folder exists
if ~isfolder(dataFolder)
    error('Data folder does not exist: %s', dataFolder);
end

%% Position settings: use original string list for exact filenames
posStrings = {'0.01','0.03','0.05','0.1','0.3','0.6','0.9','1.2','1.5','2.0','2.5','3.0','4.0'};
posListMM  = str2double(posStrings);  % numeric for interpolation
nPos       = numel(posStrings);

%% Load simulation data with existence checks
for i = 1:nPos
    % use posStrings to match exact .dat names
    fname = strcat(dataFolder, '/forcevsSpeedDisplacement', posStrings{i}, 'MMSameDir.dat');
    if ~isfile(fname)
        % debug: list available files
        tmp = dir(fullfile(dataFolder, '*.dat'));
        fprintf('Available .dat files in %s:\n', dataFolder);
        fprintf(' %s\n', tmp.name);
        error('File not found: %s', fname);
    end
    T = readtable(fname);
    if i == 1
        vels   = T.Speed;               % velocity vector (m/s)
        nV     = numel(vels);
        accV   = zeros(nV, nPos);       % acceleration along v-dir (µm/ms^2)
        popExc = zeros(nV, nPos);       % excited-state fraction
    end
    accV(:,i)   = T.av;               % acceleration in v-direction
    popExc(:,i) = T.PFeHigh;          % excited-state population
end

%% Symmetrize and interpolate setup
negIdx = vels < 0;
accV(negIdx, :) = -accV(negIdx, :);
accVfull = [-flipud(fliplr(accV)), accV];
popFull  = [flipud(fliplr(popExc)), popExc];
sortedPos = sort([-posListMM, posListMM]);
oneAxisAccel = @(d, v) interp2(sortedPos, vels, accVfull/1e3, d, v, 'spline');

%% Heat map ROI settings
d_max = 1.0; v_max = 2.0; step = 0.02;
xs = -d_max:step:d_max;
vs = -v_max:step:v_max;
[Xs, Vs] = meshgrid(xs, vs);
heatMap = arrayfun(@(d,v) oneAxisAccel(d,v), Xs, Vs);

%% Plot 2D heat map
figure;
imagesc(xs, vs, heatMap);
axis xy;          % ensure correct orientation
axis square;      % force square shape
axis tight;       % fit axes limits to data
caxis([min(heatMap(:)), max(heatMap(:))]);
cb = colorbar;
cb.Title.String = 'a_{d} (mm/ms^{2})';
xlabel('d (mm)');
ylabel('v_{d} (m/s)');
title('Heat Map (|d|≤0.5 mm, |v|≤1 m/s)');

%% Compute average acceleration curves
dMaxInt = 1.0; vMaxInt = 2.0;
dIdx = abs(xs)<=dMaxInt;
vIdx = abs(vs)<=vMaxInt;
accVsPos = trapz(vs(vIdx), heatMap(vIdx,:),1)/(2*vMaxInt);
accVsVel = trapz(xs(dIdx), heatMap(:,dIdx),2)/(2*dMaxInt);

%% Plot avg accel vs position
figure;
plot(xs, accVsPos,'LineWidth',1.5);
grid on;
xlabel('d (mm)');
ylabel('a_{d} (mm/ms^{2})');
title('Average Accel vs Position');

%% Plot avg accel vs velocity
figure;
plot(vs, accVsVel,'LineWidth',1.5);
grid on;
xlabel('v_{d} (m/s)');
ylabel('a_{d} (mm/ms^{2})');
title('Average Accel vs Velocity');

%% Optional trap-dynamics simulation
if simTrapDynamics
    meanExcPop = mean(popFull(abs(vels)<=vMaxInt, abs(sortedPos)<=dMaxInt),'all');
    mol = def_molecule(moleculeName);
    scatterRate = meanExcPop * mol.gam * 1e-3;
    tKick = 1/scatterRate; vKick = mol.hbar*mol.kA/mol.mass;
    maxTime = 20;
    r = 0; v = 0.3;
    for k=1:round(maxTime/tKick)
        phi1 = 2*pi*rand; phi2 = 2*pi*rand;
        v(end+1) = v(end) + vKick*(cos(phi1)+cos(phi2)) + oneAxisAccel(r(end),v(end))*tKick;
        r(end+1) = r(end) + v(end)*tKick;
    end
    figure;
    plot(r, v,'LineWidth',2);
    grid on;
    xlabel('d (mm)');
    ylabel('v_{d} (mm/ms)');
    title('Trajectory with Photon Scatter');
    legend('Location','best');
end