% Optimized Sn Capture Red (and compressed red) MOT Analysis Script with CSV exports (timestamped filenames)
% Exports: heat map, a(v), a(d), selected capture trajectories in one CSV,
% and one stochastic sample trajectory.

%% Setup and Parameters
clearvars; close all; clc;

% Simulation flags
simTrapDynamics = true;   % enable MOT temperature & size calculation
nTraj = 10;               % number of trajectories for temp/size

% Molecule and trap parameters
molecule = def_molecule("Sn");
gam    = molecule.gam;       % natural linewidth (1/ms)
kA     = molecule.kA;        % photon momentum (kg·mm/ms)
mass   = molecule.mass;      % atomic mass (kg)
hbar   = 1.0545718e-34;      % J·s
bkB    = 1.380649e-23;       % J/K

% Trap parameters
laser_waist   = 7;                               % mm
cap_threshold = sqrt(2)*laser_waist;            % capture volume radius (mm)

% ---- Timestamp from your dataset ----
dataStamp   = '20250831_2359';
gradField   = '25';
baseFolder  = sprintf('SnData/SnRedMOTBGradGPerCM%sWaistMM%.0fDate%s%s', ...
                      gradField, laser_waist, dataStamp, '_WithRandPhase');

% % ---- Export folder at the same level as this script ----
% thisFile  = mfilename('fullpath');
% if isempty(thisFile)
%     scriptDir = pwd;        % running as a section
% else
%     scriptDir = fileparts(thisFile);
% end
% exportFolder = fullfile(scriptDir, 'AllCSVAndMatDataForSn');
% if ~exist(exportFolder, 'dir'), mkdir(exportFolder); end

%% Load Simulation Data (capture red MOT case)
pos_mms  = (1:14)';                      % grid positions in mm
pos_grid = [-flipud(pos_mms); pos_mms];  % symmetric position grid
dataFiles = arrayfun(@(d) sprintf('%s/forcevsSpeedDisplacement%.0fMMSameDir.dat', ...
                   baseFolder, d), pos_mms, 'UniformOutput', false);

% % for the case of compressed red MOT
% % Compressed red MOT positional grid (mm)
% pos_mms = [0.01, 0.03, 0.08, 0.15, 0.3, 0.5, 0.7, 0.9, ...
%            1.2, 1.5, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0];
% 
% % Symmetric numeric grid about zero (keep as numeric for computations)
% pos_grid = [-fliplr(pos_mms), pos_mms];
% 
% % Build the numeric strings once
% s_lt1 = arrayfun(@(d) regexprep(sprintf('%.2f', d), '\.?0+$', ''), ...
%                  pos_mms, 'UniformOutput', false);   % e.g., 0.30->0.3, 0.08->0.08
% s_ge1 = arrayfun(@(d) sprintf('%.1f', d), pos_mms, 'UniformOutput', false); % e.g., 2.0
% 
% % Choose formatting based on value
% numStrs = s_ge1;                      % default ≥1: one decimal
% numStrs(pos_mms < 1) = s_lt1(pos_mms < 1);  % <1: strip trailing zeros
% 
% % Build filenames
% dataFiles = cellfun(@(s) sprintf('%s/forcevsSpeedDisplacement%sMMSameDir.dat', ...
%                                  baseFolder, s), ...
%                     numStrs, 'UniformOutput', false);

% Read velocity grid from first file
tbl      = readtable(dataFiles{1});
vel_grid = tbl.Speed;                 % m/s
negRows  = vel_grid < 0;

% Load acceleration and excited-state pop
Nvel = numel(vel_grid);
Npos = numel(pos_mms);
acc_v = zeros(Nvel, Npos);
popP  = zeros(Nvel, Npos);
for i = 1:Npos
    T = readtable(dataFiles{i});
    acc_v(:,i) = T.av;       % mm/ms^2
    popP(:,i)  = T.PFeHigh;  % dimensionless
end

% Flip acceleration sign for negative velocities
acc_v(negRows,:) = -acc_v(negRows,:);

% Symmetrize to ± positions
acc_v_full = [-flipud(fliplr(acc_v)), acc_v];
dupPop     = [flipud(fliplr(popP)),  popP];

%% Interpolants
% Note: divide by 1e3 -> convert accel from mm/ms^2 to m/ms^2 inside interpolant
F_acc_spline = griddedInterpolant({pos_grid, vel_grid}, acc_v_full'/1e3, 'spline');

%% 1) Capture Velocity Determination + Export Selected Trajectories
vs_list    = 0:0.1:50;    % m/s scan
figure('Name','Phase Space Trajectories','NumberTitle','off'); hold on;
plotHandles = []; legendLabels = {};
capture_v  = 0;

for v0 = vs_list
    odeFcn = @(t,p) [p(2); F_acc_spline(p(1),p(2))];
    [~, P] = ode45(odeFcn, [0 50], [min(pos_grid); v0]);  % p1=d(mm), p2=v(m/s)
    d_vals = P(:,1);
    if any(d_vals > cap_threshold) || any(isnan(d_vals)), break; end

    if mod(v0,3)==0
        h = plot(d_vals, P(:,2), 'LineWidth', 2);
        plotHandles(end+1) = h; 
        legendLabels{end+1} = sprintf('v_0 = %.0f m/s', v0); 
    end
    capture_v = v0;
end
capture_v = capture_v - 0.1;  % empirical correction

% --- Only these four export trajectories (per your spec)
extra_v0s = [3.0, 15.0, 28.5, 29.5];

% Prepare to collect all four into one CSV as row pairs:
traj_pos_rows = {};   % each cell: 1 row (1 x Nmax), padded with NaN
traj_vel_rows = {};   % each cell: 1 row (1 x Nmax), padded with NaN
maxLen = 0;

for v0_extra = extra_v0s
    odeFcn = @(t,p) [p(2); F_acc_spline(p(1),p(2))];
    [~, P_extra] = ode45(odeFcn, [0 50], [min(pos_grid); v0_extra]);
    d_extra = P_extra(:,1);  % mm
    v_extra = P_extra(:,2);  % m/s

    % truncate at first threshold/NaN to mirror capture logic
    idx_stop = find((d_extra > sqrt(2)*cap_threshold) | isnan(d_extra), 1, 'first');
    if isempty(idx_stop), idx_plot = 1:numel(d_extra);
    else,                 idx_plot = 1:max(1, idx_stop-1);
    end

    d_use = d_extra(idx_plot);
    v_use = v_extra(idx_plot);

    % Plot for visualization
    h_extra = plot(d_use, v_use, '--', 'LineWidth', 1.5);
    plotHandles(end+1) = h_extra; 
    legendLabels{end+1} = sprintf('v_0 = %.1f m/s', v0_extra);

    % Update max length for padding
    maxLen = max(maxLen, numel(d_use));

    % Store (will pad later). IMPORTANT: v_use(1) is v0 for labeling in Python.
    traj_pos_rows{end+1} = d_use(:).';  
    traj_vel_rows{end+1} = v_use(:).'; 
end

grid off;
legend(plotHandles, legendLabels, 'Location','best');
xlim([min(pos_grid) max(pos_grid)]);
xlabel('d (mm)'); ylabel('v_{d} (m/s)');
title(sprintf('Capture Velocity: %.1f m/s', capture_v));
hold off;

% ---- Pad rows to the same length with NaN and interleave pos/vel rows
nT = numel(traj_pos_rows);
allRows = cell(2*nT, 1);
for j = 1:nT
    posRow = traj_pos_rows{j};
    velRow = traj_vel_rows{j};
    if numel(posRow) < maxLen, posRow = [posRow, nan(1, maxLen - numel(posRow))]; end
    if numel(velRow) < maxLen, velRow = [velRow, nan(1, maxLen - numel(velRow))]; end
    allRows{2*j-1} = posRow;
    allRows{2*j}   = velRow;   % NOTE: velRow(1) == v0, used by your Python for labeling
end
AllCaptureMat = cell2mat(allRows); % (2*nT) x maxLen

% % Export single CSV for all capture trajectories (row-pair format)
% writematrix(AllCaptureMat, fullfile(exportFolder, ['AllCaptureTrajectories', dataStamp, '.csv']));

%% 2) MOT Dynamics Heat Map + Export
xs = min(pos_grid):0.1:max(pos_grid);
vs = min(vel_grid):0.1:max(vel_grid);
[Xg, Vg] = meshgrid(xs, vs);
H = F_acc_spline(Xg', Vg')';   % (length(vs) x length(xs)), still in (mm/ms^2) units for plotting

figure('Name','MOT Dynamics Heat Map','NumberTitle','off');
imagesc(xs, vs, H);
xlim([min(pos_grid) max(pos_grid)]);
ylim([min(vs) max(vs)]);
colorbar;
xlabel('d (mm)'); ylabel('v_{d} (m/s)');
title('Sn MOT Dynamics Heat Map');
g = colorbar; g.Title.String = 'a_{d} (mm/ms^{2})';

% % Export heat map with headers in first row/col and timestamped filename
% hdrRow   = [{'v_ms_rows__d_mm_cols'}, num2cell(xs)];
% firstCol = num2cell(vs');
% body     = num2cell(H);
% cellsOut = [hdrRow; [firstCol, body]];
% writecell(cellsOut, fullfile(exportFolder, ['heatMapData', dataStamp, '.csv']));

%% 3) Acceleration Profiles via Integration + Export
vMax = 8; xMax = 3;
[~, minC] = min(abs(xs + xMax)); [~, maxC] = min(abs(xs - xMax));
[~, minR] = min(abs(vs + vMax)); [~, maxR] = min(abs(vs - vMax));
acc_vs_pos = trapz(vs(minR:maxR), H(minR:maxR,:), 1) ./ (2*vMax);
acc_vs_vel = trapz(xs(minC:maxC), H(:,minC:maxC), 2) ./ (2*xMax);

figure('Name','a vs d','NumberTitle','off');
plot(xs, acc_vs_pos,'LineWidth',2);
xlabel('d (mm)'); ylabel('a_{d} (mm/ms^{2})');
xlim([min(pos_grid) max(pos_grid)]);
title('Acceleration vs Position, Sn redMOT');

figure('Name','a vs v','NumberTitle','off');
plot(vs, acc_vs_vel,'LineWidth',2);
xlabel('v_{d} (m/s)'); ylabel('a_{d} (mm/ms^{2})');
xlim([min(vel_grid) max(vel_grid)]);
title('Acceleration vs Velocity, Sn redMOT');

% % Acceleration vs position: [d_mm, a_d_mm_per_ms2]
% writematrix([xs(:), acc_vs_pos(:)], ...
%     fullfile(exportFolder, ['AccelVsPos', dataStamp, '.csv']));
% 
% % Acceleration vs velocity: [v_ms, a_d_mm_per_ms2]
% writematrix([vs(:), acc_vs_vel(:)], ...
%     fullfile(exportFolder, ['AccelVsVel', dataStamp, '.csv']));

drawnow;

%% 4) MOT Temperature & Size Estimation (+ export one stochastic sample)
if simTrapDynamics
    % Capture region for average pop -> scatter rate
    vMax_trap = vMax; xMax_trap = xMax;
    [~, vrMin] = min(abs(vel_grid + vMax_trap));  [~, vrMax] = min(abs(vel_grid - vMax_trap));
    [~, pMin]  = min(abs(pos_grid + xMax_trap));  [~, pMax]  = min(abs(pos_grid - xMax_trap));

    meanPop     = mean(dupPop(vrMin:vrMax, pMin:pMax), 'all');
    scatterRate = meanPop * gam * 1e-3;   % consistent with dt choice
    dt_kick     = 1/scatterRate;          % ms
    tTotal      = 50;                     % ms runtime
    tWindow     = 10;                     % averaging window (ms)

    % ---- Consistent sizing (avoid off-by-one) ----
    Nsteps  = ceil(tTotal/dt_kick);
    t_eval  = (0:Nsteps)' * dt_kick;      % (Nsteps+1)-by-1

    % Preallocate sliced variables for parfor
    R = zeros(Nsteps+1, nTraj);  % mm
    V = zeros(Nsteps+1, nTraj);  % m/s
    sigma_traj = zeros(nTraj,1);
    temp_traj  = zeros(nTraj,1);

    vKick = hbar * kA / mass;  % m/s per photon

    fprintf('Entering parfor loop with %d trajectories...\n\n', nTraj);
    tic;
    parfor traj = 1:nTraj
        rng(100000+traj, 'twister');
        r = zeros(Nsteps+1,1);
        v = zeros(Nsteps+1,1);
        for k = 1:Nsteps
            a_det = F_acc_spline(r(k), v(k));             % mm/ms^2
            phi1 = 2*pi*rand; phi2 = 2*pi*rand;
            v(k+1) = v(k) + vKick*(cos(phi1)+cos(phi2)) + a_det*dt_kick;  % m/s
            r(k+1) = r(k) + v(k)*dt_kick;                                  % mm
        end
        % Stats over last tWindow ms
        win = t_eval >= (t_eval(end) - tWindow);
        sigma_traj(traj) = sqrt(mean(r(win).^2));     % mm
        v_sigma = sqrt(mean(v(win).^2));              % m/s
        temp_traj(traj)  = mass*(v_sigma^2)/bkB;      % K
        R(:,traj) = r;  V(:,traj) = v;
    end
    elapsedTime = toc;

    sigma_avg_mm = mean(sigma_traj);
    sigma_std_mm = std(sigma_traj);
    temp_avg_mK  = 1e3*mean(temp_traj);
    temp_std_mK  = 1e3*std(temp_traj);

    fprintf('MOT Size (σ): %.2f ± %.2f mm\n\n', sigma_avg_mm, sigma_std_mm);
    fprintf('MOT Temperature: %.1f ± %.1f mK\n\n', temp_avg_mK, temp_std_mK);
    fprintf('Trajectory sim time: %.2f s\n\n', elapsedTime);

    % Plot 
    figure('Name','MOT Sample Trajectories','NumberTitle','off'); hold on;
    for traj=1:nTraj
        plot(R(:,traj), V(:,traj), 'LineWidth',1.5);
    end
    xlabel('d (mm)'); ylabel('v_{d} (m/s)');
    title('MOT Sample Trajectories with Photon Scatter'); hold off;

    % % ---- Export ONE stochastic sample (traj #1) as [d_mm, v_ms] with no header
    % keep = t_eval <= tTotal + eps;   % keep rows up to tTotal (t_eval, R, V are same length)
    % writematrix([R(keep,1), V(keep,1)], ...
    % fullfile(exportFolder, ['InMOTTrajWithScatter', dataStamp, '.csv']));

end

% disp('CSV exports completed:');
% disp(exportFolder);

