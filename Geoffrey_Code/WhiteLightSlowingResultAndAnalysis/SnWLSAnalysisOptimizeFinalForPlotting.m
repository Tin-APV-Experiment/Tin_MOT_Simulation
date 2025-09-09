%---------------------------------------------------------------
% Compute a_z(v_z, r, z) from interpolation of 
% a_z(v_z) through various r and z values as chosen
%from Julia simulation.
% This analysis assumes a converging slowing beam.
% Latest accurate version.
% Only keep P = 300 mW, 4.5 mrad convergence, MOT is 500 mm away
% alternative laser parameters (-11.75 gamma carrier detuning,
% 1.5 gamma and 3.4 gamma EOM modulation frequencies).
%---------------------------------------------------------------

moleculeName = "Sn";
molecule = def_molecule(moleculeName);
bfield = "15";
base_folder = "SnData/";

%for the case of P=300 mW and pure WLS, 4.5 mrad, alt laser parameters,
%z=500 mm (better interpolation)
zList = [0, 100, 200, 300, 400, 500, 750, 1000];
waists = {'5.25', '5.7', '6.15', '6.6', '7.05', '7.5', '8.63', '9.75'};
data_timestamps = {'20250831_2257', '20250831_2303', '20250831_2310', ...
                     '20250831_2316', '20250831_2323', '20250831_2329', '20250831_2337', '20250831_2344'};

% % ---- Export folder at the same level as this script ----
% thisFile  = mfilename('fullpath');
% if isempty(thisFile)
%     scriptDir = pwd;        % running as a section
% else
%     scriptDir = fileparts(thisFile);
% end
% exportFolder = fullfile(scriptDir, 'AllCSVAndMatDataForSn');
% if ~exist(exportFolder, 'dir'), mkdir(exportFolder); end

pos = {'0.01','0.1','0.5','1.0','2.0','3.0','5.0','7.0','9.0','10.0','12.0','14.0','17.0'};
dispList = str2double(pos);
vGrid = linspace(0,300,200);
azInterpTensor = zeros(length(dispList),length(vGrid),length(zList));  % [r,v,z]

for k = 1:length(zList)
    zval = zList(k);
    fprintf("Processing z = %d mm...\n", zval);
    currFolder = sprintf('%s%sWhiteLightSlowingBFieldGauss%sWaistMM%sDate%s', ...
        base_folder, moleculeName, bfield, waists{k}, data_timestamps{k});
    for i = 1:length(dispList)
        currFile = fullfile(currFolder, sprintf('forcevsSpeedDisplacement%sMMRandom.dat',pos{i}));
        if ~isfile(currFile)
            warning('File not found: %s',currFile);
            continue;
        end
        T = readtable(currFile);
        [vs, idx] = sort(T.LongSpeed);
        azs = T.az(idx)/1e3;  % mm/ms^2
        azInterpTensor(i,:,k) = interp1(vs, azs, vGrid, 'spline', 'extrap');
    end
end

%% Plot the 3D acceleration curve surfaces 
%--- Plot surfaces (trimmed to z ≤ 500 mm and v ≤ 225 m/s) ---------------
r_values = [0,3,6];
colors = lines(length(r_values));

% Trim in z and v
zMask        = zList <= 500;
zListPlot    = zList(zMask);
vMask        = vGrid <= 225;
vGrid_plot   = vGrid(vMask);

[VV, ZZ] = meshgrid(vGrid_plot, zListPlot);  % v on X, z on Y
figure; hold on;

for j = 1:length(r_values)
    [~, ridx] = min(abs(dispList - r_values(j)));
    slice_full = squeeze(azInterpTensor(ridx,:,:))';   % [z, v] after transpose
    slice_plot = slice_full(zMask, vMask);             % apply both masks
    surf(VV, ZZ, slice_plot, ...
         'EdgeColor','none','FaceAlpha',0.8, ...
         'DisplayName', sprintf('r=%.1f mm', r_values(j)), ...
         'FaceColor', colors(j,:));
end

xlabel('v_z (m/s)'); ylabel('z (mm)'); zlabel('a_z (mm/ms^2)');
title('a_z(v_z,z) at various r (z <= 500 mm, v <= 225 m/s)');
legend('Location','northeastoutside');
view(135,30); grid on;


% %% Export the 3D acceleration curve surfaces
% 
% %------------------------------------------------------
% % Export a_z(v_z, z) at various r values
% %------------------------------------------------------
% 
% % Export only z ≤ 500 mm
% zMask = zList <= 500;
% zList_trimmed = zList(zMask);
% 
% % Export only v ≤ 225 m/s
% vMask = vGrid <= 225;
% vGrid_trimmed = vGrid(vMask);
% 
% % Trim the acceleration tensor
% azInterpTensor_trimmed = azInterpTensor(:, vMask, zMask);
% 
% % Save the trimmed data -> into exportFolder
% save(fullfile(exportFolder, 'az_vs_z_various_r.mat'), ...
%      'vGrid_trimmed', 'dispList', 'zList_trimmed', 'azInterpTensor_trimmed', ...
%      '-v7');
% 
% fprintf('Saved MAT: %s\n', fullfile(exportFolder, 'az_vs_z_various_r.mat'));

% %% Example plot: a_z(v_z, r) averaged over z ≤ 500 mm and v ≤ 225 m/s
% 
% % Compute mean over z (dim 3), resulting in [dispList × vGrid_trimmed]
% azMeanOverZ = mean(azInterpTensor_trimmed, 3);  % size: [Nr, Nv]
% rMask = dispList <= 10;
% dispList_trimmed = dispList(rMask);
% azMeanOverZ_trimmed = azMeanOverZ(rMask, :);
% 
% % Meshgrid for surface
% [VV, RR] = meshgrid(vGrid_trimmed, dispList_trimmed);  % r on Y, v on X
% 
% figure;
% surf(VV, RR, azMeanOverZ_trimmed, ...
%     'EdgeColor', 'none', ...
%     'FaceAlpha', 0.9);
% xlabel('v_z (m/s)', 'FontSize', 12);
% ylabel('r (mm)', 'FontSize', 12);
% zlabel('a_z (mm/ms²)', 'FontSize', 12);
% title('a_z(v_z, r) Averaged over z ≤ 500 mm, v ≤ 225 m/s (r ≤ 10 mm)', 'FontSize', 14);
% colormap parula;
% view(135, 30);
% grid on;
% 
% % Optionally save the averaged surface data
% save('az_mean_surface_upto_r10mm_v225_final.mat', ...
%      'vGrid_trimmed', 'dispList_trimmed', 'azMeanOverZ_trimmed', '-v7');


%% Monte Carlo 3D trajectory simulation for atom slowing

rng(15248533, "twister");  % reproducible simulation seed

%--- Parameters ------------------------------------------------------------
numTrials   = 1e8;
z_target    = 500;        % mm
zQ          = 7.5;        % mm
v_capture   = 28.5;         % mm/ms
tspan       = [0, 30];    % ms
theta       = 4.5e-3;     % rad

vzMean  = 140; vzSigma = 32;  % mm/ms
vrSigma = 32;                 % mm/ms
R       = 1.5;                % mm
capture_r = 7.0;              % mm

%--- Build interpolant ----------------------------------------------------
Faz = griddedInterpolant({dispList, vGrid, zList}, azInterpTensor, ...
                         'linear','none');
opts = odeset('RelTol',1e-6,'AbsTol',1e-9);

%--- Initial conditions ---------------------------------------------------
vz0 = normrnd(vzMean, vzSigma, numTrials,1);
vx0 = normrnd(0,      vrSigma, numTrials,1);
vy0 = normrnd(0,      vrSigma, numTrials,1);

sigma_r   = vrSigma * zQ / vzMean;
thetaPos  = 2*pi*rand(numTrials,1);
r0_disk   = R * sqrt(rand(numTrials,1));
Ux        = r0_disk .* cos(thetaPos);
Uy        = r0_disk .* sin(thetaPos);
Gx        = normrnd(0, sigma_r, numTrials,1);
Gy        = normrnd(0, sigma_r, numTrials,1);
x0        = Ux + Gx;
y0        = Uy + Gy;

%--- Storage --------------------------------------------------------------
zFinal        = zeros(numTrials,1);
vFinal        = zeros(numTrials,1);
tAtCapture    = NaN(numTrials,1);
reachedTarget = false(numTrials,1);
trajZ         = cell(numTrials,1);  
trajVz        = cell(numTrials,1);
rAtCapture    = NaN(numTrials,1);

%--- Initial‐condition filter ---------------------------------------------
ok_vz   = (vz0 >= 50) & (vz0 <= 200);
xf      = x0 + vx0 * 5;
yf      = y0 + vy0 * 5;
ok_trans= sqrt(xf.^2 + yf.^2) <= capture_r;
ok_init = ok_vz & ok_trans;

%--- Run simulation -------------------------------------------------------
if isempty(gcp('nocreate')), parpool; end

fprintf('Running %d trials after filtering...\n', numTrials);
tic;
parfor i = 1:numTrials
    if ~ok_init(i)
        continue;
    end

    y0i = [zQ; vz0(i); x0(i); vx0(i); y0(i); vy0(i)];
    dydt = @(t,p) dynamicsClamped(t, p, Faz, theta, dispList, zList);
    try
        [tS, yS] = ode45(dydt, tspan, y0i, opts);
    catch
        continue;
    end

    idx = find(yS(:,1) >= z_target, 1, 'first');
    if isempty(idx)
        continue;
    end

    vzAtT = yS(idx,2);
    xAtT  = yS(idx,3);
    yAtT  = yS(idx,5);
    rAtT  = sqrt(xAtT^2 + yAtT^2);

    if vzAtT <= v_capture && rAtT <= capture_r
        reachedTarget(i)   = true;
        tAtCapture(i)      = tS(idx);
        zFinal(i)          = yS(end,1);
        vFinal(i)          = yS(end,2);
        trajZ{i}           = yS(:,1);
        trajVz{i}          = yS(:,2);
        rAtCapture(i)      = rAtT;
    end
end
elapsedTime = toc;

numReached = sum(reachedTarget);
fprintf('Simulation completed in %.2f s.\n', elapsedTime);
fprintf('Captured at z = %.1f mm (v_z ≤ %.1f mm/ms and r ≤ %.1f mm): %d / %d (%.2e)\n', ...
    z_target, v_capture, capture_r, numReached, numTrials, numReached/numTrials);
fprintf('Avg capture time: %.3f ms\n', mean(tAtCapture(reachedTarget), 'omitnan'));

%% Plotting
% Random 5% of successes
successIdx = find(reachedTarget);
rng(92173640, "twister");
nToPlot = max(1, round(0.05 * numel(successIdx)));
plotIdx = successIdx(randperm(numel(successIdx), nToPlot));

figure; hold on;
hSucc = gobjects(0);
for k = 1:numel(plotIdx)
    j = plotIdx(k);
    ztraj  = trajZ{j};
    vztraj = trajVz{j};
    mask   = (ztraj >= zQ) & (ztraj <= 1000);
    h = plot(ztraj(mask), vztraj(mask), 'Color',[0 0.6 0],'LineWidth',1);
    if k == 1, hSucc = h; end
end

% Obtain every 80,000th unsuccessful trajectory
sampStep = 8e4;
sampleCandidates = sampStep:sampStep:numTrials;
unsuccMask = (~ok_init(sampleCandidates)) | (ok_init(sampleCandidates) & ~reachedTarget(sampleCandidates));
unsuccIdx  = sampleCandidates(unsuccMask);

unsuccZ  = cell(numel(unsuccIdx),1);
unsuccVz = cell(numel(unsuccIdx),1);
for m = 1:numel(unsuccIdx)
    i = unsuccIdx(m);
    y0i = [zQ; vz0(i); x0(i); vx0(i); y0(i); vy0(i)];
    dydt = @(t,p) dynamicsClamped(t, p, Faz, theta, dispList, zList);
    try
        [tU, yU] = ode45(dydt, tspan, y0i, opts);
        unsuccZ{m}  = yU(:,1);
        unsuccVz{m} = yU(:,2);
    catch
        unsuccZ{m}  = [];
        unsuccVz{m} = [];
    end
end

% %Plot every 80000th unsuccessful trajectory if so desired.
% hUnsucc = gobjects(0);
% for m = 1:numel(unsuccZ)
%     if ~isempty(unsuccZ{m})
%         ztraj  = unsuccZ{m};
%         vztraj = unsuccVz{m};
%         mask   = (ztraj >= zQ) & (ztraj <= 1000);
%         h = plot(ztraj(mask), vztraj(mask), 'Color',[0.7 0 0], 'LineWidth',0.8);
%         if m == 1, hUnsucc = h; end
%     end
% end

hx = xline(z_target,  'k--','LineWidth',1.2);
hy = yline(v_capture, 'k--','LineWidth',1.2);
% legend([hSucc hx hy], ...
%        {'Success (10% sample)','z target','v capture'}, ...
%        'Location','northeastoutside');
xlabel('z (mm)'); ylabel('v_z (mm/ms)');


% %% Export longitudinal trajectories to CSV 
% % Assumes 'exportFolder' already set to .../AllCSVAndMatDataForSn
% 
% rowsSuccessLong = {};
% rowsFailLong    = {};
% 
% % Successful trajectories
% for k = 1:numel(plotIdx)
%     j = plotIdx(k);
%     ztraj  = trajZ{j};
%     vztraj = trajVz{j};
%     if isempty(ztraj), continue; end
%     ztraj  = ztraj(:);              % ensure column
%     vztraj = vztraj(:);             % ensure column
%     nPts   = numel(ztraj);
%     rowsSuccessLong{end+1} = table( ...
%         repmat(j, nPts, 1), (1:nPts)', ztraj, vztraj, ...
%         'VariableNames', {'trajID','step','z','vz'} );
% end
% 
% % Unsuccessful trajectories
% for m = 1:numel(unsuccIdx)
%     ztraj  = unsuccZ{m};
%     vztraj = unsuccVz{m};
%     if isempty(ztraj), continue; end
%     ztraj  = ztraj(:);
%     vztraj = vztraj(:);
%     nPts   = numel(ztraj);
%     rowsFailLong{end+1} = table( ...
%         repmat(unsuccIdx(m), nPts, 1), (1:nPts)', ztraj, vztraj, ...
%         'VariableNames', {'trajID','step','z','vz'} );
% end
% 
% timestamp = datestr(now, 'yyyymmdd_HHMM');
% success_long_filename = fullfile(exportFolder, ['successful_trajectories_longitudinal_' timestamp '.csv']);
% fail_long_filename    = fullfile(exportFolder, ['unsuccessful_trajectories_longitudinal_' timestamp '.csv']);
% 
% if ~isempty(rowsSuccessLong)
%     T_success_long = vertcat(rowsSuccessLong{:});
%     writetable(T_success_long, success_long_filename);
%     fprintf('Saved %s (%d rows)\n', success_long_filename, height(T_success_long));
% else
%     fprintf('No successful longitudinal trajectories to export.\n');
% end
% 
% if ~isempty(rowsFailLong)
%     T_fail_long = vertcat(rowsFailLong{:});
%     writetable(T_fail_long, fail_long_filename);
%     fprintf('Saved %s (%d rows)\n', fail_long_filename, height(T_fail_long));
% else
%     fprintf('No unsuccessful longitudinal trajectories to export.\n');
% end


%% Dynamics function
function dpdt = dynamicsClamped(~, p, Faz, theta, rAxis, zAxis)
    z  = p(1); vz = p(2);
    x  = p(3); vx = p(4);
    y  = p(5); vy = p(6);
    r  = sqrt(x^2 + y^2);
    if r >= rAxis(1) && r <= rAxis(end) && z >= zAxis(1) && z <= zAxis(end)
        az = Faz(r, vz, z);
        ax = az * tan(theta) * (x/(r + eps));
        ay = az * tan(theta) * (y/(r + eps));
    else
        az = 0; ax = 0; ay = 0;
    end
    dpdt = [vz; az; vx; ax; vy; ay];
end