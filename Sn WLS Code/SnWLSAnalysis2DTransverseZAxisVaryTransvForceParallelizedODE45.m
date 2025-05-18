moleculeName = "Sn";
molecule = def_molecule(moleculeName);
bfield = "15";
base_folder = "SnData/";

% Define z positions and associated file details
% zList = [0, 100, 200, 300, 400, 500];
% waists = {'5.25', '5.7', '6.15', '6.6', '7.05', '7.5'};
% data_timestamps = {'20250506_1355', '20250506_1349', '20250506_1343', ...
%                    '20250506_1335', '20250506_1328', '20250506_1316'};
% data_timestamps = {'20250502_1405', '20250502_1401', '20250502_1357', ...
%                    '20250502_1352', '20250502_1344', '20250430_1514'};


%----------------------------------------------------
% Compute a_z(v_z, r, z) from interpolation of 
% a_z(v_z) through various r and z values as chosen
% from Julia simulation.
%----------------------------------------------------


% %for the case of P=170 mW and pure WLS
% zList = [0, 100, 200, 300, 400, 500, 750, 1000];
% waists = {'5.25', '5.7', '6.15', '6.6', '7.05', '7.5', '8.63', '9.75'};
% data_timestamps = {'20250506_1644', '20250506_1648', '20250506_1653', ...
%                      '20250506_1659', '20250506_1704', '20250506_1708', '20250506_1713', '20250506_1718'};

%for the case of P=170 mW and push beam
zList = [0, 100, 200, 300, 400, 500, 750, 1000];
waists = {'5.25', '5.7', '6.15', '6.6', '7.05', '7.5', '8.63', '9.75'};
data_timestamps = {'20250515_1544', '20250515_1537', '20250515_1531', ...
                     '20250515_1519', '20250515_1512', '20250515_1505', '20250515_1459', '20250515_1453'};


% Define r displacements and v grid
pos = {'0.01', '0.1', '0.5', '1.0', '2.0', '3.0', '5.0', '7.0', '9.0', '10.0', '12.0', '14.0', '17.0'};
dispList = str2double(pos);
vGrid = linspace(0, 300, 200);
azInterpTensor = zeros(length(dispList), length(vGrid), length(zList));

% Loop over z positions and displacements
for k = 1:length(zList)
    zval = zList(k);
    fprintf("Processing z = %d mm...\n", zval);

    currFolder = strcat(base_folder, moleculeName, ...
        "WhiteLightSlowingBFieldGauss", bfield, ...
        "WaistMM", waists{k}, ...
        "Date", data_timestamps{k});

    for i = 1:length(pos)
        currFile = strcat(currFolder, '/forcevsSpeedDisplacement', pos{i}, 'MMRandom.dat');
        if ~isfile(currFile)
            warning('File not found: %s', currFile);
            continue;
        end

        currData = readtable(currFile);
        [LongSpeedsSorted, sortIdx] = sort(currData.LongSpeed);
        azsSorted = currData.az(sortIdx) / 1e3;  % µm/ms² → mm/ms²

        azInterpTensor(i, :, k) = interp1(LongSpeedsSorted, azsSorted, vGrid, 'spline', 'extrap');
    end
end

% Define 3D interpolation function
azInterpFcn3D = @(vz, r, z) interp3(vGrid, dispList, zList, azInterpTensor, vz, r, z, 'spline');

% Define desired r values to plot
r_values = [0, 3, 6];  % mm
colors = lines(length(r_values));  % For distinguishing each surface

% Generate meshgrid for velocity and z

%[VV, ZZ] = meshgrid(vGrid, zList);  % VV: v_z, ZZ: z
zListPlot = zList(zList <= 500);
[VV, ZZ] = meshgrid(vGrid, zListPlot);  % VV: v_z, ZZ: z
% Open figure
figure;
hold on;

% Loop through each r value
for idx = 1:length(r_values)
    r_val = r_values(idx);
    [~, r_idx] = min(abs(dispList - r_val));  % Find closest index

    %az_slice = squeeze(azInterpTensor(r_idx, :, :))';  % [z × v]
    az_slice_full = squeeze(azInterpTensor(r_idx, :, :))';  % size: [length(zList), length(vGrid)]
    az_slice = az_slice_full(zList <= 500, :);  % keep only z ≤ 500

    % Plot surface with partial transparency for clarity
    surf(VV, ZZ, az_slice, ...
         'EdgeColor', 'none', ...
         'FaceAlpha', 0.8, ...
         'DisplayName', sprintf('r = %.1f mm', r_val), ...
         'FaceColor', colors(idx, :));
end

% Labels and appearance
xlabel('v_z (m/s)', 'FontSize', 12);
ylabel('z (mm)', 'FontSize', 12);
zlabel('a_z (mm/ms^2)', 'FontSize', 12);
title('a_z(v_z, z) at Various r', 'FontSize', 14);
legend('Location', 'northeastoutside');
colormap parula;
view(135, 30);
grid on;
hold off;

% % Save vGrid, zList, r_values, and the slices of azInterpTensor for desired r
% save('az_data.mat', 'vGrid', 'zList', 'r_values', 'dispList', 'azInterpTensor');

% %------------------------------------------------------
% % Export a_z(v_z, r) averaged over z-value of beamline
% % a_z(v_z, z) at various r values
% %------------------------------------------------------
% 
% % Export only z ≤ 500 mm data for use in Python
% zMask = zList <= 500;
% zList_trimmed = zList(zMask);
% azInterpTensor_trimmed = azInterpTensor(:, :, zMask);
% 
% % Compute mean over z (dim 3), resulting in [dispList × vGrid]
% azMeanOverZ = mean(azInterpTensor_trimmed, 3);  % size: [Nr, Nv]
% rMask = dispList <= 10;
% dispList_trimmed = dispList(rMask);
% azMeanOverZ_trimmed = azMeanOverZ(rMask, :);
% 
% 
% 
% save('az_tensor_upto_500mm.mat', ...
%      'vGrid', 'dispList', 'zList_trimmed', 'azInterpTensor_trimmed', ...
%      '-v7');
% 
% % Meshgrid for surface
% [VV, RR] = meshgrid(vGrid, dispList_trimmed);  % r on Y, v on X
% 
% figure;
% surf(VV, RR, azMeanOverZ_trimmed, ...
%     'EdgeColor', 'none', ...
%     'FaceAlpha', 0.9);
% xlabel('v_z (m/s)', 'FontSize', 12);
% ylabel('r (mm)', 'FontSize', 12);
% zlabel('a_z (mm/ms²)', 'FontSize', 12);
% title('a_z(v_z, r) Averaged over z ≤ 500 mm (r ≤ 10 mm)', 'FontSize', 14);
% colormap parula;
% view(135, 30);
% grid on;
% 
% save('az_mean_surface_upto_r10mm.mat', ...
%      'vGrid', 'dispList_trimmed', 'azMeanOverZ_trimmed', '-v7');

% figure;
% 
% vTest = linspace(0, 200, 100);
% aTest = arrayfun(@(v) azInterpFcn3D(v, 25, 247), vTest);
% plot(vTest, aTest); xlabel('v_z'); ylabel('a_z'); title('Acceleration vs Velocity');

% %----------------------------------------------------------------------
% % Parallelized Monte Carlo trajectory simulation
% % Using a_z(v_z, r, z) with ode45
% % Includes transverse displacement + velocity
% % Include transverse confining force
% % Clamps a_z = 0 outside max transverse disp. and max z-value in sims.
% %----------------------------------------------------------------------
% 
% 
% % Parameters
% numTrials = 10000;
% z_target = 500;           % mm
% v_capture = 27;           % mm/ms
% tspan = [0, 30];          % ms
% initialZ = 0;             % mm
% vzMean = 140; vzSigma = 32; % mm/ms
% vrSigma = 32;             % mm/ms
% rInitMax = 1.5;           % mm
% capture_r = 7.0;          % mm
% theta = 4.5e-3;           % rad (convergence half-angle)
% 
% % Generate initial conditions
% vz0 = normrnd(vzMean, vzSigma, numTrials, 1);
% vx0 = normrnd(0, vrSigma, numTrials, 1);
% vy0 = normrnd(0, vrSigma, numTrials, 1);
% 
% thetaPos = 2 * pi * rand(numTrials, 1);
% radii = rInitMax * sqrt(rand(numTrials, 1));
% x0 = radii .* cos(thetaPos);
% y0 = radii .* sin(thetaPos);
% 
% % Preallocate
% zFinal = zeros(numTrials, 1);
% vFinal = zeros(numTrials, 1);
% tAtCapture = NaN(numTrials, 1);
% reachedTarget = false(numTrials, 1);
% trajZ = cell(numTrials, 1);
% trajVz = cell(numTrials, 1);
% trajX = cell(numTrials, 1);
% trajY = cell(numTrials, 1);
% 
% maxZ = max(zList);
% maxR = max(dispList);
% 
% fprintf('Running %d Monte Carlo trials with parfor + clamped a_z for r>r_max and z>z_max...\n', numTrials);
% 
% if isempty(gcp('nocreate'))
%     parpool;
% end
% 
% tic;
% parfor i = 1:numTrials
%     % Initial state: [z; vz; x; vx; y; vy]
%     y0i = [initialZ; vz0(i); x0(i); vx0(i); y0(i); vy0(i)];
% 
%     % Interpolated az with clamping
%     azVal = @(vz, r, z) double(r <= maxR && z <= maxZ) * ...
%         interp3(vGrid, dispList, zList, azInterpTensor, vz, r, z, 'spline');
% 
%     % ODE function
%     odefun = @(t, p) [
%         p(2);
%         azVal(p(2), sqrt(p(3)^2 + p(5)^2), p(1));
%         p(4);
%         azVal(p(2), sqrt(p(3)^2 + p(5)^2), p(1)) * tan(theta) * (p(3) / sqrt(p(3)^2 + p(5)^2 + eps));
%         p(6);
%         azVal(p(2), sqrt(p(3)^2 + p(5)^2), p(1)) * tan(theta) * (p(5) / sqrt(p(3)^2 + p(5)^2 + eps));
%     ];
% 
%     try
%         [tSol, ySol] = ode45(odefun, tspan, y0i);
%     catch
%         continue;
%     end
% 
%     zFinal(i) = ySol(end, 1);
%     vFinal(i) = ySol(end, 2);
% 
%     crossingIdx = find(ySol(:, 1) >= z_target, 1, 'first');
%     if ~isempty(crossingIdx)
%         t_cross = tSol(crossingIdx);
%         vzAtTarget = ySol(crossingIdx, 2);
%         r_target = sqrt(ySol(crossingIdx,3)^2 + ySol(crossingIdx,5)^2);
% 
%         if vzAtTarget <= v_capture && r_target <= capture_r
%             reachedTarget(i) = true;
%             tAtCapture(i) = t_cross;
%             trajZ{i} = ySol(:,1);
%             trajVz{i} = ySol(:,2);
%             trajX{i} = ySol(:,3);
%             trajY{i} = ySol(:,5);
%         end
%     end
% end
% 
% elapsedTime = toc;
% fprintf('Simulation completed in %.2f seconds.\n', elapsedTime);
% 
% numReached = sum(reachedTarget);
% fprintf('Captured at z = %.1f mm (v_z <= %.1f mm/ms and |r| < %.1f mm): %d / %d (%.2e)\n', ...
%     z_target, v_capture, capture_r, numReached, numTrials, numReached/numTrials);
% 
% avgCaptureTime = mean(tAtCapture(reachedTarget), 'omitnan');
% fprintf('Avg capture time: %.3f ms\n', avgCaptureTime);
% 
% % Plot successful trajectories
% figure;
% hold on;
% for i = 1:numTrials
%     if reachedTarget(i)
%         plot(trajZ{i}, trajVz{i}, 'Color', [0, 0.6, 0], 'LineWidth', 1.0);
%     end
% end
% xline(z_target, 'k--', 'LineWidth', 1.2);
% yline(v_capture, 'k--', 'LineWidth', 1.2);
% xlabel('z (mm)', 'FontSize', 12);
% ylabel('v_z (mm/ms)', 'FontSize', 12);
% title('Successful Monte Carlo Trajectories: v_z vs z', 'FontSize', 14);
% grid on;


%----------------------------------------------------
% Monte Carlo 3D trajectory simulation for atom slowing
% Includes longitudinal and transverse deceleration
% Uses ode45 for integration
%----------------------------------------------------

rng(12345);  % Fixed seed for reproducibility

% Parameters
numTrials = 2500000;
z_target = 500;           % mm
v_capture = 27;           % mm/ms
tspan = [0, 30];          % ms
initialZ = 0;             % mm
vzMean = 140; vzSigma = 32; % mm/ms
vrSigma = 32;             % mm/ms
rInitMax = 1.5;           % mm
capture_r = 7.0;          % mm
theta = 4.5e-3;           % rad (convergence half-angle)

% Generate initial conditions
vz0 = normrnd(vzMean, vzSigma, numTrials, 1);
vx0 = normrnd(0, vrSigma, numTrials, 1);
vy0 = normrnd(0, vrSigma, numTrials, 1);

thetaPos = 2 * pi * rand(numTrials, 1);
radii = rInitMax * sqrt(rand(numTrials, 1));
x0 = radii .* cos(thetaPos);
y0 = radii .* sin(thetaPos);

% Preallocate
zFinal = zeros(numTrials, 1);
vFinal = zeros(numTrials, 1);
tAtCapture = NaN(numTrials, 1);
reachedTarget = false(numTrials, 1);
trajZ = cell(numTrials, 1);
trajVz = cell(numTrials, 1);
trajX = cell(numTrials, 1);
trajY = cell(numTrials, 1);

maxZ = max(zList);
maxR = max(dispList);

% --- Local dynamics function used in parfor/ode45 ---
function dpdt = dynamics(~, p, vGrid, dispList, zList, azInterpTensor, theta, maxZ, maxR)
    z = p(1); vz = p(2); x = p(3); vx = p(4); y = p(5); vy = p(6);
    r = sqrt(x^2 + y^2);

    az = 0;
    if r <= maxR && z <= maxZ
        az = interp3(vGrid, dispList, zList, azInterpTensor, vz, r, z, 'spline');
    end

    ax = az * tan(theta) * (x / (r + eps));
    ay = az * tan(theta) * (y / (r + eps));

    dpdt = [vz; az; vx; ax; vy; ay];
end


fprintf('Running %d Monte Carlo trials with parfor + ode45 integrator...\n', numTrials);

if isempty(gcp('nocreate'))
    parpool;
end

tic;
parfor i = 1:numTrials
    y0i = [initialZ; vz0(i); x0(i); vx0(i); y0(i); vy0(i)];

    % Define ODE system as anonymous function calling shared function
    dydt = @(t, p) dynamics(t, p, vGrid, dispList, zList, azInterpTensor, theta, maxZ, maxR);
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

    try
        [tSol, ySol] = ode45(dydt, tspan, y0i, opts);
    catch
        % in case of solver failure, mark as not reaching target
        continue;
    end

    zFinal(i) = ySol(end, 1);
    vFinal(i) = ySol(end, 2);

    crossingIdx = find(ySol(:, 1) >= z_target, 1, 'first');
    if ~isempty(crossingIdx)
        t_cross = tSol(crossingIdx);
        vzAtTarget = ySol(crossingIdx, 2);
        r_target = sqrt(ySol(crossingIdx,3)^2 + ySol(crossingIdx,5)^2);

        trajZ{i} = ySol(:,1);
        trajVz{i} = ySol(:,2);
        trajX{i} = ySol(:,3);
        trajY{i} = ySol(:,5);

        if vzAtTarget <= v_capture && r_target <= capture_r
            reachedTarget(i) = true;
            tAtCapture(i) = t_cross;
            % trajZ{i} = ySol(:,1);
            % trajVz{i} = ySol(:,2);
            % trajX{i} = ySol(:,3);
            % trajY{i} = ySol(:,5);
        end
    end
end

elapsedTime = toc;
fprintf('Simulation completed in %.2f seconds.\n', elapsedTime);

numReached = sum(reachedTarget);
fprintf('Captured at z = %.1f mm (v_z <= %.1f mm/ms and |r| < %.1f mm): %d / %d (%.2e)\n', ...
    z_target, v_capture, capture_r, numReached, numTrials, numReached/numTrials);

avgCaptureTime = mean(tAtCapture(reachedTarget), 'omitnan');
fprintf('Avg capture time: %.3f ms\n', avgCaptureTime);

% Plot successful trajectories
figure;
hold on;
for i = 1:numTrials
    if reachedTarget(i)
        plot(trajZ{i}, trajVz{i}, 'Color', [0, 0.6, 0], 'LineWidth', 1.0);
    end
end
xline(z_target, 'k--', 'LineWidth', 1.2);
yline(v_capture, 'k--', 'LineWidth', 1.2);
xlabel('z (mm)', 'FontSize', 12);
ylabel('v_z (mm/ms)', 'FontSize', 12);
title('Successful Monte Carlo Trajectories: v_z vs z', 'FontSize', 14);
grid on;

% Combine and export all trajectories with r = sqrt(x^2 + y^2)
rowsSuccess = {};
rowsFail = {};
rowsSuccessLong = {};
rowsFailLong = {};

for i = 1:numTrials
    trajData = table();        % full trajectory (z, vz, x, y, r)
    trajDataLong = table();    % only z and vz

    if ~isempty(trajZ{i})
        nPts = numel(trajZ{i});
        trajID = repmat(i, nPts, 1);

        % Compute r = sqrt(x^2 + y^2)
        rVals = sqrt(trajX{i}.^2 + trajY{i}.^2);

        % Full trajectory
        trajData = table( ...
            trajID, ...
            (1:nPts)', ...
            trajZ{i}, ...
            trajVz{i}, ...
            trajX{i}, ...
            trajY{i}, ...
            rVals, ...
            'VariableNames', {'trajID', 'step', 'z', 'vz', 'x', 'y', 'r'} ...
        );

        % Longitudinal only
        trajDataLong = table( ...
            trajID, ...
            (1:nPts)', ...
            trajZ{i}, ...
            trajVz{i}, ...
            'VariableNames', {'trajID', 'step', 'z', 'vz'} ...
        );
    end

    if reachedTarget(i)
        if ~isempty(trajData)
            rowsSuccess{end+1} = trajData;
            rowsSuccessLong{end+1} = trajDataLong;
        end
    else
        % Export only every 1000th unsuccessful trajectory
        if ~isempty(trajData) && mod(i, 1000) == 0
            rowsFail{end+1} = trajData;
            rowsFailLong{end+1} = trajDataLong;
        end
    end
end

% Timestamp to avoid overwriting
timestamp = datestr(now, 'yyyymmdd_HHMM');
success_filename       = ['successful_trajectories_with_r_' timestamp '.csv'];
fail_filename          = ['unsuccessful_trajectories_with_r_' timestamp '.csv'];
success_long_filename  = ['successful_trajectories_longitudinal_' timestamp '.csv'];
fail_long_filename     = ['unsuccessful_trajectories_longitudinal_' timestamp '.csv'];

% Write to CSV
if ~isempty(rowsSuccess)
    T_success = vertcat(rowsSuccess{:});
    writetable(T_success, success_filename);
    fprintf('Saved %s\n', success_filename);

    T_success_long = vertcat(rowsSuccessLong{:});
    writetable(T_success_long, success_long_filename);
    fprintf('Saved %s\n', success_long_filename);
end

if ~isempty(rowsFail)
    T_fail = vertcat(rowsFail{:});
    writetable(T_fail, fail_filename);
    fprintf('Saved %s\n', fail_filename);

    T_fail_long = vertcat(rowsFailLong{:});
    writetable(T_fail_long, fail_long_filename);
    fprintf('Saved %s\n', fail_long_filename);
end



% %----------------------------------------------------
% % 1D Monte Carlo transverse trajectory simulation using a_perp = tan(theta) * a_z
% % Uses ode45 to compute x(t) for many trials with random x0, vx0
% % Evaluates probability that final x lies within capture_r at z = z_target
% %----------------------------------------------------
% 
% % Parameters
% numTrials = 1e4;
% z_target = 500;            % mm
% capture_r = 7.0;           % mm
% tspan = [0, 30];           % ms
% theta = 4.5e-3;            % rad
% vrSigma = 32;              % mm/ms (transverse vel std dev)
% rInitMax = 1.5;            % mm
% vzMean = 140;              % mm/ms
% 
% % Initial conditions
% x0 = -rInitMax + 2*rInitMax * rand(numTrials, 1);
% vx0 = normrnd(0, vrSigma, numTrials, 1);
% 
% % Preallocate
% xFinal = zeros(numTrials, 1);
% reachedCapture = false(numTrials, 1);
% trajZ = cell(numTrials, 1);
% trajX = cell(numTrials, 1);
% 
% maxZ = max(zList);
% maxR = max(dispList);
% 
% fprintf("Running 1D transverse Monte Carlo with ode45, %d trials...\n", numTrials);
% 
% if isempty(gcp('nocreate'))
%     parpool;
% end
% 
% tic;
% parfor i = 1:numTrials
%     % State vector: [x; vx; z]
%     y0 = [x0(i); vx0(i); 0];
% 
%     odefun = @(t, y) [
%         y(2);  % dx/dt = vx
%         double(abs(y(1)) <= maxR && y(3) <= maxZ) * ...
%         interp3(vGrid, dispList, zList, azInterpTensor, vzMean, abs(y(1)), y(3), 'spline') * tan(theta) * sign(y(1));
%         vzMean  % dz/dt = constant
%     ];
% 
%     try
%         [tSol, ySol] = ode45(odefun, tspan, y0);
%     catch
%         continue;
%     end
% 
%     % Interpolate x position at z_target
%     zTraj = ySol(:, 3);
%     xTraj = ySol(:, 1);
%     if any(zTraj >= z_target)
%         idx = find(zTraj >= z_target, 1, 'first');
%         xAtTarget = xTraj(idx);
%         xFinal(i) = xAtTarget;
%         if abs(xAtTarget) <= capture_r
%             reachedCapture(i) = true;
%             trajX{i} = xTraj;
%             trajZ{i} = zTraj;
%         end
%     end
% end
% 
% tElapsed = toc;
% fprintf("Simulation completed in %.2f seconds.\n", tElapsed);
% 
% probCapture = sum(reachedCapture) / numTrials;
% fprintf("Probability of capture in transverse dimension: %.4f\n", probCapture);
% 
% % % Optional histogram
% % figure;
% % histogram(xFinal(reachedCapture), 'BinWidth', 0.25);
% % xlabel('Final x (mm)'); ylabel('Count');
% % title('Histogram of Final x-positions (Captured)');
% % grid on;
% 
% % Plot successful transverse trajectories
% figure;
% hold on;
% for i = 1:numTrials
%     if reachedCapture(i)
%         plot(trajZ{i}, trajX{i}, 'Color', [0, 0.6, 0], 'LineWidth', 1.0);
%     end
% end
% xline(z_target, 'k--', 'LineWidth', 1.2);
% yline(capture_r, 'k--', 'LineWidth', 1.2);
% yline(-capture_r, 'k--', 'LineWidth', 1.2);
% xlabel('z (mm)', 'FontSize', 12);
% ylabel('x (mm)', 'FontSize', 12);
% xlim([0, 700]);
% title('Successful Monte Carlo Trajectories: x vs z', 'FontSize', 14);
% grid on;


% %----------------------------------------------------
% % 1D Monte Carlo transverse trajectory simulation with no transverse deceleration (collimated beam)
% % Uses ode45 to compute x(t) for many trials with random x0, vx0
% % Evaluates probability that final x lies within capture_r at z = z_target
% %----------------------------------------------------
% 
% % Parameters
% numTrials = 1e5;
% z_target = 500;            % mm
% capture_r = 7.0;           % mm
% tspan = [0, 30];           % ms
% theta = 4.5e-3;            % rad (unused in collimated case)
% vrSigma = 32;              % mm/ms (transverse vel std dev)
% rInitMax = 1.5;            % mm
% vzMean = 140;              % mm/ms
% 
% % Initial conditions
% x0 = -rInitMax + 2*rInitMax * rand(numTrials, 1);
% vx0 = normrnd(0, vrSigma, numTrials, 1);
% 
% % Preallocate
% xFinal = zeros(numTrials, 1);
% reachedCapture = false(numTrials, 1);
% trajZ = cell(numTrials, 1);
% trajX = cell(numTrials, 1);
% 
% fprintf("Running 1D transverse Monte Carlo with ode45 (collimated), %d trials...\n", numTrials);
% 
% if isempty(gcp('nocreate'))
%     parpool;
% end
% 
% tic;
% parfor i = 1:numTrials
%     % State vector: [x; vx; z]
%     y0 = [x0(i); vx0(i); 0];
% 
%     odefun = @(t, y) [
%         y(2);     % dx/dt = vx
%         0;        % no transverse acceleration
%         vzMean    % dz/dt = constant
%     ];
% 
%     try
%         [tSol, ySol] = ode45(odefun, tspan, y0);
%     catch
%         continue;
%     end
% 
%     % Interpolate x position at z_target
%     zTraj = ySol(:, 3);
%     xTraj = ySol(:, 1);
%     if any(zTraj >= z_target)
%         idx = find(zTraj >= z_target, 1, 'first');
%         xAtTarget = xTraj(idx);
%         xFinal(i) = xAtTarget;
%         if abs(xAtTarget) <= capture_r
%             reachedCapture(i) = true;
%             trajX{i} = xTraj;
%             trajZ{i} = zTraj;
%         end
%     end
% end
% 
% tElapsed = toc;
% fprintf("Simulation completed in %.2f seconds.\n", tElapsed);
% 
% probCapture = sum(reachedCapture) / numTrials;
% fprintf("Probability of capture in transverse dimension (collimated): %.4f\n", probCapture);
% 
% % % Optional histogram
% % figure;
% % histogram(xFinal(reachedCapture), 'BinWidth', 0.25);
% % xlabel('Final x (mm)'); ylabel('Count');
% % title('Histogram of Final x-positions (Captured, Collimated Beam)');
% % grid on;
% 
% % Plot successful transverse trajectories
% figure;
% hold on;
% for i = 1:numTrials
%     if reachedCapture(i)
%         plot(trajZ{i}, trajX{i}, 'Color', [0, 0.4, 0.9], 'LineWidth', 1.0);
%     end
% end
% xline(z_target, 'k--', 'LineWidth', 1.2);
% yline(capture_r, 'k--', 'LineWidth', 1.2);
% yline(-capture_r, 'k--', 'LineWidth', 1.2);
% xlabel('z (mm)', 'FontSize', 12);
% ylabel('x (mm)', 'FontSize', 12);
% title('Successful Transverse Trajectories: x vs z (Collimated)', 'FontSize', 14);
% xlim([0, 700]);
% grid on;