colors = [[0.5;0.2;0.8],[0.4;0.7;0.8],[0.8;0.5;0.2],[0.5;0.7;0.2],...
    [0.6;0.6;0.3],[0.6;0.3;0.6],[0.3;0.6;0.6],[0.4;0.4;0.8],[0.4;0.8;0.4],[0.8;0.4;0.4],[0.5;0.5;0.5],[0.9;0.2;0.1],...
    [0.2;0.4;0.7],[0.4;0.4;0.0],[0.0;0.4;0.2],[0.2;0.0;0.4],[0.0;0.0;0.0],[0.8;0.8;0.8]];

moleculeName = "Sn";
molecule = def_molecule(moleculeName);
bfield = "15"; %static b-field used
num_lasers_used = "1"; %number of lasers used in WLS. 1 for WLS no push, 2 for WLS w/push
data_timestamp = "20250518_1224"; %collimated slowing beam and push beam
%data_timestamp = "20250429_1913"; %just collimated slowing beam
base_folder = "SnData/";
dataFolder = strcat(base_folder, moleculeName, "WhiteLightSlowingBFieldGauss", bfield, "WaistMM7", "Date", data_timestamp);


% %---------------------------------------------
% % obtain curves for a_z vs v_z, v_z vs z, etc
% %---------------------------------------------
% 
% 
% %pos = {'0.01', '0.1', '0.5', '1.0', '2.0', '3.0', '5.0', '7.0', '9.0', '10.0'};
% pos = {'0.01'};
% 
% for i = 1:length(pos)
%     currDisp = pos{i};
%     fprintf('Processing displacement = %s mm\n', currDisp);
% 
%     % Read data
%     currFile = strcat(dataFolder, '/forcevsSpeedDisplacement', currDisp, 'MMRandom.dat');
%     currData = readtable(currFile);
% 
%     LongSpeeds = currData.LongSpeed;
%     azs = currData.az / 1e3; % Convert from µm/ms² to mm/ms²
%     excitedPop = currData.PFeHigh;
% 
%     % Sort data for interpolation
%     [LongSpeedsSorted, sortIdx] = sort(LongSpeeds);
%     azsSorted = azs(sortIdx);
% 
%     % 1D interpolation
%     LongSpeedsInterp = linspace(min(LongSpeedsSorted), max(LongSpeedsSorted), 200);
%     azsInterp = interp1(LongSpeedsSorted, azsSorted, LongSpeedsInterp, 'spline');
% 
%     % Plot acceleration vs velocity
%     figure(1);
%     plot(LongSpeedsInterp, azsInterp, '-', 'LineWidth', 2, 'DisplayName', ['disp = ' currDisp ' mm']);
%     hold on;
%     xlabel('v_z (m/s)');
%     ylabel('a_z(v_z) (mm/ms^2)');
%     title('White Light Slowing of Sn Atoms');
%     legend('show', 'Location', 'northeastoutside');
%     grid on;
% 
%     % Define interpolated acceleration function
%     slowingAccel = @(v) interp1(LongSpeedsSorted, azsSorted, v, 'spline');
% 
%     % Initial conditions
%     tspan = [0 20]; % Time in ms
%     initialPos = 0.0; % mm
%     initialVel = 140; % mm/ms (equivalent to m/s)
%     slowingLength = 400.0; % mm
%     y0 = [initialPos; initialVel];
% 
%     % ODE simulation
%     odefun = @(t, p) [p(2); slowingAccel(p(2))];
%     [t, y] = ode45(odefun, tspan, y0);
% 
%     % % --- Plot z vs t ---
%     % figure(2);
%     % h1 = plot(t, y(:,1), 'DisplayName', ['disp = ' currDisp ' mm'], 'LineWidth', 1.5);
%     % hold on;
%     % z_target = slowingLength;
%     % crossingIdx = find(y(:,1) >= z_target, 1, 'first');
%     % if ~isempty(crossingIdx)
%     %     t_cross = t(crossingIdx);
%     %     yline(z_target, 'k--', 'LineWidth', 1.2);
%     %     xline(t_cross, 'r--', 'LineWidth', 1.2);
%     %     plot(t_cross, z_target, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
%     % end
%     % xlabel('t (ms)');
%     % ylabel('z (mm)');
%     % title('z vs t Trajectory');
%     % legend('show', 'Location', 'northeastoutside');
% 
%     % % --- Plot v_z vs t ---
%     % figure(3);
%     % h2 = plot(t, y(:,2), 'DisplayName', ['disp = ' currDisp ' mm'], 'LineWidth', 1.5);
%     % hold on;
%     % zeroCrossingIdx = find(y(:,2) <= 0, 1, 'first');
%     % if ~isempty(zeroCrossingIdx)
%     %     t_cross2 = t(zeroCrossingIdx);
%     %     yline(0, 'k--', 'LineWidth', 1.2);
%     %     xline(t_cross2, 'r--', 'LineWidth', 1.2);
%     %     plot(t_cross2, 0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
%     % end
%     % xlabel('t (ms)');
%     % ylabel('v_z (m/s)');
%     % title('v_z vs t Trajectory');
%     % legend('show', 'Location', 'northeastoutside');
% 
%     % --- Plot v_z vs z ---
%     figure(4);
%     h3 = plot(y(:,1), y(:,2), 'DisplayName', ['disp = ' currDisp ' mm'], 'LineWidth', 1.5);
%     hold on;
%     z_target = slowingLength;
%     crossingIdx = find(y(:,1) >= z_target, 1, 'first');
%     if ~isempty(crossingIdx)
%         v_at_z_target = y(crossingIdx,2);
%         yline(v_at_z_target, 'k--', 'LineWidth', 1.2);
%         xline(z_target, 'r--', 'LineWidth', 1.2);
%         plot(z_target, v_at_z_target, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
%     end
%     xlabel('z (mm)');
%     ylabel('v_z (m/s)');
% 
%     title('v_z vs z Trajectory');
%     legend('show', 'Location', 'northeastoutside');
% 
% end


%---------------------------------------------------------
% compute a_z(v_z, r) [i.e. account for finite beam waist)
%---------------------------------------------------------

%pos = {'0.01', '0.1', '0.5', '1.0', '2.0', '3.0', '5.0', '7.0', '9.0', '10.0'};
pos = {'0.01', '0.1', '0.5', '1.0', '2.0', '3.0', '5.0', '7.0', '9.0', '10.0', '12.0', '14.0', '17.0'};
dispList = str2double(pos); % convert to numeric displacement values
vGrid = linspace(0, 300, 200); % common velocity grid (adjust as needed)
azInterpMatrix = zeros(length(pos), length(vGrid)); % rows: displacements, cols: velocities

for i = 1:length(pos)
    currFile = strcat(dataFolder, '/forcevsSpeedDisplacement', pos{i}, 'MMRandom.dat');
    currData = readtable(currFile);

    % Sort and interpolate
    [LongSpeedsSorted, sortIdx] = sort(currData.LongSpeed);
    azsSorted = currData.az(sortIdx) / 1e3; % Convert to mm/ms^2

    % Interpolate to common vGrid
    azInterpMatrix(i, :) = interp1(LongSpeedsSorted, azsSorted, vGrid, 'spline', 'extrap');
end

[VV, DD] = meshgrid(vGrid, dispList); % meshgrid: cols = vGrid, rows = dispList

% figure;
% surf(VV, DD, azInterpMatrix, 'EdgeColor', 'none');
% xlabel('v_z (m/s)');
% ylabel('Transverse Displacement (mm)');
% zlabel('a_z (mm/ms^2)');
% title('Acceleration vs Velocity vs Transverse Displacement');
% colorbar;
% view(135, 30); % 3D view angle
% 
% figure('Color', 'w');
% imagesc(vGrid, dispList, azInterpMatrix); % rows = y-axis, cols = x-axis
% axis xy; % Ensure y-axis increases upward (displacement direction)
% xlabel('v_z (m/s)', 'FontSize', 12);
% ylabel('Displacement (mm)', 'FontSize', 12);
% title('a_z(v_z, r)', 'FontSize', 14);
% colorbar;
% colormap parula; % or try 'parula', 'hot', 'magma', 'viridis' etc.
% set(gca, 'FontSize', 12);

azInterpFcn = @(vz, disp) interp2(vGrid, dispList, azInterpMatrix, vz, disp, 'spline');

% a = azInterpFcn(150, 3.5); % Evaluate at vz = 150 m/s, disp = 3.5 mm

%----------------------------------------------------
% Monte Carlo trajectory simulation 
% incorporate longitudinal deceleration force
% incorporate initial transverse displacement in 2-D
% incorporate initial transverse velocity in 2-D
% incorporate finite beam waist
% parallelized ode45 solver (native Matlab method)
%----------------------------------------------------

rng(12345);

% Parameters
numTrials = 2500000;
z_target = 500; % mm
v_capture = 27; % m/s
tspan = [0, 30]; % ms
initialZ = 0; % mm
vzMean = 140; vzSigma = 32; % mm/ms
rInitMax = 1.5; % mm (max transverse displacement)
vrSigma = 32.0; % m/s (std dev transverse velocity)
capture_r = 7.0; % mm

% % Random seed (optional for reproducibility)
% rng(1);

% Generate initial conditions outside loop
vz0 = normrnd(vzMean, vzSigma, numTrials, 1);  % mm/ms

thetaPos = 2 * pi * rand(numTrials, 1);
radii = rInitMax * sqrt(rand(numTrials, 1));
x0 = radii .* cos(thetaPos);  % mm
y0 = radii .* sin(thetaPos);  % mm

% initial transverse velocity distribution (indep. x and y)
vx0 = normrnd(0, vrSigma, numTrials, 1);  % m/s
vy0 = normrnd(0, vrSigma, numTrials, 1);  % m/s

% Preallocate outputs
zFinal = zeros(numTrials, 1);
vFinal = zeros(numTrials, 1);
reachedTarget = false(numTrials, 1);
tAtCapture = NaN(numTrials, 1);

% Preallocate cell array for successful trajectories (still parallel safe)
trajZ = cell(numTrials, 1);  % z(t)
trajVz = cell(numTrials, 1); % v_z(t)
% trajX = cell(numTrials, 1);
% trajY = cell(numTrials, 1);

fprintf('Running %d Monte Carlo trials with parfor...\n', numTrials);

% Enable parallel pool (if not already open)
if isempty(gcp('nocreate'))
    parpool;
end

tic;

% Main simulation loop (parallelized)
parfor i = 1:numTrials
    % Compute initial radial position
    r0_mag = sqrt(x0(i)^2 + y0(i)^2);

    % Slowing acceleration function with fixed r
    azInterpLocal = @(vz) interp2(vGrid, dispList, azInterpMatrix, vz, r0_mag, 'spline');
    odefun = @(t, p) [p(2); azInterpLocal(p(2))];
    y0i = [initialZ; vz0(i)];

    [t, y] = ode45(odefun, tspan, y0i);

    zFinal(i) = y(end, 1);
    vFinal(i) = y(end, 2);

    crossingIdx = find(y(:,1) >= z_target, 1, 'first');
    if ~isempty(crossingIdx)
        t_cross = t(crossingIdx);
        vzAtTarget = y(crossingIdx, 2);

        x_target = x0(i) + vx0(i) * t_cross;
        y_target = y0(i) + vy0(i) * t_cross;
        r_target = sqrt(x_target^2 + y_target^2);

        trajZ{i} = y(:,1);
        trajVz{i} = y(:,2);

        if vzAtTarget <= v_capture && r_target <= capture_r
            reachedTarget(i) = true;
            tAtCapture(i) = t_cross;

            % % Store successful trajectory (z, v_z)
            % trajZ{i} = y(:,1);
            % trajVz{i} = y(:,2);
        end
    end
end

elapsedTime = toc;
fprintf('Simulation completed in %.2f seconds.\n', elapsedTime);

% Report results
numReached = sum(reachedTarget);
fprintf('Captured at z = %.1f mm (v_z <= %.1f mm/ms and |r| < 7 mm): %d / %d (%.2e)\n', ...
    z_target, v_capture, numReached, numTrials, numReached/numTrials);

averageCaptureTime = mean(tAtCapture(reachedTarget), 'omitnan');
fprintf('Average time of capture = %.3f ms\n', averageCaptureTime);

% Plot only successful trajectories
figure;
hold on;
for i = 1:numTrials
    if reachedTarget(i)
        plot(trajZ{i}, trajVz{i}, 'Color', [0, 0.6, 0], 'LineWidth', 1.0); % green for success
    end
end
xline(z_target, 'k--', 'LineWidth', 1.2);
yline(v_capture, 'k--', 'LineWidth', 1.2);
xlabel('z (mm)', 'FontSize', 12);
ylabel('v_z (mm/ms)', 'FontSize', 12);
title('Successful Monte Carlo Trajectories: v_z vs z', 'FontSize', 14);
grid on;

% % Combine and export all trajectories with r = sqrt(x^2 + y^2)
% % rowsSuccess = {};
% % rowsFail = {};
% rowsSuccessLong = {};
% rowsFailLong = {};
% 
% for i = 1:numTrials
%     % trajData = table();        % full trajectory (z, vz, x, y, r)
%     trajDataLong = table();    % only z and vz
% 
%     if ~isempty(trajZ{i})
%         nPts = numel(trajZ{i});
%         trajID = repmat(i, nPts, 1);
% 
%         % % Compute r = sqrt(x^2 + y^2)
%         % rVals = sqrt(trajX{i}.^2 + trajY{i}.^2);
%         % 
%         % % Full trajectory
%         % trajData = table( ...
%         %     trajID, ...
%         %     (1:nPts)', ...
%         %     trajZ{i}, ...
%         %     trajVz{i}, ...
%         %     trajX{i}, ...
%         %     trajY{i}, ...
%         %     rVals, ...
%         %     'VariableNames', {'trajID', 'step', 'z', 'vz', 'x', 'y', 'r'} ...
%         % );
% 
%         % Longitudinal only
%         trajDataLong = table( ...
%             trajID, ...
%             (1:nPts)', ...
%             trajZ{i}, ...
%             trajVz{i}, ...
%             'VariableNames', {'trajID', 'step', 'z', 'vz'} ...
%         );
%     end
% 
%     if reachedTarget(i)
%         if ~isempty(trajData)
%             % rowsSuccess{end+1} = trajData;
%             rowsSuccessLong{end+1} = trajDataLong;
%         end
%     else
%         % Export only every 1000th unsuccessful trajectory
%         if ~isempty(trajData) && mod(i, 1000) == 0
%             % rowsFail{end+1} = trajData;
%             rowsFailLong{end+1} = trajDataLong;
%         end
%     end
% end
% 
% % Timestamp to avoid overwriting
% timestamp = datestr(now, 'yyyymmdd_HHMM');
% % success_filename       = ['successful_trajectories_with_r_' timestamp '.csv'];
% % fail_filename          = ['unsuccessful_trajectories_with_r_' timestamp '.csv'];
% success_long_filename  = ['successful_trajectories_longitudinal_' timestamp '.csv'];
% fail_long_filename     = ['unsuccessful_trajectories_longitudinal_' timestamp '.csv'];
% 
% % Write to CSV
% if ~isempty(rowsSuccess)
%     % T_success = vertcat(rowsSuccess{:});
%     % writetable(T_success, success_filename);
%     % fprintf('Saved %s\n', success_filename);
% 
%     T_success_long = vertcat(rowsSuccessLong{:});
%     writetable(T_success_long, success_long_filename);
%     fprintf('Saved %s\n', success_long_filename);
% end
% 
% if ~isempty(rowsFail)
%     % T_fail = vertcat(rowsFail{:});
%     % writetable(T_fail, fail_filename);
%     % fprintf('Saved %s\n', fail_filename);
% 
%     T_fail_long = vertcat(rowsFailLong{:});
%     writetable(T_fail_long, fail_long_filename);
%     fprintf('Saved %s\n', fail_long_filename);
% end



% % Optional: plot a summary figure
% figure;
% histogram(tAtCapture(reachedTarget), 50);
% xlabel('Capture Time (ms)');
% ylabel('Count');
% title('Histogram of Capture Times');
% grid on;





% % Parameters
% numTrials = 100000;
% z_target = 500; % mm
% v_capture = 27; % m/s
% tspan = [0, 30]; % ms
% initialZ = 0; % mm
% vzMean = 140; vzSigma = 32; % mean initial longitudinal velocity and std dev
% rInitMax = 1.5; % mm (uniform disk radius for transverse displacement)
% vrSigma = 32.0; % m/s (std dev of 2D transverse velocity)
% capture_r = 7.0; % mm (MOT capture radius)
% 
% % Generate initial conditions
% vz0 = normrnd(vzMean, vzSigma, numTrials, 1);  % initial v_z (mm/ms)
% 
% % --- 2D initial transverse displacement (uniform disk) ---
% thetaPos = 2 * pi * rand(numTrials, 1);
% radii = rInitMax * sqrt(rand(numTrials, 1)); % need sqrt to ensure uniform distribution in area
% x0 = radii .* cos(thetaPos);  % mm
% y0 = radii .* sin(thetaPos);  % mm
% 
% % --- 2D transverse velocity (Gaussian) ---
% vx0 = normrnd(0, vrSigma, numTrials, 1);  % m/s
% vy0 = normrnd(0, vrSigma, numTrials, 1);  % m/s
% 
% % Interpolator: a_z(v_z, r)
% azInterpFcn = @(vz, r) interp2(vGrid, dispList, azInterpMatrix, vz, r, 'spline');
% 
% % Preallocate results
% zFinal = zeros(numTrials, 1);
% vFinal = zeros(numTrials, 1);
% rAtTarget = zeros(numTrials, 1);
% reachedTarget = false(numTrials, 1);
% tAtCapture = NaN(numTrials, 1); % Time of capture
% 
% % Run Monte Carlo simulations
% fprintf('Running %d Monte Carlo trials...\n', numTrials);
% tic;
% 
% % figure;
% % hold on;
% 
% for i = 1:numTrials
%     % Compute initial radial distance
%     r0_mag = sqrt(x0(i)^2 + y0(i)^2);  % mm
% 
%     % Slowing function with fixed r for this atom
%     slowingAccel = @(v) azInterpFcn(v, r0_mag);
%     odefun = @(t, p) [p(2); slowingAccel(p(2))];  % p = [z; v_z]
%     inits = [initialZ; vz0(i)];
% 
%     [t, y] = ode45(odefun, tspan, inits);
% 
%     zFinal(i) = y(end, 1);
%     vFinal(i) = y(end, 2);
% 
%     % Check for capture at z_target
%     crossingIdx = find(y(:,1) >= z_target, 1, 'first');
%     if ~isempty(crossingIdx)
%         t_cross = t(crossingIdx);  % ms
%         vzAtTarget = y(crossingIdx, 2);  % mm/ms
% 
%         % Compute transverse displacement at target
%         x_target = x0(i) + vx0(i) * t_cross;  % mm
%         y_target = y0(i) + vy0(i) * t_cross;  % mm
%         r_target = sqrt(x_target^2 + y_target^2);  % mm
% 
%         if vzAtTarget <= v_capture && r_target <= capture_r
%             reachedTarget(i) = true;
%             tAtCapture(i) = t_cross;
%         end
%     end
% 
%     % Optional plotting (disable for speed)
%     %color = reachedTarget(i) * [0 0.6 0] + ~reachedTarget(i) * [0.6 0 0]; % green/red
%     %plot(y(:,1), y(:,2), 'Color', color, 'LineWidth', 1.2); % z vs vz
% end
% 
% elapsedTime = toc;
% fprintf('Monte Carlo simulation completed in %.2f seconds.\n', elapsedTime);
% 
% % Report results
% numReached = sum(reachedTarget);
% fprintf('Captured at z = %.1f mm (v_z <= %.1f mm/ms): %d / %d (%.2e)\n', ...
%     z_target, v_capture, numReached, numTrials, numReached/numTrials);
% 
% averageCaptureTime = mean(tAtCapture(reachedTarget), 'omitnan');
% fprintf('Average time of capture = %.3f ms\n', averageCaptureTime);
% 
% % % Plot reference lines
% % xline(z_target, 'k--', 'LineWidth', 1.2);
% % yline(v_capture, 'k--', 'LineWidth', 1.2);
% % 
% % xlabel('z (mm)', 'FontSize', 12);
% % ylabel('v_z (mm/ms)', 'FontSize', 12);
% % title('Monte Carlo Trajectories: v_z vs z', 'FontSize', 14);
% % grid on;
% % legend('Trajectory', 'Location', 'northeastoutside');