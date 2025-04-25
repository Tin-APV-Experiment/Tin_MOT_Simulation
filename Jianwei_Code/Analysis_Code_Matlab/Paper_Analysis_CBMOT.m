clear all;
colors = [[0.5;0.2;0.8],[0.4;0.7;0.8],[0.8;0.5;0.2],[0.5;0.7;0.2],...
    [0.6;0.6;0.3],[0.6;0.3;0.6],[0.3;0.6;0.6],[0.4;0.4;0.8],[0.4;0.8;0.4],[0.8;0.4;0.4],[0.5;0.5;0.5],[0.9;0.2;0.1],...
    [0.2;0.4;0.7],[0.4;0.4;0.0],[0.0;0.4;0.2],[0.2;0.0;0.4],[0.0;0.0;0.0],[0.8;0.8;0.8]];

simTrapDynamics =1; 
%change to 1 if you want to simulate particle trajectory, with random photon scatter, to get 'true' size and temperature
SavePlot=1; %change to 1 if you want to save the Plot

fa=1;
fb=1;

moleculeName = "Sn";
molecule = def_molecule(moleculeName);
bfield_grad = "60"; %b-field gradient used
data_timestamp = "20250412_0121";
laser_waist = 7; %beam waist radius in mm, which roughly defines effective size of trap
useRandPhase = 1; %flag set in Julia code if random phases were used for the laser

Trails= "20";
Date= "20250412";
MOTCondition = strcat("Bgrad=", bfield_grad, "Trail=", Trails); % Guass/cm %'_IR=1.0_IB=0.1_DR=0.5_DB=0.1'
LaserCondition = " ";

%here is for Sn atoms
base_folder = strcat("SnData/",Date,"/");
dataFolder = strcat(base_folder, moleculeName, "PapertestBGradGPerCM", bfield_grad, "Date", data_timestamp,"Trail",Trails);
if useRandPhase == 1
    dataFolder = strcat(dataFolder, "_WithRandPhase");
end
%default position settings in MOT simulation
%pos = {'1', '2', '3', '5', '7', '9', '11', '14', '17'};
%pos={'0.5','1.0','1.5','2.0','2.5','3.0','5.0'};
%pos = {'1', '2', '3', '5', '7', '9', '11'};
%compressed MOT pos settings
%pos = {'0.5', '1.0', '1.5', '2.0','2.5', '3.0', '5.0'};
%blueMOT pos settings
%pos = {'0.1', '0.3', '0.5', '0.7', '0.9', '1.5','2.0'};
%pos={'0.01','0.03','0.05','0.07','0.1','0.3','0.5','0.7','1.0','1.5','3.0'}; 
pos={'0.01','0.03','0.05','0.1','0.3','0.6','0.9','1.2','1.5','2.0','3.0'};
%pos={'0.01','0.03','0.05','0.1','0.3','0.6','1.0','1.5','3.0'}; 
%pos = {'0.2', '0.5', '1.0', '2.0', '3.0', '4.0', '5.0'};
%pos = {'0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6'};
%pos = {'0.01','0.03','0.08', '0.15', '0.3', '0.5', '0.9', '1.5','3.0','5.0'};
%{
%other default setting pos
pos = {'0.5', '1.5', '3.0', '4.5', '6.0', '7.5', '9.0', '10.5', '12.0', '13.5', '15.0', '17.0'};
%}
%displacement={'0.2','0.5','1.0','1.5','2.5','3.5','5.0'};
%pos = {'0.5', '1.0', '1.5', '2.0', '3.0', '4.0', '5.0', '6.0', '7.5', '9.0', '10.5', '12.0', '14.0', '16.0'};

for i=1:length(pos)
    posForPlot(i) = str2num(pos{i});
end
posInMM = posForPlot.*1e-3;
testFile = strcat(dataFolder,'/forcevsSpeedDisplacement',pos{1},'MMSameDir.dat');
testData = readtable(testFile);

accelsInVelDirection = zeros(size(testData,1),length(pos));
accelsInPosDirection = zeros(size(testData,1),length(pos));
excitedPop = zeros(size(testData,1),length(pos));
for i=1:length(pos)
    currFile = strcat(dataFolder,'/forcevsSpeedDisplacement',pos{i},'MMSameDir.dat');
    currData = readtable(currFile);
    vels = currData.Speed;
    accelsInVelDirection(:,i) = currData.av;
    accelsInPosDirection(:,i) = currData.ar;
    excitedPop(:,i) = currData.PFeHigh;
end
%here we inherently assume the symmetric nature of MOT with respect to
%displacement and velocity

%reverse sign of accel for negative velocities since it's called
%acceleration in velocity direction
for i=1:length(vels)
    if vels(i)<0
        accelsInVelDirection(i,:) = accelsInVelDirection(i,:).*-1;
    end
end

%add opposite positions
accelsInVelDirectionFull = [-flipud(fliplr(accelsInVelDirection)),accelsInVelDirection];
excitedPopFull = [flipud(fliplr(excitedPop)),excitedPop];
sortedPos = sort([-posForPlot,posForPlot]);

%simulate capture with linear interpolation (spline acts up here for some
%reason.  Probably not a big deal to use linear)

%oneAxisAccel is function defined to interpolate accelerations for any
%(d,v) within the grid defined by our choice of d's and v's from simulation
%here, the interpolation is linear (not spline) so our capture velocity
%calculation does not mess up (otherwise it doesn't run in time)

%diffEqVals is a function that returns a vector (v, a) [here p(1) is d and
%p(2) is v]
oneAxisAccel = @(d,v) interp2(sortedPos,vels,accelsInVelDirectionFull./1e3,d,v,'linear');%in mm,ms units
diffEqVals = @(t,p) [p(2);oneAxisAccel(p(1),p(2))];
vsToTryForCapture = [1:0.1:50]; %step size is 0.1 m/s
%cap_threshold = sqrt(2)*laser_waist; %set at sqrt(2)*w but can modify if necessary
cap_threshold = laser_waist;
%use ode45 (4th/5th order Runge-Kutta method) to solve for trajectory (out
%to 50 ms). 
%Also for every currV that is a multiple of 3, plot its
%trajectory (i.e. plot multiple phase space trajectories on same plot for
%the speeds below capture velocity).


%search for capture velocity and plot lower velocity phase plots (Fig. 1)
%{
for i=1:length(vsToTryForCapture)
    currV = vsToTryForCapture(i);
    [ts2,ps2] = ode45(diffEqVals,[0 50],[min(sortedPos);currV]);
    if any(ps2(:,1) > cap_threshold) || isnan(ps2(2,1))
        break;
    end
end
capVel = currV-0.3;
%}

%{
% Plot multiple phase space trajectories on same plot
figure;
hold on; % This ensures all plots are on the same figure

% Initialize an array to store plot handles
plotHandles = gobjects(currV, 1); % Preallocate for better performance

% Initialize a cell array to store legend labels
legendLabels = cell(currV, 1);

% Loop from currV down to 1, in steps of -3
for v = currV:-3:1
    % Solve the ODE with the current velocity
    [ts2, ps2] = ode23(diffEqVals, [0 50], [min(sortedPos); v]);
    
    % Plot the entire trajectory and store the plot handle
    plotHandles(v) = plot(ps2(:, 1), ps2(:, 2), 'LineWidth',2);
    
    % Store the legend label for the current plot
    legendLabels{v} = ['v_0 = ' num2str(v) ' m/s'];
end

% Set labels and title
xlabel('d (mm)');
ylabel('v_{d} (m/s)');
title('Particle Trajectories for Different Initial Velocities');

% Create the legend
legend(plotHandles(currV:-3:1), legendLabels(currV:-3:1), 'Location', 'best');
%saveas(gcf, 'phase_space_MOT_capture_highDPI.png');
%print(gcf, 'phase_space_MOT_capture_highDPI', '-dpng', '-r300'); % 300 DPI

hold off;
%}
%{
print('HighDPIPlot', '-dpng', '-r300'); % Saves as 'HighDPIPlot.png'


%}


%plot the highest value of v for which we have capture (evolve out to t=100 ms)
%[ts2,ps2] = ode23(diffEqVals,[0 100],[min(sortedPos);currV-0.2]);
%[ts2,ps2] = ode23(diffEqVals,[0 100],[-9;currV-0.2]); %start from z=-9 mm

%{

%Plot RedCapture Trajectory

[ts2,ps2] = ode45(diffEqVals,[0 50],[min(sortedPos); capVel]);
figure(1);
hold on;
plot(ps2(:,1),ps2(:,2),'LineWidth',2,'DisplayName', strcat(MOTCondition,'  v_{Cap}=', num2str(currV-0.2), ' m/s',LaserCondition));
xlabel('d (mm)', fontsize=12);
ylabel('v_{d} (m/s)', fontsize=12);
legend('Location', 'best');
title(strcat('Sn Atom Trajectory for Capture Velocity'), fontsize=14)

output_folder_capvel = strcat('SnBluePlot','/CaptureVelucity/','/',Date,'/');
if ~exist(output_folder_capvel, 'dir')
    mkdir(output_folder_capvel);
end
%title(strcat('Sn Atom Trajectory for Capture Velocity =',num2str(currV-0.2),' m/s'), fontsize=14)
%saveas(gcf, strcat('SnPlot','/CaptureVelucity/','BGradCompare.png'));
print(gcf, strcat(output_folder_capvel,data_timestamp,MOTCondition), '-dpng', '-r300'); % 300 DPI
%print(gcf, strcat('SnPlot','/CaptureVelucity/','BGradCompareHDPI'), '-depsc','-r300'); % 300 DPI
%print(gcf, strcat('SnPlot','/CaptureVelucity/','BGradCompareHDPI'), '-dpng', '-r300')); % 300 DPI

output_folder = 'SnMOTTrials';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

writematrix(ps2, fullfile(output_folder, strcat('CaptureVelocityTrajectory', data_timestamp, '.csv')));

%write to csv
writematrix(ps2, strcat('SnMOTTrials/CaptureVelocityTrajectory', data_timestamp, '.csv'));

%}

%{
%plot lower velocity trajectories as well, and for 50 ms instead of 20 ms
[ts2,ps2] = ode23(diffEqVals,[0 50],[min(sortedPos);currV-5]);
figure(1);

% Filter the time points between t = 4 and t = 50
timeIndices = (ts2 >= 4) & (ts2 <= 50);

% Plot the filtered data
plot(ps2(timeIndices, 1), ps2(timeIndices, 2));
xlabel('z(mm)');
ylabel('v (m/s)')
title('Particle Trajectory from t = 4 to t = 50');

%plot trajectory for a speed lower than capture velocity
%{
plot(ps2(:,1),ps2(:,2));
xlabel('x(mm)');
ylabel('v (m/s)')
title(strcat('particle trajectory for v_{Cap}=',num2str(currV-8),' m/s'))
%}
%}


%now plot a(x,v) heat map
oneAxisAccel = @(d,v) interp2(sortedPos,vels,accelsInVelDirectionFull./1e3,d,v,'spline');%in mm,ms units


%make heat map
vsForHeatMap = [min(vels):.01:max(vels)].*fa;
xsForHeatMap = [min(sortedPos):.01:max(sortedPos)].*fb;
for i=1:length(vsForHeatMap)
    for j=1:length(xsForHeatMap)
        heatMap(i,j) = oneAxisAccel(xsForHeatMap(j),vsForHeatMap(i));
    end
end
figure(1)
hold on;
imagesc(xsForHeatMap,vsForHeatMap,heatMap);

% Apply red-white-blue colormap
n = 256;
rwb = zeros(n,3);
caxis([-10 10]);
for i = 1:n
    t = (i-1)/(n-1);
    if t < 0.5
         rwb(i,:) = [t*2, t*2, 1];              % blue to white
        
    else
        rwb(i,:) = [1, 2 - 2*t, 2 - 2*t];  % white to red
       
    end
end
colormap(flipud(rwb));

colorbar
xlabel('d (mm)');
ylabel('v_{d} (m/s)')
grid on;
grid minor; 
title('Sn Convelt-Belt MOT Dynamics Heat Map')
xlim([min(xsForHeatMap) max(xsForHeatMap)])
h=colorbar;
h.Title.String = "a_{d} (mm/ms^{2})";

if SavePlot==1
    %make sure the folder exist
    output_folder_Haetplot = strcat('SnCBMOTPlot','/HeatPlot/','/',Date,'/');
    if ~exist(output_folder_Haetplot, 'dir')
        mkdir(output_folder_Haetplot);
    end
    %saveas(gcf, strcat('SnPlot','/HeatPlot/',data_timestamp,MOTCondition,'.png'))
    print(gcf, strcat(output_folder_Haetplot,data_timestamp,MOTCondition,'HDPI'),'-dpng','-r300');
end


%make x,v plots from matrix
%vMaxToInt = 4;
%xMaxToInt = 7;
% Geoffrey's code
%for capture MOT
vMaxToInt = 2;
xMaxToInt = 2;
[~,minCol] = min(abs(xsForHeatMap+xMaxToInt));
[~,maxCol] = min(abs(xsForHeatMap-xMaxToInt));
[~,minRow] = min(abs(vsForHeatMap+vMaxToInt));
[~,maxRow] = min(abs(vsForHeatMap-vMaxToInt));
accVsPosForPlot = trapz(vsForHeatMap(minRow:maxRow),heatMap(minRow:maxRow,:),1)./(2*vMaxToInt);
accVsVelForPlot = trapz(xsForHeatMap(minCol:maxCol),heatMap(:,minCol:maxCol),2)./(2*xMaxToInt);

figure(2);
hold all;
plot(xsForHeatMap,accVsPosForPlot,'Linewidth',2,'DisplayName', strcat(MOTCondition,LaserCondition));
xlabel('d (mm)');
ylabel('a_{d} (mm/ms^{2})')
legend('Location', 'best');
%xlim([0 15])
%ylim([-1 0.5])
grid on;
grid minor;
title('Acceleration vs Position, Convelt-Belt MOT')

if SavePlot==1
    output_folder_pos = strcat('SnCBMOTPlot','/AccelerationVsPosition/','/',Date,'/');
    if ~exist(output_folder_pos, 'dir')
        mkdir(output_folder_pos);
    end
    %save figure(3)
    %saveas(gcf, strcat('SnPlot','/AccelerationVsPosition/','BGradcompare.png')
    print(gcf, strcat(output_folder_pos,data_timestamp,MOTCondition),'-dpng','-r300');
end

figure(3);
hold all;
plot(vsForHeatMap,accVsVelForPlot,'LineWidth',2,'DisplayName', strcat(MOTCondition,LaserCondition));
legend('Location', 'best');
xlabel('v_{d} (m/s)');
ylabel('a_{d} (mm/ms^{2})')
grid on;
grid minor; 
%xlim([0 20])
%ylim([-5 1])
title('Acceleration vs Velocity, Convelt-Belt MOT')

if SavePlot==1
    output_folder_vel = strcat('SnCBMOTPlot','/AccelerationVsVelocity/','/',Date,'/');
    if ~exist(output_folder_vel, 'dir')
        mkdir(output_folder_vel);
    end
    print(gcf, strcat(output_folder_vel,data_timestamp,MOTCondition),'-dpng','-r300');
end

%write to csv
HeatMap_M = [NaN, xsForHeatMap(:)'; vsForHeatMap(:),heatMap];

writematrix(HeatMap_M, strcat('SnMOTPaper/Heatmap', data_timestamp, '.csv'));
writematrix([xsForHeatMap; accVsPosForPlot]', strcat('SnMOTPaper/AccelVsPos', data_timestamp, '.csv'));
writematrix([vsForHeatMap; accVsVelForPlot']', strcat('SnMOTPaper/AccelVsVel', data_timestamp, '.csv'));

%account for spontaneous scattering of photons from captured atoms in MOT
%to computer average temperature and size
if simTrapDynamics == 1
    %for the region in (d, v) space defined by +/-xMaxToInt and +/-vMaxToInt,
    %compute the average excited state population there (i.e. average excited
    %state population in MOT capture region).
    [~,minCol] = min(abs(sortedPos+xMaxToInt));
    [~,maxCol] = min(abs(sortedPos-xMaxToInt));
    [~,minRow] = min(abs(vels+vMaxToInt));
    [~,maxRow] = min(abs(vels-vMaxToInt));
    meanExcPop = mean(mean(excitedPopFull(minRow:maxRow,minCol:maxCol)));

    maxTime = 100; %100 ms time evolution after capture in MOT
    initPos = -0.0; %assume start at center of t rap
    initVel = 0; %assume start with zero velocity
    gam = molecule.gam;
    kA = molecule.kA;
    mass = molecule.mass;
    hbar = 1.05e-34;
    kb = 1.38e-23;
    %plot MOT trajectory including random photon scatter
    figure(5);
    hold on;
    % Solve the differential equation and plot the trajectory (without photon scatter)
    %[ts2,ps2] = ode23(diffEqVals,[0 100],[1;1]); %1;1 indicates initial v and r
    [ts2,ps2] = ode45(diffEqVals,[0 maxTime],[initPos;initVel]);
    plot(ps2(:,1), ps2(:,2), 'LineWidth', 2, 'DisplayName', 'Without Photon Scatter');

    % Use the trap dynamics algorithm to plot the trajectory with photon scatter
    scatterRate = meanExcPop * gam * 1e-3; % in 1/ms
    tKick = 1 / scatterRate;
    r(1) = initPos; % mm
    v(1) = initVel; % mm/ms
    vKick = hbar * kA / mass;

    for i = 1:round(maxTime / tKick)
        randPhi1 = 2 * pi * rand;
        randPhi2 = 2 * pi * rand;
        v(i + 1) = v(i) + vKick * (cos(randPhi1) + cos(randPhi2)) + oneAxisAccel(r(i), v(i)) * tKick;
        r(i + 1) = r(i) + v(i) * tKick;
    end
    
    % Plot the trajectory with photon scatter
    plot(r, v, '-', 'LineWidth', 2, 'DisplayName', 'With Photon Scatter'); % Thicker line and legend entry
    
    % Add grid, labels, and title
    grid on;
    xlabel('d (mm)');
    ylabel('v_{d} (mm/ms)');
    title('Particle Trajectory for Sn atoms ');
    
    % Create a legend
    legend('Location', 'best');
    
    simTimes = 0:(tKick):maxTime;
    startTime = maxTime/2;
    startInd = maxTime/2/tKick;
    endInd= i;
    vSq=mean(v(startInd:endInd).^2);
    subarraypos=r(startInd:endInd);
    motpos=mean(subarraypos);
    sigma=sqrt(mean((subarraypos-motpos).^2)); %in mm
    temp = vSq*mass/kb*1e3; %in mK
    
    %write to csv
    if SavePlot==1
        output_folder_traj = strcat('SnCBMOTPlot','/particle_trajectory_MOT/','/',Date,'/');
        if ~exist(output_folder_traj, 'dir')
            mkdir(output_folder_traj);
        end
        % Save the figure as a 300 DPI image
        %saveas(gcf, strcat('SnPlot','/particle_trajectory_MOT/',data_timestamp,MOTCondition,'.png'));
        print(gcf, strcat(output_folder_traj,data_timestamp,MOTCondition,'HDPI'),'-dpng', '-r300');% Save as HDPI Png
    end
    

end

