% More Particle MOT Dynamics
clear all;
SavePlot =1;


colors = [[0.5;0.2;0.8],[0.4;0.7;0.8],[0.8;0.5;0.2],[0.5;0.7;0.2],...
     [0.6;0.6;0.3],[0.6;0.3;0.6],[0.3;0.6;0.6],[0.4;0.4;0.8],[0.4;0.8;0.4],[0.8;0.4;0.4],[0.5;0.5;0.5],[0.9;0.2;0.1],...
     [0.2;0.4;0.7],[0.4;0.4;0.0],[0.0;0.4;0.2],[0.2;0.0;0.4],[0.0;0.0;0.0],[0.8;0.8;0.8]];

simTrapDynamics = 1; 
CMOT = 0; %1 -> 1mK ; 2 ->8.3mK; 0 -> initial velocity and postition will both be 0;3->Other

if CMOT ==3
    SavePlot =0;
end

initPosset=2.4;
numParticles = 100; % Number of particles to simulate

moleculeName = "Sn";
molecule = def_molecule(moleculeName);
bfield_grad = "20"; %b-field gradient used
data_timestamp = "20250421_1100";
laser_waist = 7; %beam waist radius in mm, which roughly defines effective size of trap
useRandPhase = 1; %flag set in Julia code if random phases were used for the laser

Trails= "150";
Date= "20250421";
MOTCondition = strcat("Bgrad=", bfield_grad, "Trail=", Trails); 
LaserCondition = " ";


% Data folder paths
base_folder = strcat("SnData/",Date,"/");
dataFolder = strcat(base_folder, moleculeName,"PapertestBGradGPerCM", bfield_grad, "Date", data_timestamp,"Trail",Trails);
if useRandPhase == 1
    dataFolder = strcat(dataFolder, "_WithRandPhase");
end

% Position settings for MOT simulation

%pos = {'0.1', '0.3', '0.5', '0.7', '0.9', '1.5','2.0'};
%pos = {'0.05', '0.15', '0.25', '0.45', '0.65', '0.85','1.0','1.5'};
%pos = {'0.01', '0.02', '0.03', '0.05', '0.075', '0.1','0.3','0.6','1.0','1.5'};
%pos={'0.01','0.03','0.05','0.07','0.1','0.3','0.5','0.7','1.0','1.5','3.0'}; 

%pos={'0.01','0.03','0.05','0.1','0.3','0.6','0.9','1.2','1.5','2.0','3.0'};
%pos={'0.01','0.03','0.05','0.1','0.3','0.6','1.0','1.5','3.0'}; 
pos={'0.01','0.03','0.05','0.1','0.3','0.6','0.9','1.2','1.5','2.0','2.5','3.0','4.0'};

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

% Reverse sign of acceleration for negative velocities
for i=1:length(vels)
    if vels(i)<0
        accelsInVelDirection(i,:) = accelsInVelDirection(i,:).*-1;
    end
end

accelsInVelDirectionFull = [-flipud(fliplr(accelsInVelDirection)),accelsInVelDirection];
excitedPopFull = [flipud(fliplr(excitedPop)),excitedPop];
sortedPos = sort([-posForPlot,posForPlot]);

oneAxisAccel = @(d,v) interp2(sortedPos,vels,accelsInVelDirectionFull./1e3,d,v,'spline');%in mm,ms units
oneAxisExppop = @(d,v) interp2(sortedPos,vels,excitedPopFull,d,v,'spline');%in mm,ms units

if simTrapDynamics == 1
    maxTime = 50;
    gam = molecule.gam;
    kA = molecule.kA;
    mass = molecule.mass;
    hbar = 1.05e-34;
    kb = 1.38e-23;
    
    
    scatterRate = mean(mean(excitedPopFull)) * gam * 1e-3;
    tKick = 1 / scatterRate;
    vKick = hbar * kA / mass;

    timeround = round(maxTime / tKick); 

    rend=zeros(1, numParticles);
    vend=zeros(1, numParticles);
    Zero = zeros(1, numParticles);
    finalTemperatures = zeros(1, numParticles);
    finalSizes = zeros(1, numParticles);
    population = zeros(1, numParticles);
    finalSizes_center = zeros(1,numParticles);
    rmatrix= zeros(timeround+1,numParticles);
    vmatrix = zeros(timeround+1,numParticles);
    
    
    for p = 1:numParticles
    
        if CMOT == 1
            initPos = randn * 0.25;   % Random initial position (mm)
            initVel = randn * 0.265;  % Random initial velocity (mm/ms)
        
        elseif CMOT == 2
            initPos = randn * 0.5;
            initVel = randn * 0.839;
        
        elseif CMOT == 0
            initPos = randn*0.04;
            initVel = randn*0.119;
        
        else
            if p<numParticles/2
                initPos = initPosset;
            else
                initPos=-initPosset;
            end
            initVel = randn*0.839;
            %-sign(initPos) * abs(randn * 0.767);
        end

        
        r= zeros(1, timeround+1);
        t= zeros(1, timeround+1);
        v= zeros(1, timeround+1);

        r(1) = initPos;
        v(1) = initVel;
        t(1) = 0;
        for i = 1:timeround
            randPhi1 = 2 * pi * rand;
            randPhi2 = 2 * pi * rand;
            t(i+1) = t(i)+ tKick;
            v(i + 1) = v(i) + vKick * (cos(randPhi1) + cos(randPhi2)) + oneAxisAccel(r(i), v(i)) * tKick;
            r(i + 1) = r(i) + v(i) * tKick;
        end
        
        if abs(r(i)) < 3
            figure(1)
            hold on;
            plot(r, v, '-', 'LineWidth', 2, 'DisplayName',num2str(p)); % Thicker line and legend entry
            figure(5)
            hold on;
            plot(t,r,'-', 'LineWidth', 2, 'DisplayName',num2str(p))
            figure(3)
            hold on;
            plot(t,v,'-', 'LineWidth', 2, 'DisplayName',num2str(p))
            startInd = round(maxTime*(1-1/10)/tKick);
            vSq = mean(v(startInd:end).^2);
            sigma = sqrt(mean((r(startInd:end) - mean(r(startInd:end))).^2));
            sigma_center = sqrt(mean(r(startInd:end).^2));
            temp = vSq * mass / kb * 1e3;

            rmatrix(:,p)=r;
            vmatrix(:,p)=v;        
            rend(p)=mean(r(startInd:end));
            vend(p)=sqrt(vSq);
            finalTemperatures(p) = temp;
            finalSizes(p) = sigma;
            finalSizes_center(p)= sigma_center; 
            population(p) = oneAxisExppop(sigma,sqrt(vSq));
        else
            rend(p)=NaN;
            vend(p)=NaN;
            finalTemperatures(p) = NaN;
            finalSizes(p) = NaN;
            finalSizes_center(p)= NaN; 
            population(p) = NaN;

        end
    end

    figure(1);
    grid on;
    xlabel('d (mm)');
    ylabel('v_{d} (mm/ms)');
    %legend('Location', 'best');

    figure(5);
    grid on;
    xlabel('t (ms)');
    ylabel('d(mm)');
    %legend('Location', 'best')

    figure(3);
    grid on;
    xlabel('t (ms)');
    ylabel('v_{d} (mm/ms)');
    %legend('Location', 'best');

     %write to csv
    if SavePlot==1
        figure(3);
        output_folder_traj = strcat('SnCBMOTtestPlot','/ManyParticleAnalysis_Trajectory/','/',Date,'/');
        if ~exist(output_folder_traj, 'dir')
            mkdir(output_folder_traj);
        end
        % Save the figure as a 300 DPI image
        %saveas(gcf, strcat('SnPlot','/particle_trajectory_MOT/',data_timestamp,MOTCondition,'.png'));
        print(gcf, strcat(output_folder_traj,data_timestamp,MOTCondition,'HDPI'),'-dpng', '-r300');% Save as HDPI Png

        figure(5);
        output_folder_traj = strcat('SnCBMOTtestPlot','/ManyParticleAnalysis_positiontimeplot/','/',Date,'/');
        if ~exist(output_folder_traj, 'dir')
            mkdir(output_folder_traj);
        end
        % Save the figure as a 300 DPI image
        %saveas(gcf, strcat('SnPlot','/particle_trajectory_MOT/',data_timestamp,MOTCondition,'.png'));
        print(gcf, strcat(output_folder_traj,data_timestamp,MOTCondition,'HDPI'),'-dpng', '-r300');% Save as HDPI Png

        figure(1);
        output_folder_traj = strcat('SnCBMOTtestPlot','/ManyParticleAnalysis_velocitytimeplot/','/',Date,'/');
        if ~exist(output_folder_traj, 'dir')
            mkdir(output_folder_traj);
        end
        % Save the figure as a 300 DPI image
        %saveas(gcf, strcat('SnPlot','/particle_trajectory_MOT/',data_timestamp,MOTCondition,'.png'));
        print(gcf, strcat(output_folder_traj,data_timestamp,MOTCondition,'HDPI'),'-dpng', '-r300');% Save as HDPI Png


    end

    figure(4)
    hold on;
    scatter(rend, vend, 'filled');  % 'filled' 
    xlabel('d (mm)');
    ylabel('v_{d} (m/s)');
    %title('Particle Trajectory for Sn atoms in FreeSpace SubCooling for 100 ms');


    %write to csv
    if SavePlot==1
        figure(4);
        output_folder_traj = strcat('SnCBMOTtestPlot','/particle_trajectory_FinalPosition/','/',Date,'/');
        if ~exist(output_folder_traj, 'dir')
            mkdir(output_folder_traj);
        end
        % Save the figure as a 300 DPI image
        %saveas(gcf, strcat('SnPlot','/particle_trajectory_MOT/',data_timestamp,MOTCondition,'.png'));
        print(gcf, strcat(output_folder_traj,data_timestamp,MOTCondition,'HDPI'),'-dpng', '-r300');% Save as HDPI Png
    end

    Lossatom = finalTemperatures(isnan(finalTemperatures));
    
    finalTemperatures = finalTemperatures(~isnan(finalTemperatures));
    finalSizes = finalSizes(~isnan(finalSizes));
    population = population(~isnan(population));
    finalSizes_center = finalSizes_center(~isnan(finalSizes_center));


    outliers = isoutlier(finalSizes_center,'mean');
    cleanfinalTemperatures = finalTemperatures(~outliers);
    cleanfinalSizes = finalSizes(~outliers);
    cleanpopulation = population(~outliers);
    cleanfinalSizes_center = finalSizes_center(~outliers);
    cleandata= finalSizes(outliers);

    

    fprintf('Loss Atom Number: %d\n', numel(Lossatom));

    
    fprintf('Clean sizeData: ');
    disp(cleandata)

    avgTemp = mean(cleanfinalTemperatures);
    stdTemp = std(cleanfinalTemperatures);
    avgSize_center = mean(cleanfinalSizes_center);
    stdSize_Center = std(cleanfinalSizes_center);
    maxSize_Center = max(cleanfinalSizes_center);
    avgSize = mean(cleanfinalSizes);
    maxSize = max(cleanfinalSizes);
    stdSize = std(cleanfinalSizes);
    avgexcpop = mean(population);
    stdexcpop = std(population);
    
    Particle_M=[avgTemp, stdTemp, maxSize, stdSize;rend',vend',Zero',Zero'];
    PosvsTime = [t',rmatrix];
    VelvsTime = [t',vmatrix];
    writematrix(Particle_M, strcat('SnMOTPaper/MOTTempAndSize_ManyParticle_', data_timestamp, '.csv'));
    writematrix(PosvsTime,strcat('SnMOTPaper/CBMOT_PosvsTime', data_timestamp, num2str(numParticles),'.csv'));
    writematrix(VelvsTime,strcat('SnMOTPaper/CBMOT_VelvsTime', data_timestamp, num2str(numParticles),'.csv'));

end
