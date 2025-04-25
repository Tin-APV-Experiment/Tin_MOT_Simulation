#stuff to do:

##

cd(@__DIR__);#moves julia terminal to directory where this file is.  This directory should have auxFunctions+SrF(or whatever)Variables files as well
include("SnVariables.jl")
include("Demo-auxFunctionsSnNewlatest_Remote.jl");

#2) User choices with respect to saving output files
useRandPhase = 1; #if 1, use random phase for each laser. If 0, use same phase for all lasers.
saveInRealUnits = 1;#if 1, save vel+accel in m/s, mm/ms^2.  If 0, save in normalized units (vel= v/(gam/k)), (force=1e-3*hbar*k*gam)
saveData=1;#if you want to save the data
saveDataFolderTag = "SnPapertest"; #If you want to put anything additional in "savefoldername" to tag it, see variable folderString after lasers are declared.
addHeaders=1;

#3) Non Laser Detuning/Pol Simulation Variables (B-field, beam-waist etc.)
bGradReal = 20# in units Gauss/cm.  if "Static", this becomes the static field, in Gauss\
bGrad = (1 / kA * 1e2) * bGradReal;
beam_waist = 7; #in mm. Since Sn only requires one laser the beam waist will be same for both dual freq components
numTrialsPerValueSet = 10;#number of trials per set of values (displacementsInMM,userSpeeds,longSpeeds)
velDirRelToR = 0;#-1 will choose random values for direction of v,r.  0 will force them to be same direction. 1 forces orthogonal.  2 forces opposite.
forceXY=1; #if '1', will force v, r to go along (x+y)/sqrt(2).  Simulates slowing/trapping of molecules moving along slowing axis in tandem with velDirToR = 2 or 0
if velDirRelToR==-1
    runType = string("Random");#goes into folder name.
elseif velDirRelToR==0
    runType = string("SameDir");#goes into folder name.
elseif velDirRelToR==1
    runType = string("OrthoDir");#goes into folder name.
else
    runType = string("OppDir");#goes into folder name.
end
vRound = 0.02;#nearest normalized unit to round all speeds to.  Simulation will have periodicity of 2*pi/vRound, so if you round finely, will take a while.  Usually choose 0.02
#vRound = 0.01; #for compressed redMOT
#vRound = 0.002; #for blueMOT
#4A) parameters for quick test of restoring force
#Sn gam/k = 9.63 m/s (1 unit of velocity)


#for capture redMOT
#displacementsInMM = [1,2,3,5,7,9,11，14，17];
#displacementsInMM = [1,2,3,5,7,9,11];
#userSpeeds = [-8,-6.5,-5,-3.5,-3,-2.5,-2,-1.5,-1,-.6,-.4,-.2,-.1,.1,.2,.4,.6,1,1.5,2,2.5,3,3.5,5,6.5,8]; 
#userSpeeds = [-4,-3.5,-3,-2.5,-2,-1.5,-1,-.6,-.4,-.2,-.1,.1,.2,.4,.6,1,1.5,2,2.5,3,3.5,4];
#userSpeeds = [-5,-3.5,-3,-2.5,-2,-1.5,-1,-.6,-.4,-.2,-.1,.1,.2,.4,.6,1,1.5,2,2.5,3,3.5,5];

#=for redcapture MOT
#format: [red, blue]
s0 = [1.5,0.3].*1.0;
laserEnergy = [-1,0.5];#[-1.0,0.1]
polSign = [-1,1];
waists = [beam_waist,beam_waist] .* 1e-3 .* kA;#in normalized units
=#

#=
#for compressed redMOT
displacementsInMM = [0.5,1,1.5,2.5,3.5,4.5,5.5];
userSpeeds = [-4,-3,-2.5,-2,-1.5,-1,-.75,-.5,-.25,.25,.5,.75,1,1.5,2,2.5,3,4]*.0.4;

#format: [red, blue]
s0 = [0.05, 0.01];
laserEnergy = [-0.5, 0.25];
polSign = [-1,1];
waists = [beam_waist,beam_waist] .* 1e-3 .* kA;#in normalized units
=#


#for simple blueMOT. flip detuing and polarization simultaneously.
#displacementsInMM = [1,2,3,5,7,9,11];#RedMOT Setting
#displacementsInMM=[0.1,0.3,0.5,0.7,0.9,1.5,2.0]; ##Blue MOT Setting
displacementsInMM=[0.05,0.15,0.25,0.45,0.65,0.85,1.0,1.5]; ##Blue MOT Setting
#displacementsInMM=[0.01,0.02,0.03,0.05,0.075,0.1,0.3,0.6,1.0,1.5]; ##Blue MOT Setting
#displacementsInMM=[0.01,0.03,0.05,0.1,0.3,0.6,1.0]; ##Blue MOT Setting
#displacementsInMM=[0.01,0.03,0.05,0.1,0.3,0.6,1.0,1.5,3.0]; ##Blue MOT Setting
#displacementsInMM=[0.1,0.3,0.5,0.7,0.9,1.5,2.0]; ##FreeSpace_Sub Doppler Cooling Setting
#displacementsInMM = [0.5,1,1.5,2,2.5,3,5];

#userSpeeds times 10 approximately change to m/s
#userSpeeds = [-3,-2.5,-2,-1.5,-1,-.75,-.5,-.25,.25,.5,.75,1,1.5,2,2.5,3]; ## RedMOT Setting
#userSpeeds = [-15,-10,-3,-2.4,-1.9,-1.4,-.9,-.5,-.3,-.1,-.05,-.02,.02,.05,.1,.3,.5,.9,1.4,1.9,2.4,3,5,10,15].*0.05;#BlueMOT Setting
userSpeeds = [-3,-2.5,-2,-1.5,-1,-.75,-.5,-.25,-.1,.1,.25,.5,.75,1,1.5,2,2.5,3].*0.05;#FreeSpace_Sub Doppler Cooling Setting

#format: [red, blue]
s0 = [1.5, 1.5]*1;## Blue MOT Seeting
#s0 = [0.3, 0.3]*1;## Lambda-Setting (Zero-Bfield)
#detu=[1,1]*0.5;
#laserEnergy = [-3.5, -3.5]*-1.0;## Lambda-Setting (Zero-Bfield)
laserEnergy = [3.0, 3.2]*1.0;## Blue MOT Setting
polSign = [1,-1]; ## Blue/Lambda Seeting
waists = [beam_waist,beam_waist] .* 1e-3 .* kA;#in normalized units


#other settings:
#displacementsInMM = [1,2,3,5,7,9,11,14,17];
#userSpeeds = [-5,-3.5,-3,-2.5,-2,-1.5,-1,-.6,-.4,-.2,-.1,.1,.2,.4,.6,1,1.5,2,2.5,3,3.5,5];
#displacementsInMM = [1,2,3,4,5,6,7,8,9,10,11,12];
#userSpeeds = [-5,-3.5,-3,-2.5,-2,-1.5,-1,-.6,-.4,-.2,-.1,.1,.2,.4,.6,1,1.5,2,2.5,3,3.5,5];

#=
#for compressed redMOT
displacementsInMM = [1,2,3,5,7,9,11,14,17]*0.4;
userSpeeds = [-5,-3.5,-3,-2.5,-2,-1.5,-1,-.6,-.4,-.2,-.1,.1,.2,.4,.6,1,1.5,2,2.5,3,3.5,5]*0.1;

#format: [red, blue]
s0 = [1.0,0.2].*0.1;
laserEnergy = [-0.8, 0.15];
polSign = [-1,1];
waists = [5,5] .* 1e-3 .* kA;#in normalized units
=#

(couplingMatrices,bCouplingMatrices,numZeemanStatesGround,numZeemanStatesExcited) = createCouplingTermsandLaserMasks()

lasers = Lasers(s0,laserEnergy,polSign,waists);#define lasers structure, see auxFunctions

numZeemanStatesTotal = numZeemanStatesGround + numZeemanStatesExcited;

#everything below here you is part of the initizliation block that you probably won't want to change

rInit = [0.0, 0.0, 0];#placehold not used
vInit = [3.0, 0.0, 0.1];#placehold not used
rpInit = zeros(length(s0),6);
#note: p will include a lot of pre-allocated stuff.  This is basically all of the stuff 'passed to' the obe solver, in addition to the initial condition of the density matrix defined below
#in retrospect p is not the best choice for the variable name but it's the julia house style...maybe replace later. (actually you can't.  Julia forces ODEProblem to have a variable 'p')
pPreInitialized = preInitializer(length(s0),numZeemanStatesGround,numZeemanStatesTotal)


p = [rInit, vInit, rpInit, lasers, bGrad * normalizedBohrMag,
    couplingMatrices[1], couplingMatrices[2], couplingMatrices[3], bCouplingMatrices[1], bCouplingMatrices[2], bCouplingMatrices[3]];
append!(p,pPreInitialized);

#density matrix Initialization
pStart = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);
#pStart[1,1] = 1/2;
#pStart[2,2] = 1/2;
pStart[1,1] = 1/3;
pStart[2,2] = 1/3;
pStart[3,3] = 1/3;

#initialize some 'masks' that zero out subset of population values...helpful for quick calculation of populations in various ground states
nFLow = numZeemanStatesGround;
nFHigh = numZeemanStatesExcited;

maskFHigh = zeros(numZeemanStatesTotal, numZeemanStatesTotal);
maskFLow = zeros(numZeemanStatesTotal, numZeemanStatesTotal);
maskFLow[1:nFLow, 1:nFLow] .= ones(nFLow, nFLow);
maskFHigh[(1+nFLow):(nFLow+nFHigh), (1+nFLow):(nFLow+nFHigh)] .= ones(nFHigh, nFHigh);
pHighVsSpeed = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2);
pLowVsSpeed = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2);
headers = ["Speed" "av" "del_av" "ar"  "del_ar" "PFLow" "PFeHigh"];

if saveData==1
    parent_dir = dirname(@__DIR__)  # go one level up from the script location
    date_str = Dates.format(now(), "yyyymmdd")
    time_str = Dates.format(now(), "yyyymmdd_HHMM")
    folderString = string(parent_dir,"/SnMOT_Data/Blue_MOT/Conveyor_Belt_MOT/",date_str,"/",saveDataFolderTag,"BGradGPerCM", bGradReal,time_str,"Trail",numTrialsPerValueSet)
    if useRandPhase == 1
        folderString *= "_WithRandPhase"
    end
    mkpath(folderString)
end

#below still needs to be edited (TO DO)


forceVsTime = Array{Array{ComplexF64,2},1}(undef, numTrialsPerValueSet * 2);
forceVsSpeed = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2); #a\dot v/|v|
forceVsPos = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2); #a\dot r/|r|

##

#8) Iterate over user choices for displacements and speeds

for l = 1:length(displacementsInMM)
    currDisp = displacementsInMM[l];
    for j = 1:length(userSpeeds)
        currSpeed = userSpeeds[j];
        if abs(currSpeed)<0.04
            vRound = 0.005;
        elseif abs(currSpeed)<0.1
            vRound = 0.01;
        elseif abs(currSpeed)<0.5
            vRound = 0.02;
        else
            vRound = 0.05;
        end
        #8A) Set up and solve OBEs
        (randRxs,randRys,randRzs,randVxs,randVys,randVzs,rp) = generateRandPosAndVel(numTrialsPerValueSet,velDirRelToR,currDisp,currSpeed,vRound,forceXY,length(s0));
        tForSteadyState = maximum([10 / currSpeed, 1800]);#obtained by trial and error.  Could potentially be handled more rigrorously (solve ode in steps of 'period length' until solution 'converges')
        periodLength = 2 * pi / vRound;
        saveTimes = tForSteadyState:0.1:(tForSteadyState+periodLength)#times to record obe solution for force integration
        for i = 1:(numTrialsPerValueSet*2)
            forceVsTime[i] = zeros(length(saveTimes), 3)
        end
        prob = ODEProblem(densityMatrixChangeTerms!, pStart, (0.0, tForSteadyState + periodLength), p)#set up OBE problem to solve

        function prob_func(prob, i, repeat)#change position and velocity of funtion each 'ensemble' sample based on the randomly chosen values for pos/vel vector (magnitude fixed)
            prob.p[1][1] = randRxs[i]
            prob.p[1][2] = randRys[i]
            prob.p[1][3] = randRzs[i]
            prob.p[2][1] = randVxs[i]
            prob.p[2][2] = randVys[i]
            prob.p[2][3] = randVzs[i]

            for laserVar = 1:length(s0)
                prob.p[3][laserVar,1] = rp[i][laserVar,1] .* useRandPhase;
                prob.p[3][laserVar,2] = rp[i][laserVar,2] .* useRandPhase;
                prob.p[3][laserVar,3] = rp[i][laserVar,3] .* useRandPhase;
                prob.p[3][laserVar,4] = rp[i][laserVar,4] .* useRandPhase;
                prob.p[3][laserVar,5] = rp[i][laserVar,5] .* useRandPhase;
                prob.p[3][laserVar,6] = rp[i][laserVar,6] .* useRandPhase;
            end

            remake(prob)
        end

        #these two lines here actually handle the parallized runs of the ode solver
        ens_prob = EnsembleProblem(prob, prob_func=prob_func)#solve obe problem for various initial conditions re-set by 'prob_func' each iteration
        @time sol = solve(ens_prob, Tsit5(), EnsembleThreads(); trajectories=numTrialsPerValueSet * 2, saveat=saveTimes)#parallelized OBE solver, runs on amount of threads made available by CPU (Threads.nthreads())
        
        #8B) calculate forces (f\dot r/|r|, etc.) for each random R, V trial..
        @time for i = 1:(numTrialsPerValueSet*2)
            currSol = sol[i]
            makeForceVsTime!(forceVsTime[i], currSol.t, currSol.u, lasers,
            couplingMatrices, [randRxs[i], randRys[i], randRzs[i]], [randVxs[i], randVys[i], randVzs[i]], rp[i])
            
            forceVsSpeed[j, i] = (randVxs[i] * trapz(currSol.t, forceVsTime[i][:, 1]) +
            randVys[i] * trapz(currSol.t, forceVsTime[i][:, 2]) +
            randVzs[i] * trapz(currSol.t, forceVsTime[i][:, 3])) / 1e-3 / sqrt(randVxs[i] .^ 2 + randVys[i] .^ 2 + randVzs[i] .^2) / (currSol.t[end] - currSol.t[1])
            forceVsPos[j, i] = (randRxs[i] * trapz(currSol.t, forceVsTime[i][:, 1]) +
            randRys[i] * trapz(currSol.t, forceVsTime[i][:, 2]) +
            randRzs[i] * trapz(currSol.t, forceVsTime[i][:, 3])) / 1e-3 / sqrt(randRxs[i] .^ 2 + randRys[i] .^ 2 + randRzs[i] .^2) / (currSol.t[end] - currSol.t[1])

            pHighVsSpeed[j, i] = mean(real(tr.([maskFHigh .* v for v in currSol.u])))
            pLowVsSpeed[j, i] = mean(real(tr.([maskFLow .* v for v in currSol.u])))
        end#for all trials
    end#for speeds

    #8C) for given set of speeds, for current choices of longSpeed and displacement, average a\dot v, a\dot r, populations,etc. over runs
    forceVsSpeedAvg = mean(forceVsSpeed, dims=2)
    forceVsSpeedAvg = dropdims(forceVsSpeedAvg, dims=(2))#converts to vector
    forceVsSpeedUnc = std(forceVsSpeed, dims=2) ./ sqrt(numTrialsPerValueSet * 2)
    forceVsSpeedUnc = dropdims(forceVsSpeedUnc, dims=(2))

    forceVsPosAvg = mean(forceVsPos, dims=2)
    forceVsPosAvg = dropdims(forceVsPosAvg, dims=(2))
    forceVsPosUnc = std(forceVsPos, dims=2) ./ sqrt(numTrialsPerValueSet * 2)
    forceVsPosUnc = dropdims(forceVsPosUnc, dims=(2))

    pHighVsSpeedAvg = mean(pHighVsSpeed, dims=2)
    pHighVsSpeedAvg = dropdims(pHighVsSpeedAvg, dims=(2))
    pLowVsSpeedAvg = mean(pLowVsSpeed, dims=2)
    pLowVsSpeedAvg = dropdims(pLowVsSpeedAvg, dims=(2))

    #8D) convert to real units if applicable and save data
    (forceVsSpeedAvgSaveVals,forceVsSpeedUncSaveVals,forceVsPosAvgSaveVals,forceVsPosUncSaveVals) = (forceVsSpeedAvg,forceVsSpeedUnc,forceVsPosAvg,forceVsPosUnc).*(accelFactor*saveInRealUnits+1*(1-saveInRealUnits));
    userSpeedsSaveVals = userSpeeds.*(velFactor*saveInRealUnits+1*(1-saveInRealUnits));
    if saveData==1
        open(string(folderString,"/forceVsSpeedDisplacement",displacementsInMM[l],"MM",runType,".dat"),"a") do io
            if addHeaders==1
                writedlm(io, [headers; hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                            forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pLowVsSpeedAvg, pHighVsSpeedAvg)])
                    

            else #if you've already added headers/don't want them, just append the current forceVsSpeed to the relevant file (so, if you have different longSpeeds, they'll all show up in same file since file is distinguished by displacement)
                writedlm(io, hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                            forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pLowVsSpeedAvg, pHighVsSpeedAvg))
            end
        end
    end

end#for displacements

laserVarHeaders = ["s0" "energy" "polSign"]
if saveData ==1
    open(string(folderString,"/laserVariables.dat"),"w") do io
        writedlm(io,[laserVarHeaders ; hcat(s0,laserEnergy,polSign)]);
    end
end

