#1)
#move Julia terminal to directory where this file is located.
#the directory should also contain auxFunctionsSnLatest and SnVariables files.
cd(@__DIR__);
include("SnVariables.jl")
include("auxFunctionsSnMOTLatestModifiedFinal.jl");

using Statistics
using LinearAlgebra
using SharedArrays

# --- m_J labeling: override basis order if needed (default assumes -J:1:J) ---
mJ_order_ground  = nothing           # e.g., [-1, 0, +1]
mJ_order_excited = nothing           # e.g., [-1, 0, +1]

#2) User choices with respect to saving output files
useRandPhase = 1; #if 1, use random phase for each laser. If 0, use same phase for all lasers.
saveInRealUnits = 1;#if 1, save vel+accel in m/s, mm/ms^2.  If 0, save in normalized units (vel= v/(gam/k)), (force=1e-3*hbar*k*gam)
saveData=1;#if you want to save the data
saveDataFolderTag = "SnCompressedRedMOT"; #If you want to put anything additional in "savefoldername" to tag it, see variable folderString after lasers are declared.
addHeaders=1;

#3) Non Laser Detuning/Polarization Simulation Variables (B-field, beam-waist etc.)
bGradReal = 75;# radial (weak axis) gradient in units Gauss/cm.  if "Static", this becomes the static field, in Gauss.
bGrad = (1 / kA * 1e2) * bGradReal; 
beam_waist = 7; #in mm. Since Sn only requires one laser the beam waist will be same for both dual freq components
numTrialsPerValueSet = 250;#number of trials per set of values (displacementsInMM,userSpeeds,longSpeeds)
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

#vRound = 0.02;#nearest normalized unit to round all speeds to.  Simulation will have periodicity of 2*pi/vRound, so if you round finely, will take a while.  Usually choose 0.02
vRound = 0.01; #for compressed redMOT
#vRound = 0.001; #for blueMOT

#4A) parameters for quick test of restoring force
#Sn gam/k = 9.63 m/s (1 unit of velocity)

#=
#for capture redMOT
displacementsInMM = [1,2,3,4,5,6,7,8,9,10,11,12,13,14];
userSpeeds = [-5,-4,-3,-2.5,-2,-1.5,-1,-.75,-.5,-.25,.25,.5,.75,1,1.5,2,2.5,3,4,5];


#format: [red, blue]
s0 = [1.5, 0.3];
laserEnergy = [-1.0, 0.5];
polSign = [-1, 1];
waists = [beam_waist,beam_waist] .* 1e-3 .* kA;#in normalized units
=#


#for compressed redMOT
displacementsInMM = [0.01, 0.03, 0.08, 0.15, 0.3, 0.5, 0.7, 0.9, 1.2, 1.5, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0];
userSpeeds = [-2.5,-2,-1.5,-1,-.75,-.5,-.4,-.3,-.2,-.1,-.05,.05,.1,.2,.3,.4,.5,.75,1,1.5,2,2.5].*0.4;

#format: [red, blue]
s0 = [0.045, 0.009];
#s0 = [0.15, 0.03];
laserEnergy = [-0.5, 0.25];
#laserEnergy = [-1.0, 0.5];
polSign = [-1,1];
waists = [beam_waist,beam_waist] .* 1e-3 .* kA;#in normalized units


#=
#for CB blue MOT
displacementsInMM = [0.01, 0.03, 0.05, 0.1, 0.3, 0.6, 0.9, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0];
userSpeeds = [-0.5, -0.25, -0.15, -0.095, -0.07, -0.05, -0.025, -0.015, -0.005, 0.005, 0.015, 0.025, 0.05, 0.07, 0.095, 0.15, 0.25, 0.5];

#format: [red, blue]
s0 = [1.5, 1.5];
laserEnergy = [3.0, 3.2];
polSign = [1,-1];
waists = [beam_waist,beam_waist] .* 1e-3 .* kA;#in normalized units
=#


#=
#for blueMOT. flip detuning and polarization simultaneously.
displacementsInMM = [0.01, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0];
userSpeeds = [-1.0, -0.8, -0.6, -0.4, -0.2, -0.15, -0.1, -0.07, -0.04, -0.02, -0.01, -0.005, 0.005, 0.01, 0.02, 0.04, 0.07, 0.1, 0.15, 0.2, 0.4, 0.6, 0.8, 1.0];

#format: [blue, red]
s0 = [1.0, 0.2];
laserEnergy = [3.0, -2.0];
polSign = [1,-1];
waists = [beam_waist,beam_waist] .* 1e-3 .* kA;#in normalized units
=#
#=
#for blue molasses (zero B-field)
displacementsInMM = [0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0];
userSpeeds = [-0.8, -0.4, -0.2, -0.15, -0.1, -0.07, -0.04, -0.02, -0.01, -0.005, 0.005, 0.01, 0.02, 0.04, 0.07, 0.1, 0.15, 0.2, 0.4, 0.8];

#format: [blue, red]
s0 = [0.5];
laserEnergy = [3.0];
polSign = [1];
waists = [beam_waist,beam_waist] .* 1e-3 .* kA;#in normalized units
=#

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

headersBase = ["Speed" "av" "del_av" "ar"  "del_ar" "PFLow" "PFeHigh"];

#save data folder path name
if saveData==1
    folderString = string(@__DIR__,"/SnData/",saveDataFolderTag,"BGradGPerCM", bGradReal,"WaistMM", beam_waist,"Date",Dates.format(now(),"yyyymmdd_HHMM"))
    if useRandPhase == 1
        folderString *= "_WithRandPhase"
    end
    mkpath(folderString)
end

forceVsTime = Array{Array{ComplexF64,2},1}(undef, numTrialsPerValueSet * 2);
forceVsSpeed = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2); #a\dot v/|v|
forceVsPos = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2); #a\dot r/|r|

# --- per-sublevel population storage ---
# Dimensions: (speed_index, run_index, sublevel_index)
popGroundVsSpeed  = Array{Float64,3}(undef, length(userSpeeds), numTrialsPerValueSet * 2, nFLow)
popExcitedVsSpeed = Array{Float64,3}(undef, length(userSpeeds), numTrialsPerValueSet * 2, nFHigh)

# --- m_J label helpers ---
Jg = (nFLow  - 1) / 2
Je = (nFHigh - 1) / 2
mJg = isnothing(mJ_order_ground)  ? collect(-Jg:1.0:Jg) : mJ_order_ground
mJe = isnothing(mJ_order_excited) ? collect(-Je:1.0:Je) : mJ_order_excited
@assert length(mJg) == nFLow  "mJ_order_ground length must equal numZeemanStatesGround = $nFLow"
@assert length(mJe) == nFHigh "mJ_order_excited length must equal numZeemanStatesExcited = $nFHigh"

label_m(m) = (abs(m - round(m)) < 1e-9) ? string(Int(round(m))) : string(Int(round(2*m)),"/2")

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
        tForSteadyState = maximum([10 / currSpeed, 1800]);#obtained by trial and error
        periodLength = 2 * pi / vRound;
        saveTimes = tForSteadyState:0.1:(tForSteadyState+periodLength)#times to record obe solution for force integration
        for i = 1:(numTrialsPerValueSet*2)
            forceVsTime[i] = zeros(length(saveTimes), 3)
        end
        prob = ODEProblem(densityMatrixChangeTerms!, pStart, (0.0, tForSteadyState + periodLength), p)#set up OBE problem to solve

        function prob_func(prob, i, repeat)
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

        # Solve ensemble
        ens_prob = EnsembleProblem(prob, prob_func=prob_func)
        @time sol = solve(ens_prob, Tsit5(), EnsembleThreads(); trajectories=numTrialsPerValueSet * 2, saveat=saveTimes)
        
        #8B) calculate forces and populations for each random R, V trial..
        @time for i = 1:(numTrialsPerValueSet*2)
            currSol = sol[i]
            makeForceVsTime!(forceVsTime[i], currSol.t, currSol.u, lasers,
                             couplingMatrices, [randRxs[i], randRys[i], randRzs[i]],
                             [randVxs[i], randVys[i], randVzs[i]], rp[i])

            forceVsSpeed[j, i] = (randVxs[i] * trapz(currSol.t, forceVsTime[i][:, 1]) +
                                  randVys[i] * trapz(currSol.t, forceVsTime[i][:, 2]) +
                                  randVzs[i] * trapz(currSol.t, forceVsTime[i][:, 3])) / 1e-3 /
                                  sqrt(randVxs[i] .^ 2 + randVys[i] .^ 2 + randVzs[i] .^2) /
                                  (currSol.t[end] - currSol.t[1])

            forceVsPos[j, i] = (randRxs[i] * trapz(currSol.t, forceVsTime[i][:, 1]) +
                                randRys[i] * trapz(currSol.t, forceVsTime[i][:, 2]) +
                                randRzs[i] * trapz(currSol.t, forceVsTime[i][:, 3])) / 1e-3 /
                                sqrt(randRxs[i] .^ 2 + randRys[i] .^ 2 + randRzs[i] .^2) /
                                (currSol.t[end] - currSol.t[1])

            pHighVsSpeed[j, i] = mean(real(tr.([maskFHigh .* v for v in currSol.u])))
            pLowVsSpeed[j, i]  = mean(real(tr.([maskFLow  .* v for v in currSol.u])))

            # per-sublevel populations (time-averaged over saved window)
            @inbounds for k = 1:nFLow
                popGroundVsSpeed[j, i, k] = mean([real(ρ[k, k]) for ρ in currSol.u])
            end
            @inbounds for k = 1:nFHigh
                popExcitedVsSpeed[j, i, k] = mean([real(ρ[nFLow + k, nFLow + k]) for ρ in currSol.u])
            end
        end # i
    end # speeds

    #8C) averages and uncertainties over runs
    forceVsSpeedAvg = dropdims(mean(forceVsSpeed, dims=2), dims=2)
    forceVsSpeedUnc = dropdims(std(forceVsSpeed,  dims=2), dims=2) ./ sqrt(numTrialsPerValueSet * 2)

    forceVsPosAvg   = dropdims(mean(forceVsPos, dims=2), dims=2)
    forceVsPosUnc   = dropdims(std(forceVsPos,  dims=2), dims=2) ./ sqrt(numTrialsPerValueSet * 2)

    pHighVsSpeedAvg = dropdims(mean(pHighVsSpeed, dims=2), dims=2)
    pLowVsSpeedAvg  = dropdims(mean(pLowVsSpeed,  dims=2), dims=2)

    # per-sublevel means across runs
    nRuns = numTrialsPerValueSet * 2
    popGroundVsSpeedAvg  = dropdims(mean(popGroundVsSpeed,  dims=2), dims=2)   # (Nspeed, nFLow)
    popExcitedVsSpeedAvg = dropdims(mean(popExcitedVsSpeed, dims=2), dims=2)   # (Nspeed, nFHigh)

    # --- commented out: per-sublevel standard errors over runs ---
    # popGroundVsSpeedSE   = dropdims(std(popGroundVsSpeed,   dims=2), dims=2) ./ sqrt(nRuns)
    # popExcitedVsSpeedSE  = dropdims(std(popExcitedVsSpeed,  dims=2), dims=2) ./ sqrt(nRuns)

    #8D) convert to real units if applicable and save data
    (forceVsSpeedAvgSaveVals,forceVsSpeedUncSaveVals,forceVsPosAvgSaveVals,forceVsPosUncSaveVals) =
        (forceVsSpeedAvg,forceVsSpeedUnc,forceVsPosAvg,forceVsPosUnc).*(accelFactor*saveInRealUnits+1*(1-saveInRealUnits));
    userSpeedsSaveVals = userSpeeds.*(velFactor*saveInRealUnits+1*(1-saveInRealUnits));

    if saveData==1
        # m_J-labeled headers (means only)
        pgHeaders  = ["P_g(mJ=$(label_m(m)))"  for m in mJg]
        peHeaders  = ["P_e(mJ'=$(label_m(m)))" for m in mJe]

        # --- commented out: uncertainty headers ---
        # dpgHeaders = ["del_P_g(mJ=$(label_m(m)))"  for m in mJg]
        # dpeHeaders = ["del_P_e(mJ'=$(label_m(m)))" for m in mJe]

        headersFull = [headersBase permutedims(pgHeaders) permutedims(peHeaders)]
        # headersFull = [headersBase permutedims(pgHeaders) permutedims(peHeaders) permutedims(dpgHeaders) permutedims(dpeHeaders)]

        open(string(folderString,"/forceVsSpeedDisplacement",displacementsInMM[l],"MM",runType,".dat"),"a") do io
            if addHeaders==1
                writedlm(io, [headersFull; hcat(userSpeedsSaveVals,
                                                forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                                                forceVsPosAvgSaveVals,   forceVsPosUncSaveVals,
                                                pLowVsSpeedAvg,          pHighVsSpeedAvg,
                                                popGroundVsSpeedAvg,     popExcitedVsSpeedAvg
                                                # , popGroundVsSpeedSE,      popExcitedVsSpeedSE
                                                )])
            else
                writedlm(io, hcat(userSpeedsSaveVals,
                                  forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals,
                                  forceVsPosAvgSaveVals,   forceVsPosUncSaveVals,
                                  pLowVsSpeedAvg,          pHighVsSpeedAvg,
                                  popGroundVsSpeedAvg,     popExcitedVsSpeedAvg
                                  # , popGroundVsSpeedSE,      popExcitedVsSpeedSE
                                  ))
            end
        end
    end

end # displacements

laserVarHeaders = ["s0" "energy" "polSign"]
if saveData ==1
    open(string(folderString,"/laserVariables.dat"),"w") do io
        writedlm(io,[laserVarHeaders ; hcat(s0,laserEnergy,polSign)]);
    end
end