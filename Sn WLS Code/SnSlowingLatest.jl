#1)
#move Julia terminal to directory where this file is located.
#the directory should also contain auxFunctionsSnLatest and SnVariables files.
cd(@__DIR__);
include("SnVariables.jl")
include("auxFunctionsSnLatestModifiedSlowing.jl");

#2) User choices with respect to saving output files
useRandPhase = 0; #if 1, use random phase for each laser. If 0, use same phase for all lasers.
saveInRealUnits = 1;#if 1, save vel+accel in m/s, mm/ms^2.  If 0, save in normalized units (vel= v/(gam/k)), (force=1e-3*hbar*k*gam)
saveData=1;#if you want to save the data
saveDataFolderTag = "SnWhiteLightSlowing"; #If you want to put anything additional in "savefoldername" to tag it, see variable folderString after lasers are declared.
addHeaders=1;

#3) Non Laser Detuning/Polarization Simulation Variables (B-field, beam-waist etc.)
bGradReal = 15;# radial (weak axis) gradient in units Gauss/cm.  if "Static", this becomes the static field, in Gauss.
beam_waist = 7; #in mm. 
numTrialsPerValueSet = 1;#number of trials per set of values (displacementsInMM,userSpeeds,longSpeeds)
velDirRelToR = -1;#-1 will choose random values for direction of v,r.  0 will force them to be same direction. 1 forces orthogonal.  2 forces opposite. 
#note that won't matter for pure WLS + push beam simulation (longitudinal) since for those we don't assume any transverse velocity component.
forceXY=1; #if '1', will force v, r to go along (x+y)/sqrt(2).  Simulates slowing/trapping of molecules moving along slowing axis in tandem with velDirToR = 0
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

#4A) parameters for quick test of restoring force
#Sn gam/k = 9.63 m/s (1 unit of velocity)

#4F) typical choices for simulating slowing

longSpeeds = vcat(-0.5, -0.25, 0.25, 0.5, 1:35);#longitudinal speed for 2d force profile (normalized units v/(gam/k))
userSpeeds = [0.]; #userSpeeds is for transverse velocity component in TwoD force profile. For slowing, let's assume no mean transverse velocity (make it purely 1-D) and incorporate transverse spread later.
displacementsInMM = [0.01, 0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 9.0, 10.0, 12.0, 14.0, 17.0]; #sets transverse circle radius (centered about beam axis) from which to select random initial positions (positions are on this circle, not inside). can't be zero.
#displacementsInMM = [0.01];
forceProfile = "TwoD";
bFieldSetting = "Static";

if forceProfile == "ThreeD"
    headers = ["Speed" "av" "del_av" "ar"  "del_ar" "PFLow" "PFeHigh"];
else
    #if we care about transverse cooling simulations/nonzero transverse velocity component in simulation here
    #headers = ["Speed" "av" "del_av" "ar"  "del_ar"  "LongSpeed" "az" "del_az" "PFLow" "PFeHigh"];
    #if we only care about longitudinal slowing simulation (WLS + push) and study transverse pluming in separate Monte Carlo simulations
    headers = ["LongSpeed" "az" "del_az" "PFLow" "PFeHigh"];
end

#5) User choices for laser parameters (detuning, polarization, etc) example laser values (these all work for SrF).

#= polType can be "3D" (sig +/-, with z-axis (quadrupole coil axis) reversed wrt other axes), "2DSS" (sig +/- but lasers only in x,y direction.  if \sig+ along +x then \sig- along +y).  
"2DPar"(lasers in x,y direction both polarized along z).  "2DPerp" (x laser polarized along y, y polarized along z).  "Slower" (z laser linearly polarized along x) =#

#5G) Slowing with push (second beam is push beam. two stacked EOMs by default).

s0 = [2.0, 0.1];
laserEnergy = [-11.0, 0.25];
polSign = [1,1];#doesn't matter here
polType = ["Slower", "Push"];
sidebandFreqs = [1.5, 0.];#units of \Gamma
sidebandAmps = [1.5, 0.];#radians. 2.4 rad is enough to null out the carrier.
sidebandFreqs2 = [3.1,0.];#second frequency for stacked EOM configuration
sidebandAmps2 = [1.5,0.];#second amplitude for stacked EOM configuration
waists = [beam_waist * 1e-3 * kA, beam_waist * 1e-3 * kA];#in mm. 
staticBFieldAngle = 54.7;#angle between static B-field (for dark state remixing) and slowing laser polarization, in degrees.


#5H) Slowing without push. two stacked EOMs by default.
#=
s0 = [0.77];
laserEnergy = [-11.0];
polSign = [1];#doesn't matter here
polType = ["Slower"];
sidebandFreqs = [1.5];#units of \Gamma
sidebandAmps = [1.5];#radians. 2.4 rad is enough to null out the carrier.
sidebandFreqs2 = [3.1];#second frequency for stacked EOM configuration
sidebandAmps2 = [1.5];#second amplitude for stacked EOM configuration
waists = [beam_waist * 1e-3 * kA];#in mm.  
staticBFieldAngle = 54.7;#angle between static B-field (for dark state remixing) and slowing laser polarization, in degrees.
=#
#5I) Slowing without push. only one EOM.
#=
s0 = [2.];
laserEnergy = [-9.0];
polSign = [1];#doesn't matter here
polType = ["Slower"];
sidebandFreqs = [0.8];#units of \Gamma
sidebandAmps = [1.5];#radians. 2.4 rad is enough to null out the carrier.
sidebandFreqs2 = [0.];#second frequency for stacked EOM configuration
sidebandAmps2 = [0.];#second amplitude for stacked EOM configuration
waists = [beam_waist * 1e-3 * kA];#in mm.  
staticBFieldAngle = 54.7;#angle between static B-field (for dark state remixing) and slowing laser polarization, in degrees.
=#
(couplingMatrices,bCouplingMatrices,numZeemanStatesGround,numZeemanStatesExcited) = createCouplingTermsandLaserMasks()

lasers = Lasers(s0,laserEnergy,polSign,waists,polType,sidebandFreqs,sidebandFreqs2,sidebandAmps,sidebandAmps2,staticBFieldAngle);#define lasers structure, see auxFunctions
#set bGrad (units Gauss/wavenumber) (or make "bGrad" static, in units Gauss)
if bFieldSetting!="Static"
    bGrad = (1 / kA * 1e2) * bGradReal; #  prefactor converts Gauss/cm to Gauss/wavenumber (or Gauss if "Static")
else
    bGrad = bGradReal;
end

numZeemanStatesTotal = numZeemanStatesGround + numZeemanStatesExcited;

#everything below here you is part of the initizliation block that you probably won't want to change

rInit = [0.0, 0.0, 0.0];#placehold not used
vInit = [0.0, 0.0, 0.0];#placehold not used
rpInit = zeros(length(s0),6);
#note: p will include a lot of pre-allocated stuff.  This is basically all of the stuff 'passed to' the obe solver, in addition to the initial condition of the density matrix defined below
#in retrospect p is not the best choice for the variable name but it's the julia house style...maybe replace later. (actually you can't.  Julia forces ODEProblem to have a variable 'p')
pPreInitialized = preInitializer(length(s0),numZeemanStatesGround,numZeemanStatesTotal)


p = [rInit, vInit, rpInit, lasers, bGrad * normalizedBohrMag,
    couplingMatrices[1], couplingMatrices[2], couplingMatrices[3], bCouplingMatrices[1], bCouplingMatrices[2], bCouplingMatrices[3]];
append!(p,pPreInitialized);
push!(p,bFieldSetting);

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

#save data folder path name
if saveData==1
    folderString = string(@__DIR__,"/SnData/",saveDataFolderTag,"BFieldGauss", bGradReal,"WaistMM", beam_waist,"Date",Dates.format(now(),"yyyymmdd_HHMM"))
    if useRandPhase == 1
        folderString *= "_WithRandPhase"
    end
    mkpath(folderString)
end

forceVsTime = Array{Array{ComplexF64,2},1}(undef, numTrialsPerValueSet * 2);
forceVsSpeed = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2); #a\dot v/|v|
forceVsPos = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2); #a\dot r/|r|
if forceProfile == "TwoD"
    forceVsLong = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2);#az
end

#8) Iterate over user choices for displacements and speeds

for l = 1:length(displacementsInMM)
    currDisp = displacementsInMM[l];
    for k=1:length(longSpeeds)
        currLongSpeed = longSpeeds[k];
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
            (randRxs,randRys,randRzs,randVxs,randVys,randVzs,rp) = generateRandPosAndVel(forceProfile,numTrialsPerValueSet,velDirRelToR,currDisp,currSpeed,vRound,currLongSpeed,forceXY,length(s0));
            #tForSteadyState = maximum([10 / currSpeed, 1800]);#obtained by trial and error.  Could potentially be handled more rigrorously (solve ode in steps of 'period length' until solution 'converges')
            tForSteadyState = 1800; #only do this for pure WLS when assuming no transverse velocity
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
                if forceProfile=="TwoD"
                    forceVsSpeed[j, i] = (randVxs[i] * trapz(currSol.t, forceVsTime[i][:, 1]) + randVys[i] * trapz(currSol.t, forceVsTime[i][:, 2])) / 1e-3 / sqrt(randVxs[i] .^ 2 + randVys[i] .^ 2) / (currSol.t[end] - currSol.t[1])
                    forceVsPos[j, i] = (randRxs[i] * trapz(currSol.t, forceVsTime[i][:, 1]) + randRys[i] * trapz(currSol.t, forceVsTime[i][:, 2])) / 1e-3 / sqrt(randRxs[i] .^ 2 + randRys[i] .^ 2) / (currSol.t[end] - currSol.t[1])
                    forceVsLong[j,i] = trapz(currSol.t, forceVsTime[i][:, 3]) / 1e-3 / (currSol.t[end] - currSol.t[1]);
                else
                    forceVsSpeed[j, i] = (randVxs[i] * trapz(currSol.t, forceVsTime[i][:, 1]) +
                    randVys[i] * trapz(currSol.t, forceVsTime[i][:, 2]) +
                    randVzs[i] * trapz(currSol.t, forceVsTime[i][:, 3])) / 1e-3 / sqrt(randVxs[i] .^ 2 + randVys[i] .^ 2 + randVzs[i] .^2) / (currSol.t[end] - currSol.t[1])
                    forceVsPos[j, i] = (randRxs[i] * trapz(currSol.t, forceVsTime[i][:, 1]) +
                    randRys[i] * trapz(currSol.t, forceVsTime[i][:, 2]) +
                    randRzs[i] * trapz(currSol.t, forceVsTime[i][:, 3])) / 1e-3 / sqrt(randRxs[i] .^ 2 + randRys[i] .^ 2 + randRzs[i] .^2) / (currSol.t[end] - currSol.t[1])
                end
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
        if forceProfile=="TwoD"
            forceVsLongAvg = mean(forceVsLong, dims=2)
            forceVsLongAvg = dropdims(forceVsLongAvg, dims=(2))
            forceVsLongUnc = std(forceVsLong, dims=2) ./ sqrt(numTrialsPerValueSet * 2)
            forceVsLongUnc = dropdims(forceVsLongUnc, dims=(2))
        end
        pHighVsSpeedAvg = mean(pHighVsSpeed, dims=2)
        pHighVsSpeedAvg = dropdims(pHighVsSpeedAvg, dims=(2))
        pLowVsSpeedAvg = mean(pLowVsSpeed, dims=2)
        pLowVsSpeedAvg = dropdims(pLowVsSpeedAvg, dims=(2))

        #8D) convert to real units if applicable and save data
        (forceVsSpeedAvgSaveVals,forceVsSpeedUncSaveVals,forceVsPosAvgSaveVals,forceVsPosUncSaveVals) = (forceVsSpeedAvg,forceVsSpeedUnc,forceVsPosAvg,forceVsPosUnc).*(accelFactor*saveInRealUnits+1*(1-saveInRealUnits));
        userSpeedsSaveVals = userSpeeds.*(velFactor*saveInRealUnits+1*(1-saveInRealUnits));
        if forceProfile=="TwoD"
            (forceVsLongAvgSaveVals,forceVsLongUncSaveVals) = (forceVsLongAvg,forceVsLongUnc).*(accelFactor*saveInRealUnits+1*(1-saveInRealUnits));
            currLongSpeedSaveVals = currLongSpeed.*(velFactor*saveInRealUnits+1*(1-saveInRealUnits));
        end
        if saveData==1
            open(string(folderString,"/forceVsSpeedDisplacement",displacementsInMM[l],"MM",runType,".dat"),"a") do io
                if addHeaders==1 && k==1
                    if forceProfile=="TwoD"
                        #this is for nonzero transverse velocities (e.g. transverse cooling simulations)
                        #writedlm(io,[headers ; hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals, forceVsPosAvgSaveVals, forceVsPosUncSaveVals,fill(currLongSpeedSaveVals,length(userSpeeds)),forceVsLongAvgSaveVals,forceVsLongUncSaveVals,pLowVsSpeedAvg, pHighVsSpeedAvg)]);
                        #this is for zero transverse velocities (e.g. pure longitudinal WLS + push simulations)
                        writedlm(io,[headers ; hcat(fill(currLongSpeedSaveVals,length(userSpeeds)),forceVsLongAvgSaveVals,forceVsLongUncSaveVals,pLowVsSpeedAvg, pHighVsSpeedAvg)]);
                    else
                        writedlm(io,[headers ; hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals, forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pLowVsSpeedAvg, pHighVsSpeedAvg)]);
                    end

                else #if you've already added headers/don't want them, just append the current forceVsSpeed to the relevant file (so, if you have different longSpeeds, they'll all show up in same file since file is distinguished by displacement)
                    if forceProfile=="TwoD"
                        #this is for nonzero transverse velocities (e.g. transverse cooling simulations)
                        #writedlm(io,hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals, forceVsPosAvgSaveVals, forceVsPosUncSaveVals,fill(currLongSpeedSaveVals,length(userSpeeds)),forceVsLongAvgSaveVals,forceVsLongUncSaveVals,pLowVsSpeedAvg, pHighVsSpeedAvg));
                        #this is for zero transverse velocities (e.g. pure longitudinal WLS + push simulations)
                        writedlm(io,hcat(fill(currLongSpeedSaveVals,length(userSpeeds)),forceVsLongAvgSaveVals,forceVsLongUncSaveVals,pLowVsSpeedAvg, pHighVsSpeedAvg));
                    else
                        writedlm(io,hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals, forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pLowVsSpeedAvg, pHighVsSpeedAvg));
                    end
                end
            end
        end
    end#for longitudinal speeds
end #for displacements

laserVarHeaders = ["s0" "energy" "polSign" "polType" "sidebandFreqs" "sidebandAmps" "sidebandFreqs2" "sidebandAmps2" "staticBFieldAngle"]
if saveData ==1
    if length(s0) == 1 #only one WLS laser
        open(string(folderString,"/laserVariables.dat"),"w") do io
            writedlm(io,[laserVarHeaders ; hcat(s0,laserEnergy,polSign,polType,sidebandFreqs,sidebandAmps,sidebandFreqs2,sidebandAmps2,staticBFieldAngle)]);
        end
    else #WLS laser + push beam
        open(string(folderString,"/laserVariables.dat"),"w") do io
            writedlm(io,[laserVarHeaders ; hcat(s0,laserEnergy,polSign,polType,sidebandFreqs,sidebandAmps,sidebandFreqs2,sidebandAmps2,fill(staticBFieldAngle, 2))]);
        end
    end
end