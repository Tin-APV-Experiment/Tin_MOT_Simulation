using DifferentialEquations
using BenchmarkTools
using Octavian
using LinearAlgebra
using SharedArrays
using WignerSymbols
using Trapz
using DelimitedFiles
using Statistics
using Dates

#=
function landeGJ(J,L,S)
    gj=1+(J*(J+1)+S*(S+1)-L*(L+1))/(2*J*(J+1));
    return gj;
end

function landeGF(F,J,L,I,S)#no hyperfine structure in Sn
    gf=landeGJ(J,L,S);
    return gf;
end
=#

struct Lasers{T1<:Vector{Float64},T2<:Vector{Int64},T3<:Vector{String},T4<:Float64} #structure 'defining' a laser
    s0::T1;#saturation intensity at laser center (single pass)
    laserEnergy::T1;#energy of laser (note: zero energy defined to be energy of transition from s to p_F', where F'=J+I
    polSign::T2;#polarization sign (for configurations using \sigma+/- light.  Defines if x-axis, say, is +\sigma or -\sigma (and corresponding changes to other axes...))
    waists::T1;
    polType::T3;#= polType can be "3D" (sig +/-, with z-axis (quadrupole coil axis) reversed wrt other axes), "2DSS" (sig +/- but lasers only in x,y direction.  if \sig+ along +x then \sig- along +y).  
    "2DPar"(lasers in x,y direction both polarized along z).  "2DPerp" (x laser polarized along y, y polarized along z).  "Slower" (z laser linearly polarized along x) =#
    sidebandFreqs::T1;#frequency at which sidebands are driven
    sidebandFreqs2::T1;#frequency at which sidebands are driven for second EOM (if applicable)
    sidebandAmps::T1;#phase modulation depth in radians
    sidebandAmps2::T1;#phase modulation depth in radians for second EOM (if applicable)
    staticBFieldAngle::T4;#angle of static B field with respect to slowing laser polarization (x-axis) in degrees, only used for WLS simulation
end

function preInitializer(numLasers,numZeemanStatesGround,numZeemanStatesTotal)#initializes a bunch of stuff used in the OBE solver.  Julia likes things pre-initialized if possible

    #holds the modified coupling matrices used in decay terms

    coupleMatEff1 = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);
    coupleMatEff2 = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);
    coupleMatEff3 = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);

    #convenient for fast evaluation of terms used in the 'decay' term of the density matrix evolution (second term in eq 1 of main writeup)

    decayMaskAllButTopLeft = zeros(Float64, numZeemanStatesTotal, numZeemanStatesTotal);
    decayMaskAllButTopLeft[(numZeemanStatesGround+1):numZeemanStatesTotal, (numZeemanStatesGround+1):numZeemanStatesTotal] .= -1;
    decayMaskAllButTopLeft[1:numZeemanStatesGround, (numZeemanStatesGround+1):numZeemanStatesTotal] .= -1 / 2;
    decayMaskAllButTopLeft[(numZeemanStatesGround+1):numZeemanStatesTotal, 1:numZeemanStatesGround] .= -1 / 2;
    decayMaskForCalcTopLeft = zeros(Int64, numZeemanStatesTotal, numZeemanStatesTotal);
    decayMaskForCalcTopLeft[(numZeemanStatesGround+1):numZeemanStatesTotal, (numZeemanStatesGround+1):numZeemanStatesTotal] .= 1;

    #now we make a bunch of initializations.  This makes the julia code run much faster at the cost of some readability...
    r = Array{Float64,1}(undef, 3)

    #fieldTerms[i] are the projections of the light field for laser[i] at a given position on the \sigma^-,\pi,\sigma^+ basis
    fieldTerms = Array{Array{ComplexF64,1},1}(undef, numLasers);
    for i = 1:numLasers
        fieldTerms[i] = vec(zeros(ComplexF64, 1, 3))
    end

    #will eventually 'hold' the atom-light matrix term of the hamiltonian during the diff-eq solver (see densityMatrixChangeTerms! in auxFunctions)
    atomLightTerm = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);

    #bField terms are the projections of magnetic field at a given position on the \sigma^-,\pi,\sigma^+ basis. bFieldTermFull basically holds the 'mu' tensor
    bFieldTerms = Array{ComplexF64,1}(undef, 3);
    bFieldTermFull = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);

    #will eventually hold the -\mu\cdotB (and hermitian conjugate) terms
    uProdBField = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);
    bFieldProdU = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);

    #initializations of some matrices used to speed up the decay term calculation
    decayFull = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);
    pOnlyExcitedStates = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);
    pTopLeft1PreMult = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);
    pTopLeft2PreMult = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);
    pTopLeft3PreMult = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);
    pTopLeft1 = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);
    pTopLeft2 = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);
    pTopLeft3 = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);

    #return pre-initialized stuff from here to join the rest of the pre-initialized stuff in the 'main' program
    pPreInitialized = [coupleMatEff1,coupleMatEff2,coupleMatEff3,decayMaskAllButTopLeft, 
    decayMaskForCalcTopLeft, r, fieldTerms, atomLightTerm, bFieldTerms, bFieldTermFull, uProdBField, bFieldProdU, decayFull,
    pOnlyExcitedStates, pTopLeft1PreMult, pTopLeft2PreMult, pTopLeft3PreMult, pTopLeft1, pTopLeft2, pTopLeft3];
    return pPreInitialized;
end

function generateRandPosAndVel(forceProfile,numTrialsPerSpeed,velDirRelToR,currDisp,currSpeed,vRound,longSpeed,forceXY,numLasers);
    #Function generates a set of random positions and 'pseudo'-random velocities (direction determined by 'velDirRelToR' + whether 'force profile' is 2D or 3D.)
    #if forceProfile is TwoD: z position is assumed to not matter, z velocity is fixed to longSpeed, and direction of velocity relative to random choice of \phi where x=disp*(cos(\phi)), etc. determined by velDirRelToR
    #if forceProfile is ThreeD: longSpeed isn't used, and direction of velocity chosen relative to random x,y,z direction of position is determined by velDirRelToR
    #velDirRelToR=-1 gives random orientation
    #velDirRelToR=0 forces v parallel to r
    #velDirRelToR=1 forces v perpendicular to r (or, for 2D, perpendicular in the xy plane at least)
    #velDirRelToR=2 forces v anti-parallel to r
    #=
    rp = zeros(numLasers,6);
    for laserVar=1:numLasers
        
        
        #comment out this block if don't want random phase to affect makeField terms
        rp[laserVar,1] = rand()*2*pi*1.0;
        rp[laserVar,2] = rand()*2*pi*1.0;
        rp[laserVar,3] = rand()*2*pi*1.0;
        rp[laserVar,4] = rand()*2*pi*1.0;
        rp[laserVar,5] = rand()*2*pi*1.0;
        rp[laserVar,6] = rand()*2*pi*1.0;
        
        
        #=
        #use the below block if you want to force a specific phase
        rp[laserVar,1] = 0;
        rp[laserVar,2] = 0;
        rp[laserVar,3] = 0;
        rp[laserVar,4] = 0;
        rp[laserVar,5] = 0;
        rp[laserVar,6] = 0;
        =#
    end
    =#
    
    # Random phases now have an extra axis for trials: (2*numTrialsPerSpeed × numLasers × 6)
    rp = [zeros(numLasers, 6) for _ in 1:2*numTrialsPerSpeed]
    for trial = 1:2*numTrialsPerSpeed
        randphase1 = rand()*2*pi*1.0;
        randphase2 = rand()*2*pi*1.0;
        randphase3 = rand()*2*pi*1.0;
        randphase4 = rand()*2*pi*1.0;
        randphase5 = rand()*2*pi*1.0;
        randphase6 = rand()*2*pi*1.0;
        for laserVar = 1:numLasers
            #=
            # If we want the dual frequencies to maintain coherence (for lambda cooling or just in general)
            rp[trial][laserVar, 1] = randphase1
            rp[trial][laserVar, 2] = randphase2
            rp[trial][laserVar, 3] = randphase3
            rp[trial][laserVar, 4] = randphase4
            rp[trial][laserVar, 5] = randphase5
            rp[trial][laserVar, 6] = randphase6
            =#
            #=
            # Assign a unique random phase to each laser for each trial. Use this block if you want random phases.
            rp[trial][laserVar, 1] = rand() * 2 * pi
            rp[trial][laserVar, 2] = rand() * 2 * pi
            rp[trial][laserVar, 3] = rand() * 2 * pi
            rp[trial][laserVar, 4] = rand() * 2 * pi
            rp[trial][laserVar, 5] = rand() * 2 * pi
            rp[trial][laserVar, 6] = rand() * 2 * pi
            =#
            
            # Assign a fixed phase to each laser for each trial. Use this block if you want to force a specific phase.
            rp[trial][laserVar, 1] = 0
            rp[trial][laserVar, 2] = 0
            rp[trial][laserVar, 3] = 0
            rp[trial][laserVar, 4] = 0
            rp[trial][laserVar, 5] = 0
            rp[trial][laserVar, 6] = 0
            
        end
    end
    if forceProfile=="TwoD"
        #randomize position direction
        randPhisPos = rand(numTrialsPerSpeed, 1) * 2 * pi
        randRxs = cos.(randPhisPos) .* currDisp .* 1e-3 .* kA
        randRys = sin.(randPhisPos) .* currDisp .* 1e-3 .* kA
        randRzs = rand(numTrialsPerSpeed, 1) * 2 * pi

        if velDirRelToR==-1#randomize phi
            randPhisVels = rand(numTrialsPerSpeed, 1) * 2 * pi
        else#adjust phis from positions to force v either parallel, orthogonal, or anti-parallel to the r choice
            randPhisVels = randPhisPos .+ pi/2*velDirRelToR;#pi/2 for ortho, pi for total reversal, 0 for same.
        end
        # allow for random transverse velocity component, e.g. if we care about transverse cooling simulations
        #=
        randVxs = round.(currSpeed .* cos.(randPhisVels) ./ vRound) .* vRound
        randVys = round.(currSpeed .* sin.(randPhisVels) ./ vRound) .* vRound
        randVzs = round.(currSpeed .* cos.(randPhisVels) ./ vRound) .* 0 .+ longSpeed
        =#
        
        #assume mean transverse velocity is zero so don't allow transverse velocity. The simulation is purely 1-D, longitudinal velocity vs longitudinal deceleration.
        randVxs = zeros(numTrialsPerSpeed)
        randVys = zeros(numTrialsPerSpeed)
        randVzs = fill(longSpeed, numTrialsPerSpeed)
        
    else #if 3D
        #random position direction
        randX = randn(numTrialsPerSpeed, 1);
        randY = randn(numTrialsPerSpeed, 1);
        randZ = randn(numTrialsPerSpeed, 1);
        normTerms = sqrt.(randX.^2 .+ randY.^2 .+ randZ.^2);
        randRxs = randX ./ normTerms .* currDisp .* 1e-3 .* kA;
        randRys = randY ./ normTerms .* currDisp .* 1e-3 .* kA;
        randRzs = randZ ./ normTerms .* currDisp .* 1e-3 .* kA;
        #forceXY forces position to be along (x+y)/sqrt(2) (e.g., entering from slower) (if forceXY==2, then it forces along Z)
        if forceXY == 1
            randRxs = 1 ./ sqrt(2) .* currDisp .* 1e-3 .* kA .+ 2 .* pi .* randn(numTrialsPerSpeed, 1);
            randRys = 1 ./ sqrt(2) .* currDisp .* 1e-3 .* kA .+ 2 .* pi .* randn(numTrialsPerSpeed, 1);
            randRzs = 2 .* pi .* randn(numTrialsPerSpeed, 1);
            randX = randRxs ./ sqrt.(randRxs.^2 .+ randRys.^2 .+ randRzs.^2);
            randY = randRys ./ sqrt.(randRxs.^2 .+ randRys.^2 .+ randRzs.^2);
            randZ = randRzs ./ sqrt.(randRxs.^2 .+ randRys.^2 .+ randRzs.^2);
            normTerms = sqrt.(randX.^2 .+ randY.^2 .+ randZ.^2);
        elseif forceXY == 2
            randRxs = 2 .* pi .* randn(numTrialsPerSpeed, 1);
            randRys = 2 .* pi .* randn(numTrialsPerSpeed, 1);
            randRzs = currDisp .* 1e-3 .* kA .+ 2 .* pi .* randn(numTrialsPerSpeed, 1);
            randX = randRxs ./ sqrt.(randRxs.^2 .+ randRys.^2 .+ randRzs.^2);
            randY = randRys ./ sqrt.(randRxs.^2 .+ randRys.^2 .+ randRzs.^2);
            randZ = randRzs ./ sqrt.(randRxs.^2 .+ randRys.^2 .+ randRzs.^2);
            normTerms = sqrt.(randX.^2 .+ randY.^2 .+ randZ.^2);
        end
        if velDirRelToR == -1#random velocity direction as wel
            randX = randn(numTrialsPerSpeed, 1);#re-roll
            randY = randn(numTrialsPerSpeed, 1);
            randZ = randn(numTrialsPerSpeed, 1);
            normTerms = sqrt.(randX.^2 .+ randY.^2 .+ randZ.^2);
            randVxs = randX ./ normTerms .* currSpeed;
            randVys = randY ./ normTerms .* currSpeed;
            randVzs = randZ ./ normTerms .* currSpeed;
        elseif velDirRelToR == 0#same dir
            randVxs = randX ./ normTerms .* currSpeed;
            randVys = randY ./ normTerms .* currSpeed;
            randVzs = randZ ./ normTerms .* currSpeed;
        elseif velDirRelToR == 1#ortho dir
            randX2 = randn(numTrialsPerSpeed,1);
            randY2 = randn(numTrialsPerSpeed,1);
            randZ2 = randn(numTrialsPerSpeed,1);
        for i=1:length(randX2)
            (randX2[i],randY2[i],randZ2[i]) =[randX2[i],randY2[i],randZ2[i]]-dot([randX[i],randY[i],randZ[i]],[randX2[i],randY2[i],randZ2[i]])./dot([randX[i],randY[i],randZ[i]],[randX[i],randY[i],randZ[i]]) .* [randX[i],randY[i],randZ[i]];
        end
            normTerms = sqrt.(randX2.^2 .+ randY2.^2 .+ randZ2 .^2);
            randVxs = randX2 ./ normTerms .* currSpeed;
            randVys = randY2 ./ normTerms .* currSpeed;
            randVzs = randZ2 ./ normTerms .* currSpeed;
        elseif velDirRelToR == 2#negative dir
            randVxs = -randX ./ normTerms .* currSpeed;
            randVys = -randY ./ normTerms .* currSpeed;
            randVzs = -randZ ./ normTerms .* currSpeed;
        end
        randVxs = round.(randVxs ./ vRound) .* vRound;
        randVys = round.(randVys ./ vRound) .* vRound;
        randVzs = round.(randVzs ./ vRound) .* vRound;
    end
    #run for both +/- r and +/- v (better statistics)
    randRxs = [randRxs; -randRxs];
    randRys = [randRys; -randRys];
    randRzs = [randRzs; -randRzs];
    randVxs = [randVxs; -randVxs];
    randVys = [randVys; -randVys];
    if forceProfile=="ThreeD"
        randVzs = [randVzs; -randVzs];
    else
        randVzs = [randVzs; randVzs];#for 2D, Vz is always "longSpeed"
    end
    #rp = [rp; rp]; #uncomment if want you want the same random phase for both +/- r and +/- v
    #=
    # Print the random positions and velocities as a check
    println(rp)
    println(randRxs)
    =#

    for i=1:length(randVzs)#velocity along any dimension cannot be zero (particle should have x,y,z all change throughout OBE evolution to ensure periodicity)
        if randVxs[i]==0
            randVxs[i] = vRound .* sign.(randn(Float64));
        end
        if randVys[i]==0
            randVys[i] = vRound .* sign.(randn(Float64));
        end
        if randVzs[i]==0
            randVzs[i] = vRound .* sign.(randn(Float64));
        end
    end
    
    return randRxs,randRys,randRzs,randVxs,randVys,randVzs,rp;

end

function createCouplingTermsandLaserMasks()
    numZeemanStatesGround = 3;
    numZeemanStatesExcited = 1;
    numZeemanStatesTotal = numZeemanStatesGround + numZeemanStatesExcited;
    couplingMatrices = Matrix[zeros(numZeemanStatesTotal, numZeemanStatesTotal), zeros(numZeemanStatesTotal, numZeemanStatesTotal), zeros(numZeemanStatesTotal, numZeemanStatesTotal)];

    makeCouplingMatrices!(couplingMatrices);

    bCouplingMatrices = Matrix[zeros(numZeemanStatesTotal, numZeemanStatesTotal), zeros(numZeemanStatesTotal, numZeemanStatesTotal), zeros(numZeemanStatesTotal, numZeemanStatesTotal)];

    makeBCouplingMatrices!(bCouplingMatrices)

    return couplingMatrices,bCouplingMatrices,numZeemanStatesGround,numZeemanStatesExcited;
    
end

function makeCouplingMatrices!(couplingMatrices)
    JGround = 1;
    JExcited = 0;
    nJGround = 2*JGround+1; #number of ground mJ states
    nJExcited = 2*JExcited+1; #number of excited mJ states
    #create column vectors of J and M values for ground and excited states
    JsExc = repeat([JExcited],nJExcited,1);
    MsExc = vcat((-JExcited:1:JExcited));
    JsGround = repeat([JGround],nJGround,1);
    MsGround = vcat((-JGround:1:JGround));

    for pol=-1:1:1
        for row=1:length(JsGround)
            for col=1:length(JsExc)
                if abs(JsExc[col]-JsGround[row])<=1
                    couplingMatrices[pol+2][row,col+length(JsGround)] = sqrt(2*JExcited+1)*(-1)^(JsGround[row]- MsGround[row])*wigner3j(JsGround[row],1,JsExc[col],-MsGround[row],pol,MsExc[col]);
                    #=
                    # Print the value of the coupling matrix element
                    println("couplingMatrices[", pol + 2, "][", row, ", ", col + length(JsGround), "] = ", couplingMatrices[pol + 2][row, col + length(JsGround)]);
                    =#
                end
            end
        end
    end

    couplingMatrices[1][isnan.(couplingMatrices[1])] .= 0
    couplingMatrices[2][isnan.(couplingMatrices[2])] .= 0
    couplingMatrices[3][isnan.(couplingMatrices[3])] .= 0
    #=
    # Print the final coupling matrices
    println("Final couplingMatrices[1]:\n", couplingMatrices[1])
    println("Final couplingMatrices[2]:\n", couplingMatrices[2])
    println("Final couplingMatrices[3]:\n", couplingMatrices[3])
    =#
end

function makeBCouplingMatrices!(bCouplingMatrices)
    nFLow = 3;
    FLow = 1;
    FHigh = 0;
    
    #ground state first
    for m=1:(2*FLow+1)
        for n=1:(2*FLow+1)
            #gGround = landeGF(FLow,JGround,LGround,I,S); this gives 1.5 but is only accurate for 1-electron system with no relativistic correction
            gGround = 1.502;
            wigCalcM = m-mean(1:(2*FLow+1));
            wigCalcN = n-mean(1:(2*FLow+1));
            bCouplingMatrices[1][m,n] = gGround*(-1)^(FLow-wigCalcM-1)*sqrt(FLow*(FLow+1)*(2*FLow+1))*wigner3j(FLow,1,FLow,-wigCalcM,1,wigCalcN);        
            bCouplingMatrices[2][m,n] = gGround*(-1)^(FLow-wigCalcM)*sqrt(FLow*(FLow+1)*(2*FLow+1))*wigner3j(FLow,1,FLow,-wigCalcM,0,wigCalcN);   
            bCouplingMatrices[3][m,n] = gGround*(-1)^(FLow-wigCalcM+1)*sqrt(FLow*(FLow+1)*(2*FLow+1))*wigner3j(FLow,1,FLow,-wigCalcM,-1,wigCalcN);   
        end
    end
    
    #next add excited state
    for m=1:(2*FHigh+1)
        for n=1:(2*FHigh+1)
            gExc = 0; #Sn 3P0 has no gJ factor. 
            wigCalcM = m-mean(1:(2*FHigh+1));
            wigCalcN = n-mean(1:(2*FHigh+1));
            bCouplingMatrices[1][m+nFLow,n+nFLow] = gExc*(-1)^(FHigh-wigCalcM-1)*sqrt(FHigh*(FHigh+1)*(2*FHigh+1))*wigner3j(FHigh,1,FHigh,-wigCalcM,1,wigCalcN);        
            bCouplingMatrices[2][m+nFLow,n+nFLow] = gExc*(-1)^(FHigh-wigCalcM)*sqrt(FHigh*(FHigh+1)*(2*FHigh+1))*wigner3j(FHigh,1,FHigh,-wigCalcM,0,wigCalcN);   
            bCouplingMatrices[3][m+nFLow,n+nFLow] = gExc*(-1)^(FHigh-wigCalcM+1)*sqrt(FHigh*(FHigh+1)*(2*FHigh+1))*wigner3j(FHigh,1,FHigh,-wigCalcM,-1,wigCalcN);   
        end
    end

    #replace NaN with zero 
    bCouplingMatrices[1][isnan.(bCouplingMatrices[1])] .= 0
    bCouplingMatrices[2][isnan.(bCouplingMatrices[2])] .= 0
    bCouplingMatrices[3][isnan.(bCouplingMatrices[3])] .= 0
end

function propR!(r, rInit::Array{Float64,1}, v::Array{Float64,1}, t::Float64)
    r[1] = rInit[1] + v[1] * t;
    r[2] = rInit[2] + v[2] * t;
    r[3] = rInit[3] + v[3] * t;
end

function makeFieldTerms!(fieldTerms, r::Array{Float64,1}, polSign::Array{Int64,1}, polType::Array{String,1}, waists::Array{Float64,1},rp::Matrix{Float64}) #These are different for 2D MOT
    #returns 'field terms' for all lasers.  Field terms depend on the polarization type (and, if \sigma +/-, the sign)
    #This is basically the field for laser [i] due to the 6 (if 3D), or 1 (for Slower/push), or 4 (for all 2D lasers) passes of the beam expressed in the standard \sigma^- ([i][1]), \pi ([i][2]) and \sigma^+ ([i][3])
    #This is calculated in the way illustrated in JOSAB 6(11) 2023-2045 (1989) by Cohen-Tannoudji + Dalibard section 2.  See also Eq15-16 and subsequent expressions in my writeup for the 3D example.

    #IMPORTANT CAVEAT: all terms are 'pre-conjugated' since only the complex conjugate of this term is ever used (Eq 21 of my writeup).  Better to just express it pre-conjugated instead of 
    #repeatedly taking conjugates in the diff-eq solver

    #note also that while 3D MOT fields account for finite beam waist, 1D slower/push field does not (assume uniformly distributed intensity)
    for i=1:length(polSign)#iterate through all lasers
        #=
        #the original way in code with no random phase
        fieldTerms[i][1] = cos(r[3])*exp(-2*(r[1]^2+r[2]^2)/waists[i]^2) .+ polSign[i] * sin(r[1])*exp(-2*(r[2]^2+r[3]^2)/waists[i]^2) .-
        im * (polSign[i] * sin(r[3])*exp(-2*(r[1]^2+r[2]^2)/waists[i]^2) .- cos(r[2])*exp(-2*(r[1]^2+r[3]^2)/waists[i]^2));
        
        fieldTerms[i][2] = sqrt(2) * im * (cos(r[1])*exp(-2*(r[2]^2+r[3]^2)/waists[i]^2) .+
        polSign[i] * sin(r[2])*exp(-2*(r[1]^2+r[3]^2)/waists[i]^2));
        
        fieldTerms[i][3] = cos(r[3])*exp(-2*(r[1]^2+r[2]^2)/waists[i]^2) .+ polSign[i] * sin(r[1])*exp(-2*(r[2]^2+r[3]^2)/waists[i]^2) .+
        im * (polSign[i] * sin(r[3])*exp(-2*(r[1]^2+r[2]^2)/waists[i]^2) - cos(r[2])*exp(-2*(r[1]^2+r[3]^2)/waists[i]^2));
        =#
        if polType[i]=="Slower"#polarization vs position for a given laser depends on whether it a "3D\sig\sig", "2D\sig\sig", slower, etc.
            fieldTerms[i][1] = 1/sqrt(2) .* (cos(r[3]) + im * sin(r[3])) .* exp(-(r[1]^2+r[2]^2)/waists[i]^2);
            fieldTerms[i][2] = 0;
            fieldTerms[i][3] = -1/sqrt(2) .* (cos(r[3]) + im * sin(r[3])) .* exp(-(r[1]^2+r[2]^2)/waists[i]^2);
        elseif polType[i]=="Push"
            fieldTerms[i][1] = 1/sqrt(2) .* (cos(r[3]) - im * sin(r[3])) .* exp(-(r[1]^2+r[2]^2)/waists[i]^2);
            fieldTerms[i][2] = 0;
            fieldTerms[i][3] = -1/sqrt(2) .* (cos(r[3]) - im * sin(r[3])) .* exp(-(r[1]^2+r[2]^2)/waists[i]^2);
        else #3D MOT configuration
            xTerm = -polSign[i] .* 1.0/sqrt(2.0) .* exp(-(r[1]^2+r[2]^2)/waists[i]^2) .* 
            (cos(r[3]+rp[i,2]) - cos(r[3]+rp[i,1]) - im * sin(r[3]+rp[i,1]) - im * sin(r[3]+rp[i,2])) .+
            1.0/sqrt(2.0) .* exp(-(r[1]^2+r[3]^2)/waists[i]^2) .* 
            (-im * cos(r[2]+rp[i,4]) -im * cos(r[2]+rp[i,3]) + sin(r[2]+rp[i,3]) - sin(r[2]+rp[i,4])); 

            yTerm = polSign[i] .* 1.0/sqrt(2.0) .* exp(-(r[2]^2+r[3]^2)/waists[i]^2) .* 
            (cos(r[1]+rp[i,6]) - cos(r[1]+rp[i,5]) - im * sin(r[1]+rp[i,5]) - im * sin(r[1]+rp[i,6])) .+
            1.0/sqrt(2.0) .* exp(-(r[2]^2+r[1]^2)/waists[i]^2) .* 
            (-im * cos(r[3]+rp[i,2]) -im * cos(r[3]+rp[i,1]) + sin(r[3]+rp[i,1]) - sin(r[3]+rp[i,2])); 

            zTerm = polSign[i] .* 1.0/sqrt(2.0) .* exp(-(r[3]^2+r[1]^2)/waists[i]^2) .* 
            (cos(r[2]+rp[i,4]) - cos(r[2]+rp[i,3]) - im * sin(r[2]+rp[i,3]) - im * sin(r[2]+rp[i,4])) .+
            1.0/sqrt(2.0) .* exp(-(r[3]^2+r[2]^2)/waists[i]^2) .* 
            (-im * cos(r[1]+rp[i,6]) -im * cos(r[1]+rp[i,5]) + sin(r[1]+rp[i,5]) - sin(r[1]+rp[i,6])); 
            
            fieldTerms[i][3] = conj(-1.0/sqrt(2.0) .* xTerm + im * 1.0 /sqrt(2.0) .* yTerm);
            fieldTerms[i][2] = conj(zTerm);
            fieldTerms[i][1] = conj(1.0/sqrt(2.0) .* xTerm + im * 1.0 /sqrt(2.0) .* yTerm);
        end
    end#end polSign (e.g. end iteration through lasers)
end

    
function makeBFieldTerms!(bFieldTerms,r::Array{Float64,1},bFieldSetting::String,theta::Float64)
    #expresses B field at position r in the \sigma^+/-, pi basis
    # r=zeros(1,3);

    if bFieldSetting=="TwoD" #this is for 2D quadrupole B-field gradient
        bFieldTerms[1] = 1 / sqrt(2) * (r[1] - im * r[2]);
        bFieldTerms[2] = -0 * (r[3]);
        bFieldTerms[3] = 1 / sqrt(2) * (-r[1] - im * r[2]);
    elseif bFieldSetting=="ThreeD" #this is for 3D quadrupole B-field gradient
        bFieldTerms[1] = 1 / sqrt(2) * (r[1] + im * r[2]);
        bFieldTerms[2] = -2 * (r[3]); #true quadrupole field
        bFieldTerms[3] = 1 / sqrt(2) * (-r[1] + im * r[2]);
    else#if Static, primarily for white light slowing
        
        #theta degrees along xy-axis (w.r.t. x-axis which is the slowing laser polarization) for dark state remixing during WLS
        bFieldTerms[1] = 1/sqrt(2) * (cos(theta*pi/180) + im * sin(theta*pi/180));
        bFieldTerms[2] = 0;
        bFieldTerms[3] = 1/sqrt(2) * (-cos(theta*pi/180) + im * sin(theta*pi/180));
        
        #=
        #special case of 45 degrees which Thomas used in all of his WLS code.
        bFieldTerms[1] = (im+1)/2;
        bFieldTerms[2] = 0;
        bFieldTerms[3] = (-1+im)/2;
        =#

        #=
        #special case of parallel to x-axis (stationary dark state), i.e. theta = 0 degrees
        bFieldTerms[1] = 1/sqrt(2);
        bFieldTerms[2] = 0;
        bFieldTerms[3] = -1/sqrt(2);
        =#
        
        #=
        #special case of parallel to y-axis (stationary dark state), i.e. theta = 90 degrees
        bFieldTerms[1] = im/sqrt(2);
        bFieldTerms[2] = 0;
        bFieldTerms[3] = im/sqrt(2);
        =#
        
        #=
        #special case of 54.7 degrees wrt x-axis (magic angle for maximal dark state remixing)
        bFieldTerms[1] = 1/sqrt(2) * (cos(54.7*pi/180) + im * sin(54.7*pi/180));
        bFieldTerms[2] = 0;
        bFieldTerms[3] = 1/sqrt(2) * (-cos(54.7*pi/180) + im * sin(54.7*pi/180));
        =#

        #=
        #static field along z-axis (only used for test of lambda cooling in 3-D)
        bFieldTerms[1] = 0;
        bFieldTerms[2] = 1;
        bFieldTerms[3] = 0;
        =#
    end

end

function densityMatrixChangeTerms!(du, u, p, t)
    #The meat of the program.  Here's where the density matrix is actually evolved.

    # user inputs (these vary, things like initial position, velocity, laser params, etc.).  These are determined by the user-chosen parameters in the main program
    rInit = p[1]::Vector{Float64};
    v = p[2]::Vector{Float64};
    rp = p[3]::Matrix{Float64};
    lasers = p[4]::Lasers{Vector{Float64}, Vector{Int64}, Vector{String}, Float64};
    laserEnergy = lasers.laserEnergy;
    s0 = lasers.s0;
    polSign = lasers.polSign;
    polType = lasers.polType;
    sidebandFreqs = lasers.sidebandFreqs;
    sidebandAmps = lasers.sidebandAmps;
    sidebandFreqs2 = lasers.sidebandFreqs2;
    sidebandAmps2 = lasers.sidebandAmps2;
    waists = lasers.waists;
    staticBFieldAngle = lasers.staticBFieldAngle;
    bToHamConvert = p[5]::Float64;

    # coupling matrices passed by user.  
    coupleMat1 = p[6]::Matrix{Float64};
    coupleMat2 = p[7]::Matrix{Float64};
    coupleMat3 = p[8]::Matrix{Float64};
    bCoupleMat1 = p[9]::Matrix{Float64};
    bCoupleMat2 = p[10]::Matrix{Float64};
    bCoupleMat3 = p[11]::Matrix{Float64};
    # coupling matrices used in decay calc
    coupleMatEff1 = p[12]::Matrix{ComplexF64};
    coupleMatEff2 = p[13]::Matrix{ComplexF64};
    coupleMatEff3 = p[14]::Matrix{ComplexF64};

    # decay 'masks' used in calculating the decay term.  
    decayMaskAllButTopLeft = p[15]::Matrix{Float64};
    decayMaskForCalcTopLeft = p[16]::Matrix{Int64};

    # pre-cached r Array
    r = p[17]::Vector{Float64};
    # pre-cached matrices for atom light term.  
    fieldTerms = p[18]::Vector{Vector{ComplexF64}};
    atomLightTerm = p[19]::Matrix{ComplexF64};
    atomLightTerm = zeros(ComplexF64, size(coupleMat1,1), size(coupleMat1,2));
    
    # pre-cached matrices for b field term. 
    bFieldTerms = p[20]::Vector{ComplexF64};
    bFieldTermFull = p[21]::Matrix{ComplexF64} ;
    uProdBField = p[22]::Matrix{ComplexF64} ;
    bFieldProdU = p[23]::Matrix{ComplexF64} ;

    # pre-cached matrices for decay term
    decayFull = p[24]::Matrix{ComplexF64};
    pOnlyExcitedStates = p[25]::Matrix{ComplexF64};
    pTopLeft1PreMult = p[26]::Matrix{ComplexF64};
    pTopLeft2PreMult = p[27]::Matrix{ComplexF64};
    pTopLeft3PreMult = p[28]::Matrix{ComplexF64};
    pTopLeft1 = p[29]::Matrix{ComplexF64};
    pTopLeft2 = p[30]::Matrix{ComplexF64};
    pTopLeft3 = p[31]::Matrix{ComplexF64};


    #1) evolve position
    propR!(r, rInit, v, t);
    #2) Calculate field terms at new position
    makeFieldTerms!(fieldTerms, r, polSign,polType,waists, rp);
    #3)calculate -E dot D term (see Eq 21 of writeup)
    for i=1:length(s0)
        atomLightTerm .= atomLightTerm.+ sqrt(s0[i] / 8) .* -exp(1im * laserEnergy[i] * t + 1im * sidebandAmps[i]*sin(sidebandFreqs[i] * t) + 1im * sidebandAmps2[i]*sin(sidebandFreqs2[i] * t)) .* ((fieldTerms[i][1] .* coupleMat1) .+
         (fieldTerms[i][2] .* coupleMat2) .+ (fieldTerms[i][3] .* coupleMat3));
    end
    atomLightTerm .= atomLightTerm .+ atomLightTerm';#needed here because, the way coupleMat is defined, 'atomLightTerm' up til now only has the top right half of the hermitian coupling matrix

    #4) calculate -mu dot B term (see Eq 32 of writeup)
    #makeBFieldTerms!(bFieldTerms, r, bFieldSetting, staticBFieldAngle[1]);
    makeBFieldTerms!(bFieldTerms, r, bFieldSetting, staticBFieldAngle);
    bFieldTermFull .= bToHamConvert .* (bFieldTerms[1] .* bCoupleMat1 .+ bFieldTerms[2] .* bCoupleMat2 .+ bFieldTerms[3] .* bCoupleMat3).+atomLightTerm; #'bTermFull' also sums the -mu dot B term with the calculated -D dot E term
    #5) take commutator of [H,u] where u is density matrix and H = -D dot E + -mu dot B
    mul!(uProdBField, u, bFieldTermFull);#uH
    mul!(bFieldProdU, bFieldTermFull, u);#Hu

    #6) Take decay (aka 'coupling to reservoir') into account (Eq 46 of writeup)
    pOnlyExcitedStates .= u .* decayMaskForCalcTopLeft;

    #6A) these next 6 lines calculate the last term in eq 46 of writeup
    coupleMatEff1 = coupleMat1 #.* exp.(-1im*t.*stateEnergyMatrix);
    mul!(pTopLeft1PreMult, coupleMatEff1, pOnlyExcitedStates);
    mul!(pTopLeft1, pTopLeft1PreMult, coupleMatEff1')

    coupleMatEff2 = coupleMat2 #.* exp.(-1im*t.*stateEnergyMatrix);
    mul!(pTopLeft2PreMult, coupleMatEff2, pOnlyExcitedStates);
    mul!(pTopLeft2, pTopLeft2PreMult, coupleMatEff2')

    coupleMatEff3 = coupleMat3 #.* exp.(-1im*t.*stateEnergyMatrix);
    mul!(pTopLeft3PreMult, coupleMatEff3, pOnlyExcitedStates);
    mul!(pTopLeft3, pTopLeft3PreMult, coupleMatEff3')

    decayFull .= (u .* decayMaskAllButTopLeft) .+ pTopLeft1 .+ pTopLeft2 .+ pTopLeft3;#u.*decayMask term represents 1st and 2nd term of eq 46 in writeup

    du .= 1im .* (uProdBField .- bFieldProdU) .+ decayFull;#finally, add the 'Liouville' term and the decay term (Eq 1 of writeup) to step the density matrix
end

#everything else here is used to calculate a force given a density matrix, see section 1.6 of writeup

function makeDFieldTerms!(dFieldTerms, r::Array{Float64,1}, polSign::Array{Int64,1}, polType::Array{String,1}, waists::Array{Float64,1},rp::Matrix{Float64})
    # dfieldterms (dE/dr) will be 3x3 matrix, first element is xyz second is sig+ pi sig-
    # fieldTerms = zeros(ComplexF64,3,1);
    # pre conjugated, just like 'makeFieldTerms'
    for i=1:length(polSign) #iterates through all lasers
        #=
        #original way in code with no random phase
        dFieldTerms[i][1,1] = polSign[i] * cos(r[1])*exp(-2*(r[2]^2+r[3]^2)/waists[i]^2);
        dFieldTerms[i][1,2] = -sqrt(2) * im * sin(r[1])*exp(-2*(r[2]^2+r[3]^2)/waists[i]^2);
        dFieldTerms[i][1,3] = polSign[i] * cos(r[1])*exp(-2*(r[2]^2+r[3]^2)/waists[i]^2) ;

        dFieldTerms[i][2,1] = -im * sin(r[2])*exp(-2*(r[1]^2+r[3]^2)/waists[i]^2);
        dFieldTerms[i][2,2] = sqrt(2) * im * (polSign[i] * cos(r[2])*exp(-2*(r[1]^2+r[3]^2)/waists[i]^2));
        dFieldTerms[i][2,3] = im * (sin(r[2])*exp(-2*(r[1]^2+r[3]^2)/waists[i]^2));

        dFieldTerms[i][3,1] = (-sin(r[3])*exp(-2*(r[1]^2+r[2]^2)/waists[i]^2) - im * (polSign[i] * cos(r[3])*exp(-2*(r[1]^2+r[2]^2)/waists[i]^2)));
        dFieldTerms[i][3,2] = 0;
        dFieldTerms[i][3,3] = (-sin(r[3])*exp(-2*(r[1]^2+r[2]^2)/waists[i]^2) + im * (polSign[i] * cos(r[3])*exp(-2*(r[1]^2+r[2]^2)/waists[i]^2)));
        =#

        if polType[i]=="Slower" #polarization vs position for a given laser depends on whether it a "3D\sig\sig", "2D\sig\sig", slower, etc.
            dFieldTerms[i][1,1] = 0;
            dFieldTerms[i][1,2] = 0;
            dFieldTerms[i][1,3] = 0;

            dFieldTerms[i][2,1] = 0;
            dFieldTerms[i][2,2] = 0;
            dFieldTerms[i][2,3] = 0;

            dFieldTerms[i][3,1] = 1/sqrt(2) * (-sin(r[3]) + im * cos(r[3])) .* exp(-(r[1]^2+r[2]^2)/waists[i]^2);
            dFieldTerms[i][3,2] = 0;
            dFieldTerms[i][3,3] = -1/sqrt(2) * (-sin(r[3]) + im * cos(r[3])) .* exp(-(r[1]^2+r[2]^2)/waists[i]^2);

        elseif polType[i]=="Push"
            dFieldTerms[i][1,1] = 0;
            dFieldTerms[i][1,2] = 0;
            dFieldTerms[i][1,3] = 0;

            dFieldTerms[i][2,1] = 0;
            dFieldTerms[i][2,2] = 0;
            dFieldTerms[i][2,3] = 0;

            dFieldTerms[i][3,1] = 1/sqrt(2) * (-sin(r[3]) - im * cos(r[3])) .* exp(-(r[1]^2+r[2]^2)/waists[i]^2);
            dFieldTerms[i][3,2] = 0;
            dFieldTerms[i][3,3] = -1/sqrt(2) * (-sin(r[3]) - im * cos(r[3])) .* exp(-(r[1]^2+r[2]^2)/waists[i]^2);
        else #3D MOT configuration
            dFieldTerms[i][1,3] = conj(1.0/2.0 .* exp(-(r[2]^2+r[3]^2)/waists[i]^2) .* im .* polSign[i] .*
            (-im * cos(r[1]+rp[i,5]) - im * cos(r[1]+rp[i,6]) + sin(r[1]+rp[i,5]) - sin(r[1]+rp[i,6])));

            dFieldTerms[i][1,2] = conj(1.0/sqrt(2.0) .* exp(-(r[2]^2+r[3]^2)/waists[i]^2) .* 
            (cos(r[1]+rp[i,5]) - cos(r[1]+rp[i,6]) + im * sin(r[1]+rp[i,5]) + im*sin(r[1]+rp[i,6])));

            dFieldTerms[i][1,1] = conj(1.0/2.0 .* exp(-(r[2]^2+r[3]^2)/waists[i]^2) .* im .* polSign[i] .*
            (-im * cos(r[1]+rp[i,5]) - im * cos(r[1]+rp[i,6]) + sin(r[1]+rp[i,5]) - sin(r[1]+rp[i,6])));


            dFieldTerms[i][2,3] = conj(1.0/2.0 .* exp(-(r[3]^2+r[1]^2)/waists[i]^2) .*
            (-cos(r[2]+rp[i,3]) + cos(r[2]+rp[i,4]) - im * sin(r[2]+rp[i,3]) - im * sin(r[2]+rp[i,4])));

            dFieldTerms[i][2,2] = conj(1.0/sqrt(2.0) .* exp(-(r[3]^2+r[1]^2)/waists[i]^2) .* polSign[i] .*
            (-im*cos(r[2]+rp[i,3]) - im*cos(r[2]+rp[i,4]) + sin(r[2]+rp[i,3]) - sin(r[2]+rp[i,4])));

            dFieldTerms[i][2,1] = conj(-1.0/2.0 .* exp(-(r[3]^2+r[1]^2)/waists[i]^2) .*
            (-cos(r[2]+rp[i,3]) + cos(r[2]+rp[i,4]) - im * sin(r[2]+rp[i,3]) - im * sin(r[2]+rp[i,4])));

            dFieldTerms[i][3,3] = conj(1.0/2.0 .* im .* exp(-(r[1]^2+r[2]^2)/waists[i]^2) .* (cos(r[3]+rp[i,1]) - cos(r[3]+rp[i,2]) + im*sin(r[3]+rp[i,1]) + im*sin(r[3]+rp[i,2])) .-
            1.0/2.0 .* polSign[i] .* exp(-(r[1]^2+r[2]^2)/waists[i]^2) .* (im*cos(r[3]+rp[i,1]) + im*cos(r[3]+rp[i,2]) - sin(r[3]+rp[i,1]) + sin(r[3]+rp[i,2])));

            dFieldTerms[i][3,2] = 0;

            dFieldTerms[i][3,1] = conj(1.0/2.0 .* im .* exp(-(r[1]^2+r[2]^2)/waists[i]^2) .* (cos(r[3]+rp[i,1]) - cos(r[3]+rp[i,2]) + im*sin(r[3]+rp[i,1]) + im*sin(r[3]+rp[i,2])) .+
            1.0/2.0 .* polSign[i] .* exp(-(r[1]^2+r[2]^2)/waists[i]^2) .* (im*cos(r[3]+rp[i,1]) + im*cos(r[3]+rp[i,2]) - sin(r[3]+rp[i,1]) + sin(r[3]+rp[i,2])));
        end
    end#end polSign

    # return dfieldTerms
end

function forceCalc!(force, dFieldTerms::Vector{Matrix{ComplexF64}}, rho::Matrix{ComplexF64}, 
    lasers::Lasers{Vector{Float64}, Vector{Int64}, Vector{String}, Float64}, couplingMatrices::Vector{Matrix}, t::Float64)
    #calculates force given position and lasers (used to calculate dFieldTerms) and density matrix \rho.  
    #Both r(t) and \rho(t) are recorded vs time by the OBE solver, so this runs afterwords to calculate what forces the particle experienced over the trajectory
    s0 = lasers.s0;
    sidebandFreqs = lasers.sidebandFreqs;
    sidebandFreqs2 = lasers.sidebandFreqs2;
    sidebandAmps = lasers.sidebandAmps;
    sidebandAmps2 = lasers.sidebandAmps2;
    laserEnergy = lasers.laserEnergy;

    #force pre-factor is calculated for each laser [i].  Has the rotating-frame frequency exponent + phase modulation term + intensity term \sqrt(s0/8).  Hyperfine energies are subtracted later
    forcePrefactor = zeros(ComplexF64,1,length(laserEnergy));
    for i=1:length(laserEnergy)
        forcePrefactor[i] = sqrt(s0[i] / 8) * exp(1im * laserEnergy[i] * t + 1im * sidebandAmps[i]*sin(sidebandFreqs[i] * t) + 1im * sidebandAmps2[i]*sin(sidebandFreqs2[i] * t));
    end

    #calculate x force.  Implements Eq 48 of main writeup
    dRhoDPosCalcMatrix = zeros(ComplexF64, size(rho,1), size(rho,2));
    dRhoDPosTimesDensityMatContainer = zeros(ComplexF64, size(rho,1), size(rho,2));
    for i=1:length(laserEnergy)
        dRhoDPosCalcMatrix .= dRhoDPosCalcMatrix .+ forcePrefactor[i] * (dFieldTerms[i][1,1] * couplingMatrices[1] .+ dFieldTerms[i][1,2] * couplingMatrices[2] .+ dFieldTerms[i][1,3] * couplingMatrices[3]);
    end
    dRhoDPosCalcMatrix .= dRhoDPosCalcMatrix .+ dRhoDPosCalcMatrix';
    mul!(dRhoDPosTimesDensityMatContainer, rho, dRhoDPosCalcMatrix);#multiplies dp_{x}/dt by density matrix \rho
    force[1] = real(tr(dRhoDPosTimesDensityMatContainer));#takes trace of \rho*dp_{x}/dt to determine average force over enemble (Eq 49 of writeup)

    #similarly, calculate y and z force
    dRhoDPosCalcMatrix = zeros(ComplexF64, size(rho,1), size(rho,2));
    for i=1:length(laserEnergy)
        dRhoDPosCalcMatrix .= dRhoDPosCalcMatrix .+ forcePrefactor[i] * (dFieldTerms[i][2,1] * couplingMatrices[1] .+ dFieldTerms[i][2,2] * couplingMatrices[2] .+ dFieldTerms[i][2,3] * couplingMatrices[3]);
    end
    dRhoDPosCalcMatrix .= dRhoDPosCalcMatrix .+ dRhoDPosCalcMatrix';
    mul!(dRhoDPosTimesDensityMatContainer, rho, dRhoDPosCalcMatrix);
    force[2] = real(tr(dRhoDPosTimesDensityMatContainer));

    dRhoDPosCalcMatrix = zeros(ComplexF64, size(rho,1), size(rho,2));
    for i=1:length(laserEnergy)
        dRhoDPosCalcMatrix .= dRhoDPosCalcMatrix .+ forcePrefactor[i] * (dFieldTerms[i][3,1] * couplingMatrices[1] .+ dFieldTerms[i][3,2] * couplingMatrices[2] .+ dFieldTerms[i][3,3] * couplingMatrices[3]);
    end
    dRhoDPosCalcMatrix .= dRhoDPosCalcMatrix .+ dRhoDPosCalcMatrix';
    mul!(dRhoDPosTimesDensityMatContainer, rho, dRhoDPosCalcMatrix);
    force[3] = real(tr(dRhoDPosTimesDensityMatContainer));
end#forceCalc!

function makeForceVsTime!(forceVsTime, times, rhos, lasers, couplingMatrices, rInit, v, rp)
    #given a set of times, an initial position and velocity, the lasers used, and \rho(t), calculate force vs t
    polSign = lasers.polSign;
    polType = lasers.polType;
    waists = lasers.waists;
    #initialize some stuff
    dFieldContainer = Array{Array{ComplexF64, 2},1}(undef,length(polSign));
    for i=1:length(polSign)
        dFieldContainer[i]=zeros(ComplexF64,3,3)
    end#end polSign
    forceCalcContainer = zeros(ComplexF64, 3, 1);
    r = Array{Float64,1}(undef, 3)
    #iterate through time, propegating r in the same way done in the OBEs.  Then determine force experienced given r(t), \rho(t), and the lasers used
    for i = 1:length(times)
        propR!(r, rInit, v, times[i]);
        makeDFieldTerms!(dFieldContainer, r, polSign,polType, waists,rp)
        forceCalc!(forceCalcContainer, dFieldContainer, rhos[i], lasers, couplingMatrices,times[i]);
        forceVsTime[i,:] = forceCalcContainer;
    end#end times

end#makeForceVsTime!
