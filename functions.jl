function solveBGP(p)


  eqvInit = [1.9167531502801942e+00,
             5.5035518427253116e-01,
             6.2816712824469115e-02,
            7.1398401530217703e-01,
            8.3481281988279321e-01,
            1.7306237409854832e+00]

  eqvInit = [ 1.9165753642822532,
            0.5503998860506523,
             0.06279518700611891,
             0.7138607955216772,
             0.8346597073502943,
             1.7306541656053274]


 function objFnc(x)
   out, _ = eqfunc(x,p);
   return out
 end

 res = nlsolve(objFnc, eqvInit, method = :trust_region,inplace = false);  # to-do: make inplace true version

 if !res.f_converged
     print("👎")
   eq = EqObj();  # return empty EqObj, type stability
 else
     print("👍")
     err,eq = eqfunc(res.zero,p);
 end


 return eq,res


end

######################################################

function eqfunc(eqv,p)
  # reject out of bounds
  if (any(eqv.<0.0))
    eqnd = -1000.0*ones(size(eqv)[1])

  else
    eq = EqObj();
    # load in guesses for solution
    eq.ws      = eqv[1]
    eq.cactiv  = eqv[2:3]
    eq.eyq     = eqv[4:5]
    eq.qbar    = eqv[6]


    # flow of free products
    eq.cactivtot = sum(eq.cactiv)
    eq.cinac     = 1.0 - eq.cactivtot


    # find updated eqvars
    innovation!(eq,p)
    qualityDist!(eq,p)
    qbarActive!(eq,p)
    calcey!(eq,p)
    labordem!(eq,p)


    # equilibrium differences
    eqnd = -1000.0*ones(size(eqv)[1])
    eqnd[1]   = eq.skilled_lab - p.Ls                # skilled labor market clearing
    eqnd[2:3] = eq.cactiv      - eqv[2:3]       # product shares
    eqnd[4:5] = eq.eyq         - eqv[4:5]       # innovated product value
    eqnd[6]   = eq.qbarAct     - (p.ε/(p.ε-1.0))^(p.ε - 1)


  end
  return eqnd,eq

end


###############################################################################################
###############################################################################################

function innovation!(eq,p)

  # innovation rates
  eq.x    = p.𝚯.* (((1.0 - p.γ)/eq.ws)*eq.eyq).^((1.0 - p.γ)/p.γ)
  eq.xout = p.θE.*(((1.0 - p.γE)/eq.ws).*sum(p.𝛂.*eq.eyq)).^((1.0 - p.γE)/p.γE)

  # option value of R&D
  eq.optval = p.γ*p.𝚯*((1.0 - p.γ)/eq.ws)^((1.0 - p.γ)/p.γ).*eq.eyq.^(1.0/p.γ)

  # aggregate creative destruction rates
  eq.taus = eq.cactiv.*eq.x .+ p.𝛂*eq.xout
  eq.tau  = sum(eq.taus)

  # minimum qhat
  eq.qmin = (max.(0,eq.ws*p.ϕ .- eq.optval)/p.Π).^(1.0/(p.ε-1.0))

end


###############################################################################################
###############################################################################################

function qualityDist!(eq,p)

    # -1 - Parameters for delayed differential equation
    alphaDelay = 1.0/(1.0 + (p.λ-1.0)*(1.0 - p.ω));
    betaDelay  = (p.λ - 1.0)*p.ω*eq.qbar*alphaDelay;

    # 0 - Growth rate consistent with guessed qbar
    eq.g = eq.tau*(p.λ - 1.0);

    # 1 - Overall dist.function
    # to-do: for now it is written for ω=1. Make it general
    function funcFAll!(du,u,h,p,t)  # to-do: make it inline
        du[1] = (eq.tau/(eq.g*t))*u[1] - (eq.tau/(eq.g*t))*h(p,t-betaDelay)[1] # to-do: tau/g= 1/(λ-1), we can simply use that
    end

    #prob = DDEProblem(funcFAll!,u0,h,tspan, constant_lags = lags)
    tspan = (1.0e-09,10.0)
    lags = [betaDelay]
    h(p,t)= [0.0];
    probDelay  = DDEProblem(funcFAll!,[1.0e-16],h,tspan, constant_lags = lags)
    eq.solFAll = solve(probDelay,MethodOfSteps(Tsit5()))

    eq.FAllend = eq.solFAll.u[end][1,1]  # needed for normalization


    # 2 - Gross dist.
    # different than matlab verison, we convert the problem to IVP by simply reversing time span
    function funcFH!(du,u,pp,t)
        delayT 	= min.(max.(0.0000001,t - betaDelay),10.0)
        du[1] 	= ((eq.tau + p.ν)/(eq.g*t))*u[1] - (eq.taus[2]/(eq.g*t))*eq.solFAll(delayT)[1]/eq.FAllend
    end

    probFH = ODEProblem(funcFH!,[eq.taus[2]/(eq.tau + p.ν)],(tspan[2],tspan[1]))
    eq.solFH = solve(probFH, Tsit5(),reltol=1e-8,abstol=1e-8)

    eq.PhiHG = eq.taus[2]/(eq.tau + p.ν);
    eq.PhiLG = 1.0 - eq.PhiHG;

    # 3 - Active product line dist.
    # Again same trick, solve IVP
    function funcFRest!(du,u,pp,t)
      FHt   = eq.solFH(t)[1]
      du[1] = ((eq.tau + p.ν + p.ψ)/(eq.g*t))*u[1] - (p.ψ/(eq.g.*t))*FHt
      du[2] = ((eq.tau +       p.ψ)/(eq.g*t))*u[2] - (p.ψ/(eq.g*t))*(eq.solFAll(t)[1]/eq.FAllend - FHt) - (p.ν./(eq.g*t))*u[1]
    end

    boundHRest = (p.ψ/(eq.tau + p.ν + p.ψ))*eq.PhiHG
    boundLRest = (p.ψ/(eq.tau + p.ψ))*eq.PhiLG + (p.ν/(eq.tau + p.ψ))*boundHRest

    probFRest   = ODEProblem(funcFRest!,[boundHRest,boundLRest],(tspan[2],tspan[1]))
    eq.solFRest = solve(probFRest, Tsit5(),reltol=1e-8,abstol=1e-8)

    # 4 - Update PhiL and PhiH
    FRestcut    = eq.solFRest(max.(1.0e-09,eq.qmin))
    FHcut       = eq.solFH(max.(1.0e-09,eq.qmin))
    FAllcut     = eq.solFAll(max.(1.0e-09,eq.qmin))./eq.FAllend

    eq.PhiHExo   = eq.solFRest(10.0)[1]
    eq.PhiLExo   = eq.solFRest(10.0)[2]

    eq.PhiHObso  = FHcut[1,2] - FRestcut[1,2]
    eq.PhiLObso  = (FAllcut[1,1] - FHcut[1,1]) - FRestcut[2,1]


    eq.cactiv[1] = (eq.PhiLG - eq.PhiLExo) - eq.PhiLObso
    eq.cactiv[2] = (eq.PhiHG - eq.PhiHExo) - eq.PhiHObso
    eq.cactivtot = sum(eq.cactiv);
    eq.cinac     = 1.0 - eq.cactivtot;


end

########################################

function qbarActive!(eq,p)


  #function funQbarAct(x; m1=eq.solFAll, m2 = eq.solFH, m3 = eq.solFRest,eps = p.ε,qmin = eq.qmin)
  function funQbarAct(x)
          outFAll  = eq.solFAll(x,Val{1})/eq.solFAll.u[end][1,1]  # Val{1} enables du return
          outFH    = eq.solFH(x,Val{1})
          outRest  = eq.solFRest(x,Val{1})

          outL = (outFAll[1] - outFH[1] - outRest[2])*(x^(p.ε-1))*(x>=eq.qmin[1])
          outH = (outFH[1]   - outRest[1])           *(x^(p.ε-1))*(x>=eq.qmin[2])
          out  = hcat(outL, outH)
  end
    eq.qbarAct  = sum(quadgk(funQbarAct,1.0e-08,10.0)[1])
end

#########################################

###################################################


function calcey!(eq,p)

  eq.r = p.ρ + p.σ*eq.g;


  function eyfunc(q; solDiff = eq.solFAll, qbar = eq.qbar, λ = p.λ, ω = p.ω, ν = p.ν)
    out   = solDiff(q,Val{1})/solDiff.u[end][1,1]
    qPlus = q + (λ - 1)*(ω*qbar + (1 - ω)*q)
    eyL = zfunc(qPlus,0.0,1)[1]*out
    eyH = (zfunc(qPlus,ν,2)[1] + zfunc(qPlus,0.0,1)[1] - zfunc(qPlus,ν,1)[1])*out
    eyout = vec(hcat(eyL, eyH))
            # LOW TYPE                 # HIGH TYPE
  end

  ###################################################

  function zfunc(q,x,s; tau = eq.tau, r = eq.r, ψ = p.ψ, ε = p.ε, g = eq.g, optval = eq.optval, ws = eq.ws, ϕ = p.ϕ, qmin = eq.qmin)

    psit = r + tau + ψ + x

    kappa1 = psit + (ε - 1.0)*g
    kappa2 = psit

    coeff1 = p.Π
    coeff2 = optval[s] - ws*ϕ;

    # this min means that when q <= qmin, qrat = 1.0, so the value contribution is zero
    qrat = min(1.0,qmin[s]/q);

    zout = coeff1/kappa1*q^(ε-1.0)*(1.0-qrat^(kappa1/g)) + coeff2/kappa2*(1.0-qrat^(kappa2/g))

  end

  eq.eyq = quadgk(eyfunc,1.0e-08,10.0)[1]

end
###################################################



function labordem!(eq,p)

  # managerial fixed cost
  eq.cfix = eq.cactivtot*p.ϕns

  # production labor
  eq.wu = (p.ε - 1.0)/p.ε   # unskilled wage over Q.

  # incumbent R&D labor
  eq.cx = (eq.x./(p.𝚯ns.^p.γ)).^(1.0/(1.0-p.γ))
  eq.crnd = sum(eq.cactiv.*eq.cx)

  # entrant R&D labor
  eq.cout = (eq.xout/(p.θEns^p.γE))^(1.0/(1.0-p.γE));

  # total skilled labor
  eq.skilled_lab = eq.crnd + eq.cfix + eq.cout
end
