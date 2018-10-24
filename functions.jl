function solveBGP(p)

  eqvInit = [ 2.91739277816310,
              0.350560526647831,
              0.162895991993503,
              0.513246861682734,
              0.234732865047170,
              1.632272083084665
              ]

  eq = EqObj();

  # gausslegendre for integration
  eq.qMinAll = 1.0e-8
  eq.qMaxAll = 10.0
  eq.node,eq.weight = gausslegendre(1000)
  eq.node = ((eq.qMaxAll - eq.qMinAll)/2)*eq.node .+ (eq.qMaxAll + eq.qMinAll)/2    # adjust the limits



  @inline function objFnc(eqnd,x)
     eqfunc!(eqnd,x,eq,p);
  end

  res = nlsolve(objFnc, eqvInit, method = :trust_region,inplace = true)

  if !res.f_converged
       print("👎")
       eq = EqObj();  # return empty EqObj, type stability
  else
       print("👍")
  end


  return eq,res


end

######################################################

function eqfunc!(eqnd,eqv,eq,p)
  # reject out of bounds

  if (any(eqv.<0.0))
    eqnd = -1000.0*ones(size(eqv)[1])

  else


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

    eqnd[1]   = eq.skilled_lab - p.Lˢ                # skilled labor market clearing
    eqnd[2:3] = eq.cactiv      - eqv[2:3]       # product shares
    eqnd[4:5] = eq.eyq         - eqv[4:5]       # innovated product value
    eqnd[6]   = eq.qbarAct     - (p.ε/(p.ε-1.0))^(p.ε - 1)


  end

end


###############################################################################################
###############################################################################################

function innovation!(eq,p)

  # innovation rates
  eq.x    = p.𝚯.* (((1.0 - p.γ)/eq.ws)*eq.eyq).^((1.0 - p.γ)/p.γ)
  eq.xout = p.θᴱ.*(((1.0 - p.γᴱ)/eq.ws).*sum(p.𝛂.*eq.eyq)).^((1.0 - p.γᴱ)/p.γᴱ)

  # option value of R&D
  eq.optval = p.γ*p.𝚯*((1.0 - p.γ)/eq.ws)^((1.0 - p.γ)/p.γ).*eq.eyq.^(1.0/p.γ)

  # aggregate creative destruction rates
  eq.τs = eq.cactiv.*eq.x .+ p.𝛂*eq.xout
  eq.τ  = sum(eq.τs)

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
    eq.g = eq.τ*(p.λ - 1.0);

    # 1 - Overall dist.function
    # to-do: for now it is written for ω=1. Make it general
    function funcFAll!(du,u,h,p,t)
        du[1] = (eq.τ/(eq.g*t))*u[1] - (eq.τ/(eq.g*t))*h(p,t-betaDelay)[1] # to-do: τ/g= 1/(λ-1), we can simply use that
    end

    #prob = DDEProblem(funcFAll!,u0,h,tspan, constant_lags = lags)
    tspan = (1.0e-09,10.0)
    lags = [betaDelay]
    h(p,t)= [0.0];
    probDelay  = DDEProblem(funcFAll!,[1.0e-16],h,tspan, constant_lags = lags)
    eq.solFAll = solve(probDelay,MethodOfSteps(Tsit5()),reltol=1e-8,abstol=1e-8) # to-do: we can directly get output at the nodes by using ``saveat``
    # option and pass it to integration functions. One thing to note: once saveat is used, between node results are linearly interpolated.

    eq.FAllend = eq.solFAll.u[end][1,1]  # needed for normalization

    eq.qAllDist  = vcat_nosplat(eq.solFAll(eq.node,Val{1}),1).*(1/eq.FAllend) # saving for later use


    # 2 - Gross dist.
    # different than matlab verison, we convert the problem to IVP by simply reversing time span

    eq.PhiHG = eq.τs[2]/(eq.τ + p.ν);
    eq.PhiLG = 1.0 - eq.PhiHG;

    # 3 - Active product line dist.
    # Again same trick, solve IVP
    function funcFRest!(du,u,pp,t)
      delayT 	= min.(max.(0.0000001,t - betaDelay),10.0)
      du[1] 	= ((eq.τ + p.ν)/(eq.g*t))*u[1] - (eq.τs[2]/(eq.g*t))*eq.solFAll(delayT)[1]/eq.FAllend
      du[2] = ((eq.τ + p.ν + p.ψ)/(eq.g*t))*u[2] - (p.ψ/(eq.g.*t))*u[1]
      du[3] = ((eq.τ +       p.ψ)/(eq.g*t))*u[3] - (p.ψ/(eq.g*t))*(eq.solFAll(t)[1]/eq.FAllend - u[1]) - (p.ν./(eq.g*t))*u[2]
    end

    boundHRest = (p.ψ/(eq.τ + p.ν + p.ψ))*eq.PhiHG
    boundLRest = (p.ψ/(eq.τ + p.ψ))*eq.PhiLG + (p.ν/(eq.τ + p.ψ))*boundHRest

    probFRest   = ODEProblem(funcFRest!,[eq.PhiHG,boundHRest,boundLRest],(tspan[2],tspan[1]))
    eq.solFRest = solve(probFRest, Tsit5(),reltol=1e-8,abstol=1e-8)

    # 4 - Update PhiL and PhiH
    FRestcut    = eq.solFRest(max.(1.0e-09,eq.qmin))
    FAllcut     = eq.solFAll(max.(1.0e-09,eq.qmin))./eq.FAllend

    eq.PhiHExo   = eq.solFRest(10.0)[2]
    eq.PhiLExo   = eq.solFRest(10.0)[3]

    eq.PhiHObso  = FRestcut[1,2] - FRestcut[2,2]
    eq.PhiLObso  = (FAllcut[1,1] - FRestcut[1,1]) - FRestcut[3,1]


    eq.cactiv[1] = (eq.PhiLG - eq.PhiLExo) - eq.PhiLObso
    eq.cactiv[2] = (eq.PhiHG - eq.PhiHExo) - eq.PhiHObso
    eq.cactivtot = sum(eq.cactiv);
    eq.cinac     = 1 - eq.cactivtot


end

########################################

function qbarActive!(eq,p)
  outL,outH = funQbarAct(eq,p)
  eq.qbarAct  = ((eq.qMaxAll - eq.qMinAll)/2)*(dot(outL,eq.weight) + dot(outH,eq.weight))
end
########

function funQbarAct(eq,p)


        outRest  = eq.solFRest(eq.node,Val{1})

        outL = (eq.qAllDist .- vcat_nosplat(outRest,1) .- vcat_nosplat(outRest,3)).*(eq.node.^(p.ε-1)).*(eq.node.>=eq.qmin[1])
        outH = (vcat_nosplat(outRest,1)                .- vcat_nosplat(outRest,2)).*(eq.node.^(p.ε-1)).*(eq.node.>=eq.qmin[2])
        return outL, outH
end


###################################################


function calcey!(eq,p)

  eq.r = p.ρ + p.σ*eq.g

  outL, outH = eyfunc(eq,p)

  eq.eyq = [((eq.qMaxAll - eq.qMinAll)/2)*dot(outL,eq.weight), ((eq.qMaxAll - eq.qMinAll)/2)*dot(outH,eq.weight)]


end
###################################################
function eyfunc(eq,p)

  qPlus = eq.node .+ (p.λ - 1).*(p.ω*eq.qbar .+ (1 - p.ω).*eq.node)

  eyL = zfunc(qPlus,0.0,1,eq,p).*eq.qAllDist     # LOW TYPE
  eyH = (zfunc(qPlus,p.ν,2,eq,p) .+ zfunc(qPlus,0.0,1,eq,p) .- zfunc(qPlus,p.ν,1,eq,p)).*eq.qAllDist      # HIGH TYPE

  return eyL, eyH

end

###################################################

function zfunc(q, x, s, eq, p)

  psit = eq.r + eq.τ + p.ψ + x

  kappa1 = psit + (p.ε - 1.0)*eq.g
  kappa2 = psit

  coeff1 = p.Π
  coeff2 = eq.optval[s] - eq.ws*p.ϕ;

  # this min means that when q <= qmin, qrat = 1.0, so the value contribution is zero
  qrat = min.(1.0,eq.qmin[s]./q);

  zout = coeff1/kappa1.*q.^(p.ε-1.0).*(1.0 .- qrat.^(kappa1/eq.g)) .+ coeff2/kappa2.*(1.0 .- qrat.^(kappa2/eq.g))

end

#####################


function labordem!(eq,p)

  # managerial fixed cost
  eq.cfix = eq.cactivtot*p.ϕⁿˢ

  # production labor
  eq.wu = (p.ε - 1.0)/p.ε   # unskilled wage over Q.

  # incumbent R&D labor
  eq.cx = (eq.x./(p.𝚯ⁿˢ.^p.γ)).^(1.0/(1.0-p.γ))
  eq.crnd = sum(eq.cactiv.*eq.cx)

  # entrant R&D labor
  eq.cout = (eq.xout/(p.θᴱⁿˢ^p.γᴱ))^(1.0/(1.0-p.γᴱ));

  # total skilled labor
  eq.skilled_lab = eq.crnd + eq.cfix + eq.cout
end

#########################################################
# Utility Functions

vcat_nosplat(y,l) =  eltype(y[l])[el[l] for el in y] # this is for converting weird output of diff solver
