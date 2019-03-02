
# SOCIAL PLANNER PROBLEM

function solveSocialPlanner(p)
  
  # first solve baseline 
  eqInit = [ 2.00, 0.55, 0.06,
             0.72, 0.83, 1.73]

  eq,res = solveBGP(p,eqInit) 
  
  # initial values for the solver and minimization
  socPolInit = [eq.qmin; eq.x]
  eqSocInit = [eq.qbar;eq.cactiv;eq.xout]
  
  cevMax = [calCev(eq.g, eq.cactivtot, eq.g, eq.cactivtot,p)]	
  @printf("SOCIAL PLANNER PROBLEM\n")
  @printf("-----------------------\n")
  @printf("Baseline Welfare is %3.4f (xˡ = %3.4f, xʰ = %3.4f, qminˡ = %3.4f, qminʰ = %3.4f) \n", cevMax[1],eq.x[1],eq.x[2],eq.qmin[1],eq.qmin[2])

  # optimization
  # ------------------------------------------------------------------------------------------------
  opt = Opt(:LN_NELDERMEAD, length(socPolInit))
  ftol_rel!(opt,1e-8)
  min_objective!(opt, (x, grad) -> socialPlannerFunc(x, p, eqSocInit, cevMax, eq.g, eq.cactivtot)[1])
  minf,minpol,ret = NLopt.optimize(opt, socPolInit)
  nEvals = opt.numevals
  #-------------------------------------------------------------------------------------------------
  
  if ret == :FTOL_REACHED 
  	@printf("Optimal policy found!\n\n")
    eqSoc = socialPlannerFunc(minpol,  p, eqSocInit, cevMax, eq.g, eq.cactivtot)[2]
	res = Dict(:cev=>-minf, :ret=>ret, :opt=>opt, :nEvals=>nEvals) 
  	return eqSoc, res   
  else
	@printf("Optimal policy could not solved!\n\n")
    eqSoc = EqObj()
	return eqSoc, ret
  end	
end

######################################################################################################################

function socialPlannerFunc(socPol, p, eqSocInit, cevMax, gBase, cactivtotBase)
  if (any(socPol.<=0.0))
    cev = -Inf 
    eq  = EqObj()
  else

    eq = EqObj()   # initialize equilibrium object

    # gausslegendre for integration
    eq.qMinAll = 1.0e-8
    eq.qMaxAll = 10.0
    eq.node,eq.weight = gausslegendre(1000)
    eq.node = ((eq.qMaxAll - eq.qMinAll)/2)*eq.node .+ (eq.qMaxAll + eq.qMinAll)/2    # adjust the limits


    # social planner policy variables 
    eq.qmin = socPol[1:2]
    eq.x    = socPol[3:4]
    

    res = nlsolve((eqnd,x)->socEqFunc!(eqnd,x,eq,p), eqSocInit, method = :trust_region,inplace = true,show_trace=false,iterations=50)
    
    if !res.f_converged
        #println("Couldn't solve!!!")	
        cev = -Inf
    else
        socEqFunc!(similar(res.zero),res.zero,eq,p)   # evaluate at minimizer
        eqSocInit .= res.zero
        # welfare calculation
        cev = calCev(gBase, cactivtotBase, eq.g, eq.cactivtot, p)
    end
    
    if cev>cevMax[1]
      @printf("Welfare improved to %3.4f (xˡ = %3.4f, xʰ = %3.4f, qminˡ = %3.4f, qminʰ = %3.4f) \n",
               cev,eq.x[1],eq.x[2],eq.qmin[1],eq.qmin[2])
      cevMax[1] = cev
    end
  end

  return -cev, eq


end

######################################################################################################################

function socEqFunc!(eqnd,eqv,eq,p)
  # reject out of bounds

  if (any(eqv.<=0.0))
    eqnd = -1000.0*ones(size(eqv)[1])

  else
    # load in guesses for solution
    eq.qbar   = eqv[1]
    eq.cactiv = eqv[2:3]
    eq.xout   = eqv[4]
   
    # flow of free products
    eq.cactivtot = sum(eq.cactiv)
    eq.cinac     = 1.0 - eq.cactivtot
    eq.ws        = 0.0                   # just to make things work

    # aggregate creative destruction rates
    eq.τs = eq.cactiv.*eq.x .+ p.𝛂*eq.xout
    eq.τ  = sum(eq.τs)
    # find updated eqvars
    labordem!(eq,p)
    qualityDist!(eq,p)
    qbarActive!(eq,p)


    # equilibrium differences
    eqnd[1]   = eq.qbarAct     - (p.ε/(p.ε-1.0))^(p.ε - 1)
    eqnd[2:3].= eq.cactiv     .- eqv[2:3]       # product shares
    eqnd[4]   = eq.skilled_lab - p.Lˢ           # skilled labor market clearing

  end
end
