
function policy_opt(p ;sⁱ = false, sᶠ = false, sᴱ = false)
  poltype = [sⁱ, sᶠ, sᴱ]
  npol    = sum(poltype)
  polinit = ones(npol)*.05

  
# first solve without policy  
  eqInit = [ 1.91, 0.55, 0.06,
             0.72, 0.83, 1.73]

  eq,res = solveBGP(p,eqInit) 
  eqInit.= res.zero	

  welfMax = [calWelfare(p, eq)]	
  
  @printf("Baseline Welfare is %3.4f \n\n", welfMax[1])
  
  opt = Opt(:LN_NELDERMEAD, npol)
  ftol_rel!(opt,1e-8)
  min_objective!(opt, (x, grad) -> policy_obj(x, poltype, p, eqInit, welfMax))
  minf,minpol,ret = NLopt.optimize(opt, polinit)
  nEvals = opt.numevals

  optpolval = zeros(3)
  optpolval[poltype] .= minpol

  @printf("Optimal policy found!\n\n")
  		
  return Dict(:sⁱ => optpolval[1] , :sᶠ => optpolval[2], :sᴱ => optpolval[3]), (polinit=polinit, minpol=minpol, minf=minf, ret=ret, opt=opt, nEvals=nEvals) 

end

function policy_obj(polguess, poltype, p, eqInit, welfMax)
  
  polval = zeros(3)
  polval[poltype] .= polguess

  pNew = Params(p.λ,  p.ψ,  p.ν,  p.α,  p.ϕ, 
                p.θˡ, p.θʰ, p.θᴱ, p.ε,  p.ρ, 
                p.γ,  p.γᴱ, p.σ,  p.Lˢ, p.ω, 
                sⁱ = polval[1], sᶠ = polval[2], sᴱ = polval[3])
  
  eq, res = solveBGP(pNew, eqInit)

  eqInit .= res.zero

  # welfare calculation
  welf  = (1.0/(1.0-p.σ))*((eq.cactivtot^((1.0 - p.σ)/(p.ε-1.0))/(p.ρ-(1.0-p.σ)*eq.g))-(1.0/p.ρ))
  if welf>welfMax[1]
  	@printf("Welfare improved to %3.4f (sⁱ = %2.2f, sᶠ = %2.2f, sᴱ = %2.2f)\n", welf, pNew.sⁱ, pNew.sᶠ, pNew.sᴱ)
	welfMax[1] = welf
  end
  return -welf

end

function calWelfare(p, eq)
	welf  = (1.0/(1.0-p.σ))*((eq.cactivtot^((1.0 - p.σ)/(p.ε-1.0))/(p.ρ-(1.0-p.σ)*eq.g))-(1.0/p.ρ))
end
