
function policy_opt(p ;sⁱ = false, sᶠ = false, sᴱ = false)
  poltype = [sⁱ, sᶠ, sᴱ]
  npol    = sum(poltype)
  polinit = ones(npol)*.05

  
# first solve without policy  
  eqInit = [ 2.00, 0.55, 0.06,
             0.72, 0.83, 1.73]

  eq,res = solveBGP(p,eqInit) 
  eqInit.= res.zero	

  cevMax = [calCev(eq.g, eq.cactivtot,eq.g, eq.cactivtot,p)]	
  @printf("OPTIMAL POLICY\n")
  @printf("----------------------------\n")
  @printf("Baseline Welfare is %3.4f \n", cevMax[1])

  # optimization
  # ------------------------------------------------------------------------------------------------
  opt = Opt(:LN_NELDERMEAD, npol)
  ftol_rel!(opt,1e-8)
  min_objective!(opt, (x, grad) -> policy_obj(x, poltype, p, eqInit, cevMax, eq.g, eq.cactivtot))
  minf,minpol,ret = NLopt.optimize(opt, polinit)
  nEvals = opt.numevals
  #-------------------------------------------------------------------------------------------------
  
  if ret == :FTOL_REACHED 
  	@printf("Optimal policy found!\n\n")	
  	optpolval = zeros(3)
    optpolval[poltype] .= minpol
	res = Dict(:sⁱ=> optpolval[1], :sᶠ=>optpolval[2], :sᴱ=>optpolval[3], :cev=>-minf, :ret=>ret, :opt=>opt, :nEvals=>nEvals) 
  	return res   
  else
	@printf("Optimal policy could not solved!\n\n")
	return ret
  end	
		

end

function policy_obj(polguess, poltype, p, eqInit, cevMax, gBase, cactivtotBase)
 
  polval = zeros(3)
  polval[poltype] .= polguess

  pNew = Params(p, (sⁱ = polval[1], sᶠ = polval[2], sᴱ = polval[3]))
  
  eq, res = solveBGP(pNew, eqInit)
  
  if res.f_converged		
  	eqInit .= res.zero
  	# welfare calculation
  	cev = calCev(gBase, cactivtotBase, eq.g, eq.cactivtot, pNew)
  else
    cev = -Inf
  end 
  if cev>cevMax[1]
  	@printf("Welfare improved to %3.4f (", cev)
	poltype[1] ? @printf("sⁱ = %2.2f", pNew.sⁱ) : nothing
	poltype[1]&(poltype[2] | poltype[3]) ? @printf(",") : nothing	
	poltype[2] ? @printf("sᶠ = %2.2f", pNew.sᶠ) : nothing
	poltype[2]&poltype[3] ? @printf(",") : nothing	 
	poltype[3] ? @printf("sᴱ = %2.2f", pNew.sᴱ) : nothing
	@printf(")\n") 
	cevMax[1] = cev
  end
  return -cev
end

function calCev(gBase, cactivtotBase, gPol, cactivtotPol, p)
#	welf  = (1.0/(1.0-p.σ))*((eq.cactivtot^((1.0 - p.σ)/(p.ε-1.0))/(p.ρ-(1.0-p.σ)*eq.g))-(1.0/p.ρ))
	cev  = ((p.ρ-(1.0-p.σ)*gBase)/(p.ρ-(1.0-p.σ)*gPol))^(1.0/(1.0-p.σ))*((cactivtotPol/cactivtotBase)^(1.0/(p.ε-1.0)));
end
