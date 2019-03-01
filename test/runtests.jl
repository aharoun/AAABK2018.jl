using AAABK2018, Test
λ =  1.131874350361868; ψ =  0.037175089899834; ν =  0.205622086535123
α =  0.925520992611854; ϕ =  0.216311889844990; θˡ=  1.390807772162332
θʰ=  1.751156379867318; θᴱ=  0.023600456924424; ε =  2.900000000000000
ρ =  0.020000000000000; γ =  0.500000000000000; γᴱ=  0.500000000000000
σ =  2.000000000000000; Lˢ=  0.165500000000000; ω =  1.0

p = Params(λ,ψ,ν, α, ϕ, θˡ, θʰ, θᴱ, ε, ρ, γ, γᴱ, σ, Lˢ, ω)

# initial guess for equilibirum objects
eqInit = [ 1.3, 0.3, 0.1,
           0.5, 0.6, 1.5]

# solve
eq,res = solveBGP(p,eqInit)

@test res.f_converged | res.x_converged

# solve with 5% incumbent R&D subsidy
pIncSub = Params(p, (sⁱ = 0.05,))
eqIncSub,resIncSub = solveBGP(pIncSub, eqInit)

@test resIncSub.f_converged | resIncSub.x_converged

# solve with 5% operation cost subsidy
pFixSub = Params(p, (sᶠ = 0.05,))
eqFixSub,resFixSub = solveBGP(pFixSub,eqInit)

@test resFixSub.f_converged | resFixSub.x_converged

# solve with 5% entrant R&D subsidy
pEntSub = Params(p, (sᴱ= 0.05,))
eqEntSub,resEntSub = solveBGP(pEntSub,eqInit)

@test resEntSub.f_converged | resEntSub.x_converged

@test all(eqIncSub.x .> eq.x)
@test all(eqFixSub.x .< eq.x)
@test all(eqEntSub.x .< eq.x)

# optimal policy
# one-tool policies
optPol = policy_opt(p, sⁱ = true)
@test optPol[:ret] == :FTOL_REACHED 
@test abs(optPol[:sⁱ] - 0.39) < .01    # compare with results from the paper
@test abs(optPol[:cev]/1.0122 - 1.0) < 1.0e-4    # compare with results from the paper

optPol = policy_opt(p, sᶠ = true)
@test optPol[:ret] == :FTOL_REACHED 
@test abs(optPol[:sᶠ] + 0.69) < .01    # compare with results from the paper
@test abs(optPol[:cev]/1.0142 - 1.0) < 1.0e-4    # compare with results from the paper

optPol = policy_opt(p, sᴱ = true)
@test optPol[:ret] == :FTOL_REACHED 
@test abs(optPol[:sᴱ] - 0.18) < .01    # compare with results from the paper
@test abs(optPol[:cev]/1.0004 - 1.0) < 1.0e-4    # compare with results from the paper

# two-tool policy
optPol = policy_opt(p, sⁱ = true, sᶠ = true)
@test optPol[:ret] == :FTOL_REACHED 
@test abs(optPol[:sⁱ] + 0.03) < .01    # compare with results from the paper
@test abs(optPol[:sᶠ] + 0.74) < .01    # compare with results from the paper
@test abs(optPol[:cev]/1.0142 - 1.0) < 1.0e-4    # compare with results from the paper
