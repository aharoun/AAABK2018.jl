module AAABK_Baseline

using NLsolve
using LinearAlgebra
using OrdinaryDiffEq
using QuadGK
using DelayDiffEq
using FastGaussQuadrature


struct Params
    λ :: Float64
    ψ :: Float64
    ν :: Float64
    α :: Float64
    ϕ :: Float64
    θL:: Float64
    θH:: Float64
    θE:: Float64
    ε :: Float64
    ρ :: Float64
    γ :: Float64
    γE:: Float64
    σ :: Float64
    Ls:: Float64
    ω :: Float64

    # transformed
    Π :: Any
    𝛎 :: Any
    𝛂 :: Any
    𝚯 :: Any
    #
    𝚯ns :: Any
    θEns:: Any
    ϕns :: Any

end




function Params(λ,ψ,ν, α, ϕ, θL, θH, θE, ε, ρ, γ, γE, σ, Ls, ω; si = 0.0,sf = 0.0, sE = 0.0)
    # with policy
    Π = (1/(ε-1))*((ε-1)/ε)^ε
    𝛎 = [0.0, ν]
    𝛂 = [1.0 - α, α]
    𝚯 = [θL, θH]

    𝚯ns  = 𝚯
    θEns = θE
    ϕns  = ϕ


    # R&D subsidy
    𝚯 = ((1.0/(1.0 - si))^((1.0 - γ)/γ)).*𝚯

    # fixed cost subsidy (exit tax)
    ϕ = (1.0 - sf)*ϕ;

    # entry subsidy
    θE= ((1.0/(1.0-sE))^((1.0 - γE)/γE))*θE

    Params(λ, ψ, ν, α, ϕ, θL, θH, θE, ε, ρ, γ, γE, σ, Ls, ω, Π, 𝛎, 𝛂, 𝚯, 𝚯ns, θEns, ϕns)
end

mutable struct EqObj
             ws:: Float64
         cactiv:: Array{Float64,1}
            eyq:: Array{Float64,1}
           qbar:: Float64
      cactivtot:: Float64
          cinac:: Float64
              x:: Array{Float64,1}
           xout:: Float64
         optval:: Array{Float64,1}
           taus:: Array{Float64,1}
            tau:: Float64
           qmin:: Array{Float64,1}
              g:: Float64
        solFAll:: ODESolution
          solFH:: ODESolution
          PhiHG:: Float64
          PhiLG:: Float64
       solFRest:: ODESolution
        PhiHExo:: Float64
        PhiLExo:: Float64
       PhiHObso:: Float64
       PhiLObso:: Float64
        qbarAct:: Float64
              r:: Float64
           cfix:: Float64
             wu:: Float64
             cx:: Array{Float64,1}
           crnd:: Float64
           cout:: Float64
    skilled_lab:: Float64
        FAllend:: Float64
        qMinAll:: Float64
        qMaxAll:: Float64
        node   :: Array{Float64,1}
        weight :: Array{Float64,1}

    EqObj() = new()
end

include("functions.jl")

export Params, EqObj, solveBGP,eqfunc





end
