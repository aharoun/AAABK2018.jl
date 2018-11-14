module AAABK_BaselineModel

using NLsolve
using LinearAlgebra
using OrdinaryDiffEq
using DelayDiffEq
using FastGaussQuadrature


struct Params
    λ :: Float64
    ψ :: Float64
    ν :: Float64
    α :: Float64
    ϕ :: Float64
    θˡ:: Float64
    θʰ:: Float64
    θᴱ:: Float64
    ε :: Float64
    ρ :: Float64
    γ :: Float64
    γᴱ:: Float64
    σ :: Float64
    Lˢ:: Float64
    ω :: Float64

    # transformed
    Π :: Any
    𝛎 :: Any
    𝛂 :: Any
    𝚯 :: Any
    #
    𝚯ⁿˢ :: Any
    θᴱⁿˢ:: Any
    ϕⁿˢ :: Any

end



function Params(λ,ψ,ν, α, ϕ, θˡ, θʰ, θᴱ, ε, ρ, γ, γᴱ, σ, Lˢ, ω; sⁱ = 0.0, sᶠ = 0.0, sᴱ = 0.0)
    # with policy
    Π = (1/(ε-1))*((ε-1)/ε)^ε
    𝛎 = [0.0, ν]
    𝛂 = [1.0 - α, α]
    𝚯 = [θˡ, θʰ]

    𝚯ⁿˢ  = 𝚯
    θᴱⁿˢ = θᴱ
    ϕⁿˢ  = ϕ


    # R&D subsidy
    𝚯 = ((1.0/(1.0 - sⁱ))^((1.0 - γ)/γ)).*𝚯

    # fixed cost subsidy (exit tax)
    ϕ = (1.0 - sᶠ)*ϕ;

    # entry subsidy
    θᴱ= ((1.0/(1.0-sᴱ))^((1.0 - γᴱ)/γᴱ))*θᴱ

    Params(λ, ψ, ν, α, ϕ, θˡ, θʰ, θᴱ, ε, ρ, γ, γᴱ, σ, Lˢ, ω, Π, 𝛎, 𝛂, 𝚯, 𝚯ⁿˢ, θᴱⁿˢ, ϕⁿˢ)
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
             τs:: Array{Float64,1}
              τ:: Float64
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
        qAllDist::Array{Float64,1}

    EqObj() = new()
end

include("functions.jl")

export Params, EqObj, solveBGP,eqfunc!





end
