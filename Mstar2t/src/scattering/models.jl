# *************************************************************************** #
# *                                                                         * #
# *         Mstar2t - Central Michigan University University, 2023          * #
# *                                                                         * #
# *************************************************************************** #
#  This file is part of Mstar2t.                                              #                        
#                                                                             #
#  Mstar2t is free software: you can redistribute it and/or modify it under   #
#  the terms of the GNU General Public License as published by the Free       #
#  Software Foundation, either version 3 of the License, or (at your option)  #
#  any later version.                                                         #
#                                                                             #
#  Mstar2t is distributed in the hope that it will be useful, but WITHOUT     #
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or      #
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for   #
#  more details.                                                              #
#                                                                             #
#  You should have received a copy of the GNU General Public License along    #
#  with this program. If not, see <http://www.gnu.org/licenses/>.             #
#                                                                             #
# *************************************************************************** #


@doc raw"""
    constant(A::T=1.0) where {T<:Real}

Model for the relaxation time is set to the constant relaxation time approximation.
τ model: τ = A*1e-14

Parameter:
A: value for the A constant above.

Example:
τ_model = constant()
"""
function constant(A::T=1.0) where {T<:Real}
    τ() = A*1e-14
    return Constant(τ=τ)
end


@doc raw"""
    T_fun(τ::Function)

Set the **temperature dependence** of the relaxation time from a function given by the user.  

Example:
f(T) = sqrt(T)  # τ ∝ √T
τ_model = T_fun(f)
"""
function T_fun(τ::Function)
    return TFunc(τ=τ)
end


@doc raw"""
    impurity(ϵ_im::Real, A_im::Real=1;γ::Real=1.)

Include **impurity scattering** in the simulations.

Parameters:
ϵ_im: energy of the impurity in eV
A: multiplicative constant in front of the functional expression
γ: γ-parameter three-parameter Lorentzian function (Ref: https://en.wikipedia.org/wiki/Cauchy_distribution)

Example:
ϵ_im = 0.2  # eV
τ_model = impurity(ϵ_im,A=1e-1)
"""
function impurity(ϵ_im::Real, A_im::Real=1; γ::Real=1.)

    # allow for (relatively) reasonable parameters in the interface
    A_im *= 1e13
    γ *= mu0
    ϵ_im *= mu0

    τ = function(μ::Float64)
        x = (μ-ϵ_im)
        P = A_im*γ^2/(x^2+γ^2) # scattering rate
        return 1/P
    end

    return Impurity(τ=τ)
end


@doc raw"""
    acoustic(bands_min::Union{BandStructure,Real}, A_sm::Real=1., τm_max::Real=1.; T₀::Real=50., μ_min::Real=2, μ_max::Real=2, s=2)

Include **acoustic scattering** in the simulations. Implementation of the functional expression in Wilson, Alan Herries. 1953, The theory of metals / by A. H. Wilson  Cambridge Uni. Press Cambridge, England.

Parameters:
bands_min: energy of the lowest band. If a BandStructure type is passed, the value is automatically derived. 
A_sm: multiplicative constant in front of the functional expression for the acoustic τ in semiconductors in units of 5e-20 s.
τm_max: free parameter to constraint the functional form for the acoustic τ in metals. $τ(bandmin+μ_max,T) = τm_max$ in units of 2e-12 s$ 
T₀: minimum of temperature range in which the τ is defined (note: τ ∝ 1/(T-T₀), where T is the temperature input defined by the user.
μ_min: left shift in energy from the energy of the lowest band (`ϵ_min`). It defines the first point at which τ is computed (i.e., $τ \propto 1/sqrt(μ-(ϵ_min-μ_min)$).
μ_max: right shift in energy from the energy of the band `ϵ₀`. It defines the last point at which τ is computed (i.e., $τ(ϵ₀+μ_max,T) = τm_max$).

Example:
band_1 = ParabBand([5.0, 5.0, 5.0, 0.0, 0.0, 0.0],1.0, 1,1);    # conduction band
band_2 = ParabBand([0.1, 0.5, 3.0, 0.0, 0.0, 0.0],0.5,-1,1);    # valence band
model = BandStructure(2,[band_1,band_2],0.8);   # build the two-band structure

τ_form = Scattering.acoustic(model);
"""
function acoustic(bands_min::Union{BandStructure,Real}, A_sm::Real=1., τm_max::Real=1.; T₀::Real=50., μ_min::Real=2, μ_max::Real=2, s=2)
    # allow for (relatively) reasonable parameters in the interface
    A_sm *= 5e-20
    τm_max  *= 2e-12

    local ϵ_min
    if isa(bands_min,BandStructure) 
        # set energy reference at the lowest band
        ϵ_min = get_enlowestband(bands_min)  
    elseif isa(bands_min,Real)
        ϵ_min = bands_min
    end

    # [eV] -> [J]
    ϵ_min *= mu0
    μ_min *= mu0
    μ_max *= mu0

    τ_sigmoid = function(t::Union{Float64,Int64},μ::Float64,ϵ₀::Float64,bandtype::Int64)
        
        if (t <= T₀)
            error("[ERROR] Move T₀ parameter in acoustic scattering definition.")
        end

        B_sm = 0
        # B_m is such that the two functions match at μ=ϵ₀
        B_m = A_sm/sqrt(ϵ₀-(ϵ_min-μ_min)) + B_sm
        # A_m is such that f_m(μ_max+ϵ₀) = τm_max
        A_m = (τm_max*(t-T₀)-B_m)/(μ_max)^(3/2)

        # acoustic τ for semiconductors
        f_sm(μ::Float64,t::Union{Float64,Int64}) = μ-ϵ₀ <= 0 ? (A_sm/sqrt(μ-(ϵ_min-μ_min)) + B_sm)/(t-T₀) : 0.0
        # acoustic τ for metals
        f_m(μ::Float64,t::Union{Float64,Int64}) = μ-ϵ₀ > 0 ? (A_m*(μ-ϵ₀)^(3/2)+B_m)/(t-T₀) : 0.0
        # a sigmoid is used to smooth the two functions
        sigmoid(μ::Float64) =  1/(1+exp(s*(μ-ϵ₀)))

        # conduction band
        if bandtype == 1
            return (sigmoid(μ))*f_m(μ,t) + (1-sigmoid(μ))*f_sm(μ,t)
        # in our formalism, a valence band is a conduction band with the y-axis flipped
        elseif bandtype == -1 
            return (sigmoid(ϵ₀-(μ-ϵ₀)))*f_m(ϵ₀-(μ-ϵ₀),t) + (1-sigmoid(ϵ₀-(μ-ϵ₀)))*f_sm(ϵ₀-(μ-ϵ₀),t)
        end
    end
    return Acoustic(τ=τ_sigmoid)
end


@doc raw"""
    matthiessen(τ_models::Array{ScModel};γ::Float64=-1.)
 
This function applies Matthiessen's rule to sum up different scattering mechanisms. 

Parameters:
τ_models: vector of relaxation time models
γ: exponent in the generalized Matthiessen's rule:
$τ = (\sum_i τ_i^γ)^(1/γ)$,
where each $\tau_i$ can be a function of $T$ and/or $μ$. 

Example:
ϵ_im = 0.2  # eV
τ1 = constant()
τ2 = impurity(ϵ_im,A=1e-1)
τ_model =  matthiessen([τ1,τ2])
"""
function matthiessen(τ_models::Array{ScModel};γ::Real=-1.)
    return Matthiessen(τ_models=τ_models,γ=γ)
end



# single interface to compute the relaxation time for each specific model
function compute_τ(scm::ScModel,t::Float64,μ::Float64,ϵ₀::Float64,bandtype::Int64)
    
    if scm isa Constant
        return scm.τ()
    elseif scm isa TFunc
        return scm.τ(t)
    elseif scm isa Impurity
        return scm.τ(μ)
    elseif scm isa Acoustic
        return scm.τ(t,μ,ϵ₀,bandtype)
    else
        error("[ERROR] Unknown relaxation time model.")
    end

end


# function that applies the Matthiessen's rule for a given values of relaxation time models
function matthiessen_rule(matt::Matthiessen,t::Float64,μ::Float64,ϵ₀::Float64=nothing,bandtype::Int64=nothing)
    τ_models = matt.τ_models
    γ = matt.γ
    sum_m = 0.0
    for i in eachindex(τ_models)
        sum_m += (compute_τ(τ_models[i],t,μ,ϵ₀,bandtype))^γ
    end
	return (sum_m)^(1/γ)
end
