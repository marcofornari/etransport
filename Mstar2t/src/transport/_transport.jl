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


# collection of elapsed time for each method
time = Dict{String, Float64}([("SigmaMultiBand", 0.0), 
                                ("SeebeckMultiBand", 0.0), 
                                ("CarrierConcentrationMultiBand", 0.0),
                                ("ThermalConductivityMultiBand", 0.0),
                                ("LorentzTensorMultiBand", 0.0),
                                ("HallMultiBand", 0.0),
                                ("cTensorSingle", 0.0),
                                ("EulerA", 0.0), 
                                ("EllipsC", 0.0), 
                                ("Integrate", 0.0)])


# this function performs all the steps necessary before running the calculation.
# 1. Rescale the unit of measure
# 2. Construct an iterator over all the combinations of band structures
# 3. Check temperature range
function precalculation!(num_bands::Int64,ebandmins::Union{Vector{Vector{Float64}},Vector{Any},Vector{Float64}},μs::Union{Vector{Float64},Float64},Ts::Union{Float64,Array{Float64},Vector{Int64},Int64})
    μs *= mu0           # rescale the unit [eV] -> [J]
    ebandmins *= mu0    # rescale the unit [eV] -> [J]
    
    # get iterator over the combinations for all the bands minima
    if length(ebandmins) > 1 
        it_e = combination_energy(num_bands, ebandmins)
    else
        it_e = ebandmins[1]
    end

    if length(Ts) > 1
        filter!(x -> x>9, Ts)    # filter temperature to avoid numerical instability
        Ts = convert(Array{Float64},Ts)
    end
    βs = ./(1,(kB*Ts))

    return it_e,μs,Ts,βs
end


@doc raw"""
    electrical_conductivity(model::BandStructure,Ts::Union{Vector{Float64},Float64,Vector{Int64},Int64},τ::Union{ScModel,Matthiessen}; exportasdf::Bool=false,fulltensor::Bool=false)
 
Function that by default computes the trace of the tensorial version of the electical conductivity in units of $(\Omega m)^{-1}$ for a given `bandstructure` and a choice of the relaxation time, for an array of `fermi_level` and `temperature`. The transport coefficient is returned as a `matrix` of dimensions length(T)*length(μ). The boolean variable `exportasdf` allows to return the calculations as a `DataFrame` with all the parameters included. If `fulltensor` is set to `true` the full tensor is returned in place of the trace.
Reference: "Theory of band warping and its effects on thermoelectronic transport properties", PHYSICAL REVIEW B89.
"""
function electrical_conductivity(model::BandStructure,Ts::Union{Vector{Float64},Float64,Vector{Int64},Int64},τ::Union{ScModel,Matthiessen}; exportasdf::Bool=false,fulltensor::Bool=false)

    # extract params of the calculation from the model
    cTensorParameters,ebandmins,bandtypes,degeneracies,μs = extractparams(model)    

    num_bands = length(bandtypes)   # number of bands
    it_e,μs_J,Ts,βs = precalculation!(num_bands,ebandmins,μs,Ts)   # convert units, check temp range

    # compute the electical conductivity
    σ = compute_elcond(num_bands,cTensorParameters,it_e,bandtypes,degeneracies,μs_J,Ts,βs,τ)

    # return the results
    if !fulltensor
        if !exportasdf
            return to_matrix(σ,length(it_e),length(μs),length(Ts))
        else
            return to_df(σ,"sigma",cTensorParameters,it_e,bandtypes,degeneracies,μs,Ts)
        end
    else
        return fulltensor_df(σ,"sigma",cTensorParameters,it_e,bandtypes,degeneracies,μs,Ts)
    end
end

# calculation of the electical conductivity
function compute_elcond(num_bands::Int64,cTensorParameters::Vector{Vector{Float64}},it_e::Union{Float64,Vector{Array{Float64}}},bandtypes::Vector{Int64},degeneracies::Vector{Int64},μs_J::Union{Vector{Float64},Float64},Ts::Union{Vector{Float64},Float64},βs::Union{Vector{Float64},Float64},τ::Union{ScModel,Matthiessen})
    sSOTMB = time_ns()
    σ = SharedArray{Float64}(3,3,length(βs),length(it_e),length(μs_J))
    @sync Threads.@threads for m in eachindex(μs_J)     # for each chemical potential
        μ = μs_J[m]
        for (e,en) in enumerate(it_e)   # for each combination of band minima
            for t in eachindex(βs)   # for each temperature
                β = βs[t]
                T = Ts[t]
                # compute the KinCoeff and then the electical conductivity
                L₀ = zeros(3,3)
                # sum over all bands
                for b in 1:num_bands
                    bandtype,deg,ϵ₀ = bandtypes[b],degeneracies[b],en[b]
                    pf = .*(cTensorSingle(cTensorParameters[b]),bandtype,deg)
                    L₀ = .+(L₀, pf*Ln(ϵ₀, μ, bandtype, T, β, 0, 1, τ))
                end
                σ[:,:,t,e,m] = *(qe^2, L₀)  # store the result
            end
        end
    end
    time["SigmaMultiBand"] += time_ns() - sSOTMB
    return σ
end


@doc raw"""
    seebeck_coefficient(model::BandStructure, Ts::Union{Vector{Float64},Float64,Vector{Int64},Int64},τ::Union{ScModel,Matthiessen};exportasdf::Bool=false,fulltensor::Bool=false)
 
Function that by default computes the Seebeck coefficient in units of $V/K$ for a given `bandstructure` and a choice of the relaxation time, for an array of `fermi_level` and `temperature`. The transport coefficient is returned as a `matrix` of dimensions length(T)*length(μ). The boolean variable `exportasdf` allows to return the calculations as a `DataFrame` with all the parameters included. If `fulltensor` is set to `true`, the full tensor is returned in place of the trace.
Reference: "Theory of band warping and its effects on thermoelectronic transport properties", PHYSICAL REVIEW B89.
"""
function seebeck_coefficient(model::BandStructure, Ts::Union{Vector{Float64},Float64,Vector{Int64},Int64},τ::Union{ScModel,Matthiessen};exportasdf::Bool=false,fulltensor::Bool=false)

    # extract params of the calculation from the model
    cTensorParameters,ebandmins,bandtypes,degeneracies,μs = extractparams(model)
    
    num_bands = length(bandtypes)   # number of bands
    it_e,μs_J,Ts,βs = precalculation!(num_bands,ebandmins,μs,Ts)   # convert units, check temp range

    # compute the Seebeck coefficient
    S = compute_seebeck(num_bands,cTensorParameters,it_e,bandtypes,degeneracies,μs_J,Ts,βs,τ)
    
    # return the results
    if !fulltensor
        if !exportasdf
            return to_matrix(S,length(it_e),length(μs),length(Ts))
        else
            return to_df(S,"S",cTensorParameters,it_e,bandtypes,degeneracies,μs,Ts)
        end
    else
        return fulltensor_df(S,"S",cTensorParameters,it_e,bandtypes,degeneracies,μs,Ts)
    end
end


# calculation of the Seebeck coefficient
function compute_seebeck(num_bands::Int64,cTensorParameters::Vector{Vector{Float64}},it_e::Union{Float64,Vector{Array{Float64}}},bandtypes::Vector{Int64},degeneracies::Vector{Int64},μs_J::Union{Vector{Float64},Float64},Ts::Union{Vector{Float64},Float64},βs::Union{Vector{Float64},Float64},τ::Union{ScModel,Matthiessen})
    sSMB = time_ns()
    S = SharedArray{Float64}(3,3,length(βs),length(it_e),length(μs_J))
    @sync Threads.@threads for m in eachindex(μs_J)     # for each chemical potential
        μ = μs_J[m]
        for (e,en) in enumerate(it_e)   # for each combination of band minima
            for t in eachindex(βs)   # for each temperature
                β = βs[t]
                T = Ts[t]
                # compute the KinCoeff and then the Seebeck coefficient
                L₀,L₁ = zeros(3,3),zeros(3,3)
                # sum over all bands
                for b in 1:num_bands
                    bandtype,deg,ϵ₀ = bandtypes[b],degeneracies[b],en[b]
                    pf = .*(cTensorSingle(cTensorParameters[b]),bandtype,deg)
                    L₀ = .+(L₀, pf*Ln(ϵ₀, μ, bandtype, T, β, 0, 1, τ))
                    L₁ = .+(L₁, pf*Ln(ϵ₀, μ, bandtype, T, β, 1, 1, τ))
                end
                S[:,:,t,e,m] = *((-kB*β)/qe,*(inv(L₀),L₁))  # store the result
            end
        end
    end
    time["SeebeckMultiBand"] += time_ns() - sSMB
    return S
end

@doc raw"""
    carrier_concentration(model::BandStructure, Ts::Union{Vector{Float64},Float64,Vector{Int64},Int64},τ::Union{ScModel,Matthiessen}; exportasdf::Bool=false)
 
Function that computes the carrier concentration for a given `bandstructure` and a choice of the relaxation time, and for an array of `fermi_level` and `temperature`. The transport coefficient is returned as a `matrix` of dimensions length(T)*length(μ). The boolean variable `exportasdf` allows to return the calculations as a `DataFrame` with all the parameters included.
"""
function carrier_concentration(model::BandStructure, Ts::Union{Vector{Float64},Float64,Vector{Int64},Int64},τ::Union{ScModel,Matthiessen}; exportasdf::Bool=false)

    # extract params of the calculation from the model
    cTensorParameters,ebandmins,bandtypes,degeneracies,μs = extractparams(model)
    
    num_bands = length(bandtypes)   # number of bands
    it_e,μs_J,Ts,βs = precalculation!(num_bands,ebandmins,μs,Ts)   # convert units, check temp range

    # compute the carrier concentration
    n = compute_carrierconc(num_bands,cTensorParameters,it_e,bandtypes,degeneracies,μs_J,Ts,βs,τ)

    # return the results
    if !exportasdf
        return to_matrix(n,length(it_e),length(μs),length(Ts))
    else
        return to_df(n,"n",cTensorParameters,it_e,bandtypes,degeneracies,μs,Ts)
    end
end


# calculation of the carrier concentration
function compute_carrierconc(num_bands::Int64,cTensorParameters::Vector{Vector{Float64}},it_e::Union{Float64,Vector{Array{Float64}}},bandtypes::Vector{Int64},degeneracies::Vector{Int64},μs_J::Union{Vector{Float64},Float64},Ts::Union{Vector{Float64},Float64},βs::Union{Vector{Float64},Float64},τ::Union{ScModel,Matthiessen})
    sCC = time_ns()
    # compute the Hall tensor and then cc
    n = SharedArray{Float64}(length(βs),length(it_e),length(μs_J)) 
    @sync Threads.@threads for m in eachindex(μs_J)     # for each chemical potential
        μ = μs_J[m]
        for (e,en) in enumerate(it_e)   # for each combination of band minima
            for t in eachindex(βs)   # for each temperature
                β = βs[t]
                T = Ts[t]
                # sum over all bands
                nᵢ = Vector{Float64}(undef,num_bands)
                for b in 1:num_bands
                    bandtype,deg,ϵ₀,ctparams = bandtypes[b],degeneracies[b],en[b],cTensorParameters[b]
                    cT = cTensorSingle(ctparams)
                    pf = bandtype*deg
                    L₀ = cT*pf*Ln(ϵ₀, μ, bandtype, T, β, 0, 1, τ)
                    invsigma = inv(*(qe^2, L₀))
                    σ² = .*(sspref*(-bandtype)*pf, Ln(ϵ₀, μ, bandtype, T, β, 0, 2, τ), cTensorkron(ctparams,cT)) # sspref: qe^3/(me*c)
                    Rsingle = zeros(3,3,3)  # R_{i,j,k} = σ^{-1}_{mj}⋅σ_{mnk}⋅σ_{in}
                    Einsum.@einsum Rsingle[1,2,3] = invsigma[2,n] * σ²[n,m,3] * invsigma[m,1]
                    nᵢ[b] = 1/Rsingle[1,2,3]
                end
                n[t,e,m] = sum(qe*sl*nᵢ)  # store the result
            end
        end
    end 
    time["CarrierConcentrationMultiBand"] += time_ns() - sCC
    return n
end



@doc raw"""
    thermal_conductivity(model::BandStructure, Ts::Union{Vector{Float64},Float64,Vector{Int64},Int64},τ::Union{ScModel,Matthiessen}; exportasdf::Bool=false, fulltensor::Bool=false)

Function that computes the thermal conductivity in units of $W/K$ for a given `bandstructure` and a choice of the relaxation time, for an array of `fermi_level` and `temperature`. The transport coefficient is returned as a `matrix` of dimensions length(T)*length(μ). The boolean variable `exportasdf` allows to return the calculations as a `DataFrame` with all the parameters included. If `fulltensor` is set to `true`, the full tensor is returned in place of the trace.. 
Reference: "Theory of band warping and its effects on thermoelectronic transport properties", PHYSICAL REVIEW B89.
"""
function thermal_conductivity(model::BandStructure, Ts::Union{Vector{Float64},Float64,Vector{Int64},Int64},τ::Union{ScModel,Matthiessen}; exportasdf::Bool=false, fulltensor::Bool=false)
    
    # extract params of the calculation from the model
    cTensorParameters,ebandmins,bandtypes,degeneracies,μs = extractparams(model)
    
    num_bands = length(bandtypes)   # number of bands
    it_e,μs_J,Ts,βs = precalculation!(num_bands,ebandmins,μs,Ts)   # convert units, check temp range
    
    # compute the thermal conductivity
    κ = compute_thermalcond(num_bands,cTensorParameters,it_e,bandtypes,degeneracies,μs_J,Ts,βs,τ)

    # return the results
    if !fulltensor
        if !exportasdf
            return to_matrix(κ,length(it_e),length(μs),length(Ts))
        else
            return to_df(κ,"K",cTensorParameters,it_e,bandtypes,degeneracies,μs,Ts)
        end
    else
        return fulltensor_df(κ,"K",cTensorParameters,it_e,bandtypes,degeneracies,μs,Ts)
    end
end

# calculation of the thermal conductivity
function compute_thermalcond(num_bands::Int64,cTensorParameters::Vector{Vector{Float64}},it_e::Union{Float64,Vector{Array{Float64}}},bandtypes::Vector{Int64},degeneracies::Vector{Int64},μs_J::Union{Vector{Float64},Float64},Ts::Union{Vector{Float64},Float64},βs::Union{Vector{Float64},Float64},τ::Union{ScModel,Matthiessen})
    sTCMB = time_ns()
    κ = SharedArray{Float64}(3,3,length(βs),length(it_e),length(μs_J))
    @sync Threads.@threads for m in eachindex(μs_J)     # for each chemical potential
        μ = μs_J[m]
        for (e,en) in enumerate(it_e)   # for each combination of band minima
            for t in eachindex(βs)   # for each temperature
                β = βs[t]
                T = Ts[t]
                # compute the KinCoeff and then the thermal conductivity
                L₀,L₁,L₂ = zeros(3,3),zeros(3,3),zeros(3,3)
                # sum over all bands
                for b in 1:num_bands
                    bandtype,deg,ϵ₀ = bandtypes[b],degeneracies[b],en[b]
                    pf = .*(cTensorSingle(cTensorParameters[b]),bandtype,deg)
                    L₀ = .+(L₀, pf*Ln(ϵ₀, μ, bandtype, T, β, 0, 1, τ))
                    L₁ = .+(L₁, pf*Ln(ϵ₀, μ, bandtype, T, β, 1, 1, τ))
                    L₂ = .+(L₂, pf*Ln(ϵ₀, μ, bandtype, T, β, 2, 1, τ))
                end
                var = *(inv(L₀),L₁)
                var = *(L₁, var)
                var = .*(kB*β, L₂ - var)
                κ[:,:,t,e,m] = var  # store the result
            end
        end
    end

    time["ThermalConductivityMultiBand"] += time_ns() - sTCMB
    return κ
end


@doc raw"""
    lorenz_tensor(model::BandStructure, Ts::Union{Vector{Float64},Float64,Vector{Int64},Int64},τ::Union{ScModel,Matthiessen};exportasdf::Bool=false,fulltensor::Bool=false)
 
Function that computes the Lorentz tensor for a given `bandstructure` and a choice of the relaxation time, for an array of `fermi_level` and `temperature`. The Lorentz tensor is returned as a `matrix` of dimensions length(T)*length(μ). The boolean variable `exportasdf` allows to return the calculations as a `DataFrame` with all the parameters included. If `fulltensor` is set to `true`, the full tensor is returned in place of the trace. 
The Lorentz tensor is defined as the ratio between thermal and electrical conductivity multiplied by temperature. Ref.: https://en.wikipedia.org/wiki/Wiedemann-Franz_law
"""
function lorenz_tensor(model::BandStructure, Ts::Union{Vector{Float64},Float64,Vector{Int64},Int64},τ::Union{ScModel,Matthiessen};exportasdf::Bool=false,fulltensor::Bool=false)
    sLTMB = time_ns()

    # extract params of the calculation from the model
    cTensorParameters,ebandmins,bandtypes,degeneracies,μs = extractparams(model)

    num_bands = length(bandtypes)   # number of bands
    it_e,μs_J,Ts,βs = precalculation!(num_bands,ebandmins,μs,Ts)   # convert units, check temp range

    L = SharedArray{Float64}(3,3,length(βs),length(it_e),length(μs_J))
    @sync Threads.@threads for m in eachindex(μs_J)     # for each chemical potential
        μ = μs_J[m]
        for (e,en) in enumerate(it_e)   # for each combination of band minima
            for t in eachindex(βs)   # for each temperature
                β = βs[t]
                T = Ts[t]
                # compute the KinCoeff and then the Lorent number
                L₀,L₁,L₂ = zeros(3,3),zeros(3,3),zeros(3,3)
                # sum over all bands
                for b in 1:num_bands
                    bandtype,deg,ϵ₀ = bandtypes[b],degeneracies[b],en[b]
                    pf = .*(cTensorSingle(cTensorParameters[b]),bandtype,deg)
                    L₀ = .+(L₀, pf*Ln(ϵ₀, μ, bandtype, T, β, 0, 1, τ))
                    L₁ = .+(L₁, pf*Ln(ϵ₀, μ, bandtype, T, β, 1, 1, τ))
                    L₂ = .+(L₂, pf*Ln(ϵ₀, μ, bandtype, T, β, 2, 1, τ))
                end
                inv_k0 = inv(L₀)
                var = *(inv_k0, L₁)
                var = *(L₁, var)
                var = *(inv_k0, L₂ - var)
                var = .*((kB*β/qe)^2, var)
                L[:,:,t,e,m] = var  # store the result
            end
        end
    end    
    time["LorentzTensorMultiBand"] += time_ns() - sLTMB

    # return the results
    if !fulltensor
        if !exportasdf
            return to_matrix(L,length(it_e),length(μs),length(Ts))
        else
            return to_df(L,"L",cTensorParameters,it_e,bandtypes,degeneracies,μs,Ts)
        end
    else
        return fulltensor_df(L,"L",cTensorParameters,it_e,bandtypes,degeneracies,μs,Ts)
    end
end


"""
    hall(bandstructure::BandStructure, temperature, scattering_time)
 
Function that computes the Hall tensor for a given `bandstructure` and a choice of the relaxation time, for an array of `fermi_level` and `temperature`. 
"""
function hall(model::BandStructure, Ts::Union{Vector{Float64},Float64,Vector{Int64},Int64}, τ::Union{ScModel,Matthiessen})
    sHMB = time_ns()

    # extract params of the calculation from the model
    cTensorParameters,ebandmins,bandtypes,degeneracies,μs = extractparams(model)

    num_bands = length(bandtypes)   # number of bands
    it_e,μs_J,Ts,βs = precalculation!(num_bands,ebandmins,μs,Ts)   # convert units, check temp range

    # R_{i,j,k}
    R = SharedArray{Float64}(3,3,3,length(βs),length(it_e),length(μs_J))
    @sync Threads.@threads for m in eachindex(μs_J)     # for each chemical potential
        μ = μs_J[m]
        for (e,en) in enumerate(it_e)   # for each combination of band minima
            for t in eachindex(βs)   # for each temperature
                β = βs[t]
                T = Ts[t]
                L₀,Rsingle,σ² = zeros(3,3),zeros(3,3,3),zeros(3,3,3) # R_{i,j,k}, σ_{i,j,k}
                # sum over all bands
                for b in 1:num_bands
                    bandtype,deg,ϵ₀,ctparams = bandtypes[b],degeneracies[b],en[b],cTensorParameters[b]
                    cT = cTensorSingle(ctparams)
                    pf = bandtype*deg
                    L₀ = .+(L₀, pf*Ln(ϵ₀, μ, bandtype, T, β, 0, 1, τ)*cT)
                    σ² = .+(σ², .*(pf*Ln(ϵ₀, μ, bandtype, T, β, 0, 2, τ),cTensorkron(ctparams,cT),-bandtype))
                end
                invsigma = inv(*(qe^2, L₀))
                σ² = .*(sspref, σ²) # sspref: qe^3/(me*c)
                # R_{i,j,k} = σ^{-1}_{mj}⋅σ_{mnk}⋅σ_{in}
                Einsum.@einsum Rsingle[i,j,k] = invsigma[j,n] * σ²[n,m,k] * invsigma[m,i]
                R[:,:,:,t,e,m] = Rsingle  # store the result
            end
        end
    end
    time["HallMultiBand"] += time_ns() - sHMB
    return R
end


# function that performs the product and index contraction between the Levi-Civita tensor and the mass tensor $C_{αβ}$.
# Reference: "Theory of band warping and its effects on thermoelectronic transport properties", PHYSICAL REVIEW B89.
function cTensorkron(cTensorParameters::Array{Float64,1},cTensor::Array{Float64,2})
    inv_masses = ./(1,[cTensorParameters[1],cTensorParameters[2],cTensorParameters[3]])

    # each addi is a tensor with indices α,β,γ created with tensor product between vec A_α and matrix B_βγ
    add1 = reshape(kron(.*(-ϵijk[:,:,1], inv_masses), cTensor[:,1]), (3,3,3))
    add2 = reshape(kron(.*(-ϵijk[:,:,2], inv_masses), cTensor[:,2]), (3,3,3))
    add3 = reshape(kron(.*(-ϵijk[:,:,3], inv_masses), cTensor[:,3]), (3,3,3))

    return add1 + add2 + add3
end


# compute $C_{αβ}$
# Reference: "Theory of band warping and its effects on thermoelectronic transport properties", PHYSICAL REVIEW B89.
function cTensorSingle(param::Array{Float64,1})
    scTensor = time_ns()
    out = Array{Float64,2}(undef,3,3) 
    mx,my,mz,psi,theta,phi = param
    out = eulerA(psi, theta, phi) * ellipsC(mx, my, mz) * inv(eulerA(psi, theta, phi))
    time["cTensorSingle"] += time_ns() - scTensor
    return out
end


# compute rotation matrix from euler angles
function eulerA(psi::Float64, theta::Float64, phi::Float64)
    seuler = time_ns()
    ea = [cos(phi)*cos(psi) - cos(theta)*sin(phi)*sin(psi) cos(theta)*cos(phi)*sin(psi) + cos(psi)*sin(phi) sin(theta)*sin(psi); -cos(phi)*sin(psi) - cos(theta)*cos(psi)*sin(phi) cos(theta)*cos(phi)*cos(psi) - sin(phi)sin(psi) cos(psi)*sin(theta); sin(theta)*sin(phi) -cos(phi)*sin(theta) cos(theta)]
    time["EulerA"] += time_ns() - seuler
    return ea
end


function ellipsC(mx::Float64, my::Float64, mz::Float64)
    sEllipsC = time_ns()
    B = Diagonal([8.377580409572781*sqrt((my*mz)/mx),8.377580409572781*sqrt((mx*mz)/my),8.377580409572781*sqrt((mx*my)/mz)])
    time["EllipsC"] += time_ns() - sEllipsC
    return B
end


@doc raw"""
Ln(ϵ₀::Float64, μ::Float64, bandtype::Int64, T::Float64, β::Float64, n::Int64, idx::Int64, scm::ScModel)Ln(ϵ₀, μ, bandtype, T, β, n, idx, τ)
 
This function is called when an integration is required. $n=0,1,2$ is the index of the kinetic coefficient.
Reference: "Theory of band warping and its effects on thermoelectronic transport properties", PHYSICAL REVIEW B89. 
"""
function Ln(ϵ₀::Float64, μ::Float64, bandtype::Int64, T::Float64, β::Float64, n::Int64, idx::Int64, scm::ScModel)
    sF = time_ns()
    s = bandtype*β*(ϵ₀-μ)
    Lₙ = pref*bandtype^(n+1)/(β^(n + 3/2))  # constant prefactor
    τ = compute_τ(scm,T,μ,ϵ₀,bandtype)
    Lₙ *= τ^idx * k_series_int(n, 0, s) # integral calculation
    time["Integrate"] += time_ns() - sF
    return Lₙ
end


# If τ is an array (>1 scattering mechanisms) -> combine with Matthiessen rule
function Ln(ϵ₀::Float64, μ::Float64, bandtype::Int64, T::Float64, β::Float64, n::Int64, idx::Int64, scm::Matthiessen)
    sF = time_ns()
    s = bandtype*β*(ϵ₀-μ)
    Lₙ = pref*bandtype^(n+1)/(β^(n + 3/2))  # constant prefactor
    τ = matthiessen_rule(scm,T,μ,ϵ₀,bandtype)
    Lₙ *= τ^idx * k_series_int(n, 0, s) # integral calculation
    time["Integrate"] += time_ns() - sF
    return Lₙ
end