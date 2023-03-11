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


"""
    runcalculation(input_json)
    
This function is called by the server whenever a request from the client is received. It parses the JSON file containing the parameters for the computations, check if the export path chosen by the user exists, save the input parameters for future reference, construct the relaxation time functional form, and start the calculation.  
"""
function runcalculation(input_json)
    global export_all_data
    parsed_args, local_path, export_all_data, num_bands, cTensorParameters, ebandmins, bandtypes, degeneracies, μs, Ts, τ_form, τ_ac_pars, τ_im_pars = parseargs(input_json) # parse the input json
    if !isdir(local_path)   # check path exists
        println("Export path not found.")
        println("Check input file [\"results fullpath\" value].")
        return -10
    end
    # check for empty tau
    if τ_form ∉ ("constant","acoustic","impurity")
        printstyled("ERROR: "; color = :red, bold=true)
        println("Relaxation time unknown.")
        println("Check input file.")
        return -20
    end    
    τ = set_τ(τ_form, τ_ac_pars, τ_im_pars)  # set relaxation time func form
    full_test_path = create_testfolder(local_path)    # set the correct path to export data
    save_inputfile(input_json, full_test_path)    # save input file for reference
    println("--> 2.0 Starting Calculation.")
    try
        tensor_calculation(full_test_path, parsed_args, num_bands, cTensorParameters, ebandmins, bandtypes, degeneracies, μs, Ts, τ)    # start the calculation
        return full_test_path
    catch e
        if isa(e, DomainError)
            printstyled("ERROR: "; color = :red, bold=true)
            println("Domain error in the τ function calculation.")
            println("Possible reasons: negative arguments in logs, negative radicand in radical expressions.")
            printstyled("Hint: "; color = :cyan, bold=true)
            println("shift the zero value of the chosen energy scale in the construction of the band structure.")
            return -30
        else
            throw(e)
        end
    end
end


"""
    runGUIcalculation(input_json)
 
Function similar to `runcalculation` in its purposes, it handles requests from the GUI instead of from the CLI. The results are sent back to the GUI for visualization. 
"""
function runGUIcalculation(input_json)
    global export_all_data
    tensor_name, export_all_data, num_bands, cTensorParameters, bandtypes, ebandmins, degeneracies, μs, Ts, τ_form, τ_ac_pars, τ_im_pars, τ_matth_models, τ_matth_γ = GUIparseargs(input_json) # parse the input json
    # check for empty tau
    if τ_form ∉ ("constant","acoustic","impurity","Matthiessen's rule")
        printstyled("ERROR: "; color = :red, bold=true)
        println("Relaxation time unknown.")
        println("Check input file.")
        return -20
    end  
    try
        τ = set_τ(τ_form, τ_ac_pars, τ_im_pars, τ_matth_models, τ_matth_γ) # set relaxation time func form
        data = tensor_GUIcalculation(tensor_name, num_bands, cTensorParameters, ebandmins, bandtypes, degeneracies, μs, Ts, τ)    # start the calculation
        if tensor_name == "concentration" # carrier concentration
            tensor = Array{Float64,3}(undef, length(Ts), length(μs), 1)  # array to collect results
            for t in eachindex(Ts)
                for m in eachindex(μs) 
                    tensor[t,m,1] = data[t,1,m]  # populate the array with results
                end
            end
            return Dict("T" => Ts, "mu" => μs, "data" => tensor)  # send results to GUI
        else    # other tensors, i.e. electical conductivity, seebeck, thermal conductivity
            tensor = Array{Float64,3}(undef, length(Ts), length(μs), 6)  # array to collect results
            var = Array{Float64,1}(undef, 6)
            for t in eachindex(Ts)
                for m in eachindex(μs)
                    tmp_matrix = data[:,:,t,1,m]
                    var[1] = tmp_matrix[1,1] 
                    var[2] = tmp_matrix[2,2] 
                    var[3] = tmp_matrix[3,3] 
                    var[4] = tmp_matrix[1,2] 
                    var[5] = tmp_matrix[1,3] 
                    var[6] = tmp_matrix[2,3]
                    tensor[t,m,:] = var  # populate the array with results
                end
            end
            return Dict("T" => Ts, "mu" => μs, "data" => tensor)  # send results to GUI
        end
    catch e
        if isa(e, DomainError)
            printstyled("ERROR: "; color = :red, bold=true)
            println("Domain error in the τ function calculation.")
            println("Possible reasons: negative arguments in logs, negative radicand in radical expressions.")
            printstyled("Hint: "; color = :cyan, bold=true)
            println("shift the zero value of the chosen energy scale in the construction of the band structure.")
            return -30
        elseif isa(e, MatthiessenError)
            printstyled("ERROR: "; color = :red, bold=true)
            println("Please use the checkbox to select scatterings for Matthiessen's rule.")
            return -40
        else
            throw(e)
        end
    end
end



"""
set_τ(τ_form_str, τ_ac_pars, τ_im_pars, τ_matth_models, τ_matth_γ)
    
Construct the correct functional form of the relaxation time from the input file. `τ_model` is "constant", "acoustic" or "impurity". `τ_ac_pars` are the coefficients of the acoustic relaxation time. `τ_im_pars` are the coefficients of the impurity relaxation time.
"""
function set_τ(τ_form_str::String, τ_ac_pars::Array{Any}=[], τ_im_pars::Array{Any}=[], τ_matth_models::String="000", τ_matth_γ::Float64=-1.)
    if τ_form_str == "constant"
        return constant()
    elseif τ_form_str == "acoustic"
        if length(τ_ac_pars) == 0
            error("[ERROR] Please set acoustic scattering parameters.")
        end
        bands_min, A_sm, τm_max, T₀, μ_min, μ_max = τ_ac_pars
        return acoustic(bands_min, A_sm, τm_max, T₀=T₀, μ_min=μ_min, μ_max=μ_max)
    elseif τ_form_str == "manual"
        # TODO
    elseif τ_form_str == "impurity"
        if length(τ_im_pars) == 0
            error("[ERROR] Please set impurity scattering parameters.")
        end
        ϵ_im,A_im,γ_im = τ_im_pars
        return impurity(ϵ_im,A_im,γ=γ_im)
    elseif τ_form_str == "Matthiessen's rule"
        if τ_matth_models == "000"
            throw(MatthiessenError())
            # error("[ERROR] Please use the checkbox to select scatterings for Matthiessen's rule.")
        end
        bands_min, A_sm, τm_max, T₀, μ_min, μ_max = τ_ac_pars
        ϵ_im,A_im,γ_im = τ_im_pars
        γ_mr = τ_matth_γ
        τ_models = Vector{ScModel}()
        if τ_matth_models[1] == '1'
            push!(τ_models, constant())
        end
        if τ_matth_models[2] == '1'
            push!(τ_models, acoustic(bands_min, A_sm, τm_max, T₀=T₀, μ_min=μ_min, μ_max=μ_max))
        end
        if τ_matth_models[3] == '1'
            push!(τ_models, impurity(ϵ_im,A_im,γ=γ_im))
        end
        return matthiessen(τ_models,γ=γ_mr)
    end
end


"""
    tensor_calculation(path, args, num_bands, band_masses, band_min, band_type, deg, fermi_level, temp, τ)
 
This method calls the routines that compute the transport tensors according to the requests of the user.
"""
function tensor_calculation(full_test_path::String, parsed_args, num_bands::Int64, cTensorParameters::Vector{Vector{Float64}}, ebandmins::Union{Vector{Vector{Float64}},Vector{Any},Vector{Float64}}, bandtypes::Vector{Int64}, degeneracies::Vector{Int64}, μs::Union{Vector{Float64},Float64}, Ts::Union{Float64,Array{Float64},Vector{Int64},Int64}, τ::Union{ScModel,Matthiessen})
    
    global export_all_data
    local counter = 0
    local σ, S, n, κ = repeat([[]],4)

    it_e,μs_J,Ts,βs = precalculation!(num_bands,ebandmins,μs,Ts)   # convert units, check temp range

    if parsed_args["conductivity"]
        counter += 1
        print("---> 2.", counter, " Conductivity ...\r")
        σ = compute_elcond(num_bands,cTensorParameters,it_e,bandtypes,degeneracies,μs_J,Ts,βs,τ) 
        println("---> 2.", counter, " Conductivity done.")
    end
    if parsed_args["seebeck"]
        counter += 1
        print("---> 2.", counter, " Seebeck coefficient ...\r")
        if isempty(S)
            S = compute_seebeck(num_bands,cTensorParameters,it_e,bandtypes,degeneracies,μs_J,Ts,βs,τ)
        end
        println("---> 2.", counter, " Seebeck coefficient done.")
    end
    if parsed_args["concentration"]
        counter += 1
        print("---> 2.", counter, " Carrier concentration ...\r")
        n = compute_carrierconc(num_bands,cTensorParameters,it_e,bandtypes,degeneracies,μs_J,Ts,βs,τ)
        println("---> 2.", counter, " Carrier concentration done.")
    end
    if parsed_args["thermal"]
        counter += 1
        print("---> 2.", counter, " Thermal conductivity ...\r")
        if isempty(κ)
            κ = compute_thermalcond(num_bands,cTensorParameters,it_e,bandtypes,degeneracies,μs_J,Ts,βs,τ)
        end
        println("---> 2.", counter, " Thermal conductivity done.")
    end
    println("--> 2.", counter+1, " Calculation completed.")
    save_to_disk(export_all_data, full_test_path, [σ, S, n, κ], cTensorParameters, it_e, bandtypes, degeneracies, μs, Ts)
    if parsed_args["tplot"] == true
        CLI_plotresults(full_test_path, [σ, S, n, κ], it_e, μs, Ts, type="T")
    end
    if parsed_args["muplot"] == true
        CLI_plotresults(full_test_path, [σ, S, n, κ], it_e, μs, Ts, type="μ")
    end
    total_computation = *(length(it_e),length(μs_J),length(βs))
    println("[INFO] Number of runs for each tensor: ", total_computation)
end


"""
    tensor_GUIcalculation(tensor_name, num_bands, band_masses, band_min, band_type, deg, fermi_level, temp, τ)
 
This method calls the routines that compute the transport tensors according to the requests of the GUI.
"""
function tensor_GUIcalculation(tensor_name::String, num_bands::Int64, cTensorParameters::Vector{Vector{Float64}}, ebandmins::Union{Vector{Vector{Float64}},Vector{Any},Vector{Float64}}, bandtypes::Vector{Int64}, degeneracies::Vector{Int64}, μs::Union{Vector{Float64},Float64}, Ts::Union{Float64,Array{Float64},Vector{Int64},Int64},τ::Union{ScModel,Matthiessen})

    local tensor
    it_e,μs_J,Ts,βs = precalculation!(num_bands,ebandmins,μs,Ts)   # convert units, check temp range

    if tensor_name == "conductivity"
        tensor = compute_elcond(num_bands,cTensorParameters,it_e,bandtypes,degeneracies,μs_J,Ts,βs,τ)
    elseif tensor_name == "seebeck"
        tensor = compute_seebeck(num_bands,cTensorParameters,it_e,bandtypes,degeneracies,μs_J,Ts,βs,τ)
    elseif tensor_name == "concentration"
        tensor = compute_carrierconc(num_bands,cTensorParameters,it_e,bandtypes,degeneracies,μs_J,Ts,βs,τ)
    elseif tensor_name == "thermal"
        tensor = compute_thermalcond(num_bands,cTensorParameters,it_e,bandtypes,degeneracies,μs_J,Ts,βs,τ)
    end
    return tensor
end


"""
GUItaucalculation(input_json)
 
Function to compute the relaxation time for a specific range of temperature and Fermi level position. The results are sent back to the GUI for visualization. 
"""
function GUItaucalculation(input_json)
    mu0 = 1.602176e-19 # 1 ev in Joules (J/eV) 
    bandtype, ϵ₀, μs, Ts, τ_form, τ_ac_pars, τ_im_pars, τ_matth_models, τ_matth_γ = GUIparseargs_tau(input_json) # parse the input json
    Ts = convert(Array{Float64},Ts)
    # check for empty tau
    if τ_form ∉ ("constant","acoustic","impurity","Matthiessen's rule")
        printstyled("ERROR: "; color = :red, bold=true)
        println("Relaxation time unknown.")
        println("Check input file.")
        return -20
    end  
    try
        τ = set_τ(τ_form, τ_ac_pars, τ_im_pars, τ_matth_models, τ_matth_γ) # set relaxation time func form
        local τ_fun
        if τ isa ScModel
            τ_fun = compute_τ
        elseif τ isa Matthiessen
            τ_fun = matthiessen_rule
        end
        τ_matrix = Array{Float64,2}(undef,length(Ts),length(μs))
        for t in eachindex(Ts)
            for m in eachindex(μs)
                τ_matrix[t,m] = τ_fun(τ,Ts[t],μs[m]*mu0,ϵ₀*mu0,bandtype)
            end
        end
        return Dict("mu" => μs, "data" => τ_matrix, "T" => Ts)  # send results to GUI
    catch e
        if isa(e, DomainError)
            printstyled("ERROR: "; color = :red, bold=true)
            println("Domain error in the τ function calculation.")
            println("Possible reasons: negative arguments in logs, negative radicand in radical expressions.")
            printstyled("Hint: "; color = :cyan, bold=true)
            println("shift the zero value of the chosen energy scale in the construction of the band structure.")
            return -30
        elseif isa(e, MatthiessenError)
            printstyled("ERROR: "; color = :red, bold=true)
            println("Please use the checkbox to select scatterings for Matthiessen's rule.")
            return -40
        else
            throw(e)
        end
    end
end