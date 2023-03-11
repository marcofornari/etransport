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


function save_to_disk(export_all_data::Bool, path::String, tensors, cTensorParameters::Vector{Vector{Float64}}, it_e::Union{Vector{Array{Float64}},Vector{Any},Vector{Float64},Float64}, bandtypes::Vector{Int64}, degeneracies::Vector{Int64}, μs::Union{Vector{Float64}, Float64}, Ts::Union{Vector{Float64},Float64,Vector{Int64},Int64})

    data_path = string(path, "/data")
    if !isdir(data_path) # if there is not a folder to save the text data -> create it
        mkdir(data_path)
    end

    counter_i = 3
    counter_j = 1

    tensors_names = ["El_conductivity","Seebeck","Carrier_density","Ther_conductivity"]
    tensors_names_short = ["sigma","S","n","K"]
    
    println("--> ", counter_i, ".0 Start exporting all data to csv files.")
    for (tensor,tensor_name,tensor_name_short) in zip(tensors,tensors_names,tensors_names_short)
        if !isempty(tensor)
            local data
            print("---> ", counter_i, ".", counter_j, " ", tensor_name, " ...\r")
            filepath = string(data_path, "/", tensor_name, ".csv")
            if export_all_data && name != "n"  # export all tensor components (n is a scalar)
                data = fulltensor_df(tensor, tensor_name_short, cTensorParameters, it_e, bandtypes, degeneracies, μs, Ts)
            else  # export the tensor trace as a function of T
                data = to_df(tensor, tensor_name_short, cTensorParameters, it_e, bandtypes, degeneracies, μs, Ts)
            end
            CSV.write(filepath, data) 
            println("---> ", counter_i, ".", counter_j, " ", tensor_name, " done.")
            counter_j += 1
        end
    end
    println("--> ", counter_i, ".", counter_j, " Export all data completed.")
    counter_i += 1
end


function to_df(tensor::SharedArray{Float64}, name::String, cTensorParameters::Vector{Vector{Float64}}, it_e::Union{Vector{Array{Float64}},Vector{Any},Vector{Float64},Float64}, bandtypes::Vector{Int64}, degeneracies::Vector{Int64}, μs::Union{Vector{Float64}, Float64}, Ts::Union{Vector{Float64},Float64,Vector{Int64},Int64})

    # 1. define dataframe columns
    num_bands = length(bandtypes)
    model_header = get_bandheader(num_bands)
    
    num_energies = length(it_e)
    num_mu = length(μs)
    num_temps = length(Ts)
    
    temp_header = Vector{String}(undef, num_temps)
    for t in 1:num_temps
        temp_header[t] = *("T",string(t))
    end

    tensor_header = Vector{String}(undef, num_temps)
    for t in 1:num_temps
        tensor_header[t] = *(name,"_T",string(t))
    end

    header = cat(["Model #"], model_header, temp_header, tensor_header, dims=1)
    
    # 2. populate the dataset with the corresponding data
    data = Matrix{Float64}(undef, *(num_energies,num_mu), length(header))
    
    if name != "n"  # σ,S,K,L
        norm_const = 1/3
        for m in eachindex(μs)
            for (e,en) in enumerate(it_e)
                en /= mu0
                v_tensor_traces = Vector{Float64}(undef,num_temps)
                for t in eachindex(Ts)
                    v_tensor_traces[t] = norm_const*tr(tensor[:,:,t,e,m])
                end
                model_idx = (m-1)*num_energies+e
                data[model_idx,:] = cat(Int64(model_idx),vcat(cTensorParameters...),bandtypes,degeneracies,en,μs[m],Ts,v_tensor_traces, dims=1)
            end
        end
    else
        for m in eachindex(μs)
            for (e,en) in enumerate(it_e)
                en /= mu0
                v_tensor_traces = Vector{Float64}(undef,num_temps)
                for t in eachindex(Ts)
                    v_tensor_traces[t] = tensor[t,e,m]               
                end
                model_idx = (m-1)*num_energies+e
                data[model_idx,:] = cat(Int64(model_idx),vcat(cTensorParameters...),bandtypes,degeneracies,en,μs[m],Ts,v_tensor_traces, dims=1)
            end
        end
    end
    return DataFrame(data,header)
end


function to_matrix(tensor::Union{SharedArray{Float64},Array{Float64}},num_energies::Int64,num_mus::Int64,num_temps::Int64)
    
    data = Matrix{Float64}(undef, num_temps, *(num_energies,num_mus))

    if length(size(tensor)) == 5
        norm_const = 1/3
        for m in 1:num_mus
            for e in 1:num_energies
                model_idx = (m-1)*num_energies+e
                for t in 1:num_temps
                    data[t,model_idx] = norm_const*tr(tensor[:,:,t,e,m])
                end
            end
        end
    elseif length(size(tensor)) == 3
        for m in 1:num_mus
            for e in 1:num_energies
                model_idx = (m-1)*num_energies+e
                for t in 1:num_temps
                    data[t,model_idx] = tr(tensor[t,e,m])
                end
            end
        end
    end

    return data
end


function get_bandheader(num_bands::Int64)
    bands_header = Vector{String}(undef, 6*num_bands)
    bandtype_header = Vector{String}(undef, num_bands)
    deg_header = Vector{String}(undef, num_bands)
    energy_header = Vector{String}(undef, num_bands)
    for b in 1:num_bands
        bands_header[(b-1)*6+1:b*6] = .*(["mx_", "my_", "mz_","psi_", "theta_", "phi_"],string(b))
        bandtype_header[b] = *("Type_",string(b))
        deg_header[b] = *("Deg_",string(b))
        energy_header[b] = *("E0_",string(b))
    end
    return cat(bands_header, bandtype_header, deg_header, energy_header, ["mu"], dims=1)
end


function fulltensor_df(tensor::SharedArray{Float64}, name::String, cTensorParameters::Vector{Vector{Float64}}, it_e::Union{Vector{Array{Float64}},Vector{Any},Vector{Float64},Float64}, bandtypes::Vector{Int64}, degeneracies::Vector{Int64}, μs::Union{Vector{Float64}, Float64}, Ts::Union{Vector{Float64},Float64,Vector{Int64},Int64})
    
    # 1. define dataframe columns
    num_bands = length(bandtypes)
    model_header = get_bandheader(num_bands)

    num_energies = length(it_e)
    num_mu = length(μs)
    num_temps = length(Ts)

    # electrical and thermal conductivity are simmetric matrices
    if (name == "sigma") || (name == "K")
        tensor_header = .*(repeat([name], 6),["_1_1","_1_2","_1_3","_2_2","_2_3","_3_3"])
        header = cat(model_header, ["T"], tensor_header, dims=1)
        
        data = Matrix{Float64}(undef, *(num_energies,num_mu,num_temps), length(header))

        # add the dataframe rows
        for m in eachindex(μs)
            mu = μs[m]
            for (e,en) in enumerate(it_e)
                en /= mu0
                for t in eachindex(Ts)
                    low_triang = [tensor[:,:,t,e,m][i,j] for j in 1:3 for i in j:3]
                    data[((m-1)*num_energies+(e-1))*num_temps+t,:] = cat(vcat(cTensorParameters...),bandtypes,degeneracies,en,mu,Ts[t],low_triang, dims=1)
                end
            end
        end
    elseif (name == "S") || (name == "L")
        tensor_header = .*(repeat([name], 9),["_1_1","_2_1","_3_1","_1_2","_2_2","_3_2","_1_3","_2_3","_3_3"])
        header = cat(model_header, ["T"], tensor_header, dims=1)
        
        data = Matrix{Float64}(undef, *(num_energies,num_mu,num_temps), length(header))

        # add the dataframe rows
        for m in eachindex(μs)
            μ = μs[m]
            for (e,en) in enumerate(it_e)
                en /= mu0
                for t in eachindex(Ts)
                    data[((m-1)*num_energies+(e-1))*num_temps+t,:] = cat(vcat(cTensorParameters...),bandtypes,degeneracies,en,μ,Ts[t],vec(tensor[:,:,t,e,m]), dims=1)
                end
            end
        end
    end 
    return DataFrame(data,header)
end


@doc raw"""
    savedata(path::String, data::DataFrame)

Export to disk a given tranport coefficients simulation to a csv file.
"""
function savedata(path::String, data::DataFrame)
    CSV.write(path, data)
end

function savedata(path::String, data::Matrix{Float64})
    CSV.write(path, DataFrame(data,:auto))
end

