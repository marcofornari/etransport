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


module Utils


using JSON
using Dates
using LaTeXStrings
using SharedArrays
using LinearAlgebra
using DataFrames: DataFrame
using CSV

export  ParabBand,
        BandStructure,
        extractparams,
        combination_energy,
        parseargs,
        GUIparseargs,
        GUIparseargs_tau,
        create_testfolder,
        save_inputfile,
        create_runfolder,
        get_enlowestband,
        save_to_disk,
        to_df,
        to_matrix,
        fulltensor_df,
        savedata
        
include("../other/globals.jl")
include("../other/mtypes.jl")
include("../other/export.jl")
        

struct PathNotFound <: Exception
    var::String
end
Base.showerror(io::IO, e::PathNotFound) = print(io, e.var, " not found")

# this method is a mid-layer method to simplipy the interface of Mstar2t
function extractparams(model::BandStructure)
    num_bands = model.n
    cTensorParameters = Vector{Vector{Float64}}(undef,num_bands)
    ebandmins = Vector{Any}(undef,num_bands)
    bandtypes = Vector{Int64}(undef,num_bands)
    degeneracies = Vector{Int64}(undef,num_bands)
    for i in 1:num_bands
        cTensorParameters[i] = model.bands[i].mstar
        ebandmins[i] = model.bands[i].ϵ₀
        bandtypes[i] = model.bands[i].type
        degeneracies[i] = model.bands[i].deg
    end
    μ = model.μ
    return cTensorParameters,ebandmins,bandtypes,degeneracies,μ
end


# read ranges of band minima and construct an iterator over all the combinations of band structures
function combination_energy(num_bands, ebandmins)
    it_e = ebandmins[1]
    for i in 2:num_bands
        ebands = Array{Float64}[]
        for j in Iterators.product(it_e,ebandmins[i])
            push!(ebands,[k for k in Iterators.flatten(j)])
        end
        it_e = ebands  
    end
    return it_e
end


# function to parse input file from python client
function parseargs(input_json)
    args = input_json["args"]
    local_path = input_json["# results fullpath"]
    export_all_data = lowercase(input_json["# export all data [true/false]"]) == "true" ? true : false
    num_bands = input_json["# number of bands"]
    cTensorParameters = parse_bandparams(Float64, num_bands, input_json["# bands masses and angles"])
    ebandmins = custom_parse.([Float64], input_json["# energy extrema"])
    bandtypes = parse.([Int64], input_json["# band type"])
    degeneracies = parse.([Int64], input_json["# degeneracy"])
    μs = custom_parse(Float64, input_json["# Fermi level"])
    Ts = custom_parse(Float64, input_json["# temperature"])
    τ_form = input_json["# tau model [constant/acoustic/impurity]"]
    τ_ac_pars = input_json["# tau acoustic coefficients"]
    τ_im_pars = input_json["# tau impurity coefficients"]
    return args, local_path, export_all_data, num_bands, cTensorParameters, ebandmins, bandtypes, degeneracies, μs, Ts, τ_form, τ_ac_pars, τ_im_pars
end


# function to parse input file from GUI
function GUIparseargs(input_json)
    tensor_name = input_json["tensor_name"]
    export_all_data = input_json["# export all data [true/false]"]
    num_bands = input_json["# number of bands"]
    cTensorParameters = parse_bandparams(Float64, num_bands, input_json["# bands masses and angles"])
    ebandmins = custom_parse.([Float64], input_json["# energy extrema"])
    bandtypes = parse.([Int64], input_json["# band type"])
    degeneracies = parse.([Int64], input_json["# degeneracy"])
    μs = custom_parse(Float64, input_json["# Fermi level"])
    Ts = custom_parse(Float64, input_json["# temperature"])
    τ_form = input_json["# tau model [constant/acoustic/impurity/matthiessen]"]
    τ_ac_pars = input_json["# tau acoustic coefficients"]
    τ_im_pars = input_json["# tau impurity coefficients"]
    τ_matth_models = input_json["# tau matthiessen models"]
    τ_matth_γ = parse(Float64,input_json["# tau matthiessen gamma"])
    return tensor_name, export_all_data, num_bands, cTensorParameters, bandtypes, ebandmins, degeneracies, μs, Ts, τ_form, τ_ac_pars, τ_im_pars, τ_matth_models, τ_matth_γ
end


# function to parse relaxation time parameters from GUI
function GUIparseargs_tau(input_json)
    bandtype = input_json["# band type"]
    ebandmin = input_json["# energy extremum"]
    μs = custom_parse(Float64, input_json["# Fermi level"])
    Ts = custom_parse(Float64, input_json["# temperature"])
    τ_form = input_json["# tau model [constant/acoustic/impurity/matthiessen]"]
    τ_ac_pars = input_json["# tau acoustic coefficients"]
    τ_im_pars = input_json["# tau impurity coefficients"]
    τ_matth_models = input_json["# tau matthiessen models"]
    τ_matth_γ = parse(Float64,input_json["# tau matthiessen gamma"])
    return bandtype, ebandmin, μs, Ts, τ_form, τ_ac_pars, τ_im_pars, τ_matth_models, τ_matth_γ
end


# function that parse the parameter "# bands masses and angles" in the input file
function parse_bandparams(type, num_bands, bands)
    bandmasses = Vector{Array{Float64,1}}(undef,num_bands)
    for i in eachindex(bands)
        tmp = Array{Float64,1}(undef,6)
        single_band = filter!(e->e!="",split(bands[i]," "))
        for j in eachindex(single_band)
            single_param = single_band[j]
            tmp[j] = parse(type, single_param)
        end
        bandmasses[i] = tmp
    end
    return bandmasses
end


function custom_parse(type, value)
    # if input is a single value: call parse from standard library, if it is a range: call parse_range below
    return occursin(":", value) ? parse_range(value) : parse(type, value)
end


# return an evenly spaced array within a given range.
function parse_range(input_string)
    (start, stop, step) = [parse(Float64,match(r"[-+]?[0-9]*\.?[0-9]+", num).match) for num in split(input_string,":")]
    return collect(start:step:stop)
end


# function that creates folders for each test perfomed by the Python interface
# tests are organized as follows:
# results
#   └─── date_of_the_test
#           └─── time_of_the_test
function create_testfolder(script_path)
    test_path = string(script_path, "/results")
    if !isdir(test_path) # if there is not a results folder -> create it
        mkdir(test_path)
    end
    date = Dates.format(Dates.today(),"mm-dd-YYYY")
    test_path = string(test_path, "/", date)
    if !isdir(test_path) # if there is not a today folder -> create it
        mkdir(test_path)
    end
    test_name = Dates.format(now(), "HH-MM-SS")
    test_path = string(test_path, "/", test_name)
    mkdir(test_path)

    return test_path
end


# function that save an input file for debugging purposes
function save_inputfile(input_dict, path)
    open(path*"/input.json","w") do fout
        JSON.print(fout,input_dict) 
    end
end



# get the energy of the lowerst band in a band structure
function get_enlowestband(bandmodel::BandStructure)
    lowest_ϵ₀ = bandmodel.bands[1].ϵ₀
    @views bands = bandmodel.bands[2:end]
    for i in eachindex(bands)
        if bands[i].ϵ₀ < lowest_ϵ₀
            lowest_ϵ₀ = bands[i].ϵ₀
        end
    end
    return lowest_ϵ₀
end

end


