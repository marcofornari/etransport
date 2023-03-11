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


titlesize   = 22
xlabelsize  = 20
ylabelsize  = 20
zlabelsize  = 18
annsize     = 14

const mu0   = 1.602176e-19 # 1 ev in Joules (J/eV) 

jlblue = Colors.JULIA_LOGO_COLORS.blue
jlred = Colors.JULIA_LOGO_COLORS.red
jlgreen = Colors.JULIA_LOGO_COLORS.green


@doc raw"""
    plot(num_plots::Int64, x_axis::Union{Vector{Float64},Float64,Vector{Int64},Int64}, data::Union{Vector{Float64},Matrix{Float64},Vector{Vector{Float64}},Vector{Matrix{Float64}}}, z::Union{Vector{Float64},Vector{Int64}}=Float64[]; vlines::Union{Float64,Array{Float64}}=Float64[], annotations=[], color::ColorTypes.RGBA{Float64}=ColorTypes.RGBA{Float64}(Colors.JULIA_LOGO_COLORS.blue), colorscheme=:viridis, titles::Union{LaTeXString,Vector{LaTeXString}}=L"", xlabels::Union{LaTeXString,Vector{LaTeXString}}=L"", ylabels::Union{LaTeXString,Vector{LaTeXString}}=L"", zlabel::LaTeXString=L"")
 
Plot the transport coefficients. y must be a vector if plotting a single line, or a matrix of **shape** `length(x_axis)*length(z_axis)` if plotting more than one line. 

Parameters:
num_plots: number of transport coefficients to plot in the same figure
x: x axis vector
y: y axis vector/matrix
z: z axis vector
titles: titles for each plot
xlabels: x labels for each plot
ylabels: x labels for each plot
zlabel: z labels for the shared colorbar
color: line color
colorscheme: colormap
vlines: add vertical lines to the plot
annotations: add text annotation to the plot

Examples:
1. Single line (z is empty). Plot of electrical conductivity, Seebeck and carrier concentration as functions of the band gap.

titles  = [L"$σ$ vs gap, $μ$ = %$μ, $τ = const$", L"$S$ vs gap, $μ$ = %$μ, $τ = const$", L"$n$ vs gap, $μ$ = %$μ, $τ = const$"]
xlabels = [L"$Gap$ $[eV]$" for i in 1:3];
ylabels = [L"$\sigma$ $[(\Omega m)^{-1}]$", L"$S$ $[\mu VK^{-1}]$", L"n"];
annotations = Dict(L"v1" => [(0.11,6e5),(0.11,0),(0.11,5e6)], L"v2" => [(0.01,6e5-5e3),(0.01,0),(0.01,5e6-5e4)]);

fig = plot(3, gap, [σ,S,n], titles=titles, xlabels=xlabels, ylabels=ylabels, vlines=[0.0,0.1], annotations=annotations);

2. Multi-lines. Plot of electrical conductivity, Seebeck and carrier concentration as functions of the band gap and temperature (z axis). 

zlabel  = L"$T$ $[K]$";

fig = plot(3, gap, [σ,S,n], T, titles=titles, xlabels=xlabels, ylabels=ylabels, zlabel=zlabel, vlines=[0.0,0.1], annotations=annotations);

"""
function plot(num_plots::Int64, x_axis::Union{Vector{Float64},Float64,Vector{Int64},Int64}, data::Union{Vector{Float64},Matrix{Float64},Vector{Vector{Float64}},Vector{Matrix{Float64}}}, z::Union{Vector{Float64},Vector{Int64}}=Float64[]; vlines::Union{Float64,Array{Float64}}=Float64[], annotations=[], color::ColorTypes.RGBA{Float64}=ColorTypes.RGBA{Float64}(Colors.JULIA_LOGO_COLORS.blue), colorscheme=:viridis, titles::Union{LaTeXString,Vector{LaTeXString}}=L"", xlabels::Union{LaTeXString,Vector{LaTeXString}}=L"", ylabels::Union{LaTeXString,Vector{LaTeXString}}=L"", zlabel::LaTeXString=L"")
    
    fig = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98), resolution = (650*num_plots, 500))
    
    if isempty(z)   # data is a Vector or a Vector{Vector}
        if num_plots == 1   # just one plot
            ax = Axis(fig[1,1], xlabel=xlabels, ylabel=ylabels, xlabelsize=xlabelsize, ylabelsize=ylabelsize)
            Label(fig[1,1,Top()], titles, padding = (0, 0, 25, 10), fontsize=titlesize)
            lines!(ax, x_axis, vec(data), linewidth=1.5, color=color)
            for value in vlines
                vlines!(ax, value, linewidth=1, linestyle = :dash)
            end
            for ann in annotations
                text!(fig[1,1], ann.first, position = ann.second[1], color = :black, fontsize = annsize)
            end
        else
            for i in 1:num_plots
                ax = Axis(fig[1,i], xlabel=xlabels[i], ylabel=ylabels[i], xlabelsize=xlabelsize, ylabelsize=ylabelsize)
                Label(fig[1,i,Top()], titles[i], padding = (0, 0, 25, 10), fontsize=titlesize)
                lines!(ax, x_axis, vec(data[i]), linewidth=1.5, color=color)
                for value in vlines
                    vlines!(ax, value, linewidth=1, linestyle = :dash)
                end
                for ann in annotations
                    text!(fig[1,i], ann.first, position = ann.second[i], color = :black, fontsize = annsize)
                end
            end
        end

    else   # data is a Matrix or a Vector{Matrix}
        perm_idx = sortperm(z)
        z = z[perm_idx]
        palette = cgrad(colorscheme, length(z), categorical=true)
        if num_plots == 1   # just one plot
            c_i = 1
            ax = Axis(fig[1,1], xlabel=xlabels, ylabel=ylabels, xlabelsize=xlabelsize, ylabelsize=ylabelsize)
            Label(fig[1,1,Top()], titles, padding = (0, 0, 25, 10), fontsize=titlesize)
            if size(data,1) != length(x_axis)
                data = transpose(data)
            end
            data = data[:,perm_idx]
            for col in eachcol(data)
                lines!(ax, x_axis, col, linewidth=1.5, color=palette[c_i])
                c_i += 1
            end
            for value in vlines
                vlines!(ax, value, linewidth=1, linestyle = :dash)
            end
            for ann in annotations
                text!(fig[1,1], ann.first, position = ann.second[1], color = :black, fontsize = annsize)
            end
        else
            for i in 1:num_plots
                c_i = 1
                ax = Axis(fig[1,i], xlabel=xlabels[i], ylabel=ylabels[i], xlabelsize=xlabelsize, ylabelsize=ylabelsize)
                Label(fig[1,i,Top()], titles[i], padding = (0, 0, 25, 10), fontsize=titlesize)
                y_axis = nothing
                if size(data[i],1) != length(x_axis)
                    y_axis = transpose(data[i])
                else
                    y_axis = data[i]
                end
                y_axis = y_axis[:,perm_idx]
                for col in eachcol(y_axis)
                    lines!(ax, x_axis, col, linewidth=1.5, color=palette[c_i])
                    c_i += 1
                end
                for value in vlines
                    vlines!(ax, value, linewidth=1, linestyle = :dash)
                end
                for ann in annotations
                    text!(fig[1,i], ann.first, position = ann.second[i], color = :black, fontsize = annsize)
                end
            end
        end
        Colorbar(fig[1,num_plots+1], limits = (z[1],z[end]), colormap = colorscheme, label = zlabel, labelsize=zlabelsize)
    end

    return fig
end


# function to modify a fig returned by `plot` and include other plots. 
function plot!(fig::Figure, num_plots::Int64, x_axis::Union{Vector{Float64},Float64,Vector{Int64},Int64}, data::Union{Vector{Float64},Matrix{Float64},Vector{Vector{Float64}},Vector{Matrix{Float64}}}, z::Union{Vector{Float64},Vector{Int64}}=Float64[]; color::ColorTypes.RGBA{Float64}=ColorTypes.RGBA{Float64}(Colors.JULIA_LOGO_COLORS.blue), colorscheme=:viridis)
        
    if isempty(z)   # data is a Vector or a Vector{Vector}
        if num_plots == 1   # just one plot
            lines!(fig[1,1], x_axis, data[:,1], linewidth=1.5, color=color)
        else
            for i in 1:num_plots
                lines!(fig[1,i], x_axis, data[i][:,1], linewidth=1.5, color=color)
            end
        end

    else   # data is a Matrix or a Vector{Matrix}
        perm_idx = sortperm(z)
        z = z[perm_idx]
        palette = cgrad(colorscheme, length(z), categorical=true)
        if num_plots == 1   # just one plot
            c_i = 1
            data = data[perm_idx, :]
            for col in eachcol(data)
                lines!(fig[1,1], x_axis, col, linewidth=1.5, color=palette[c_i])
                c_i += 1
            end
        else
            for i in 1:num_plots
                c_i = 1
                y_axis = data[i][perm_idx, :]
                for col in eachcol(y_axis)
                    lines!(fig[1,i], x_axis, col, linewidth=1.5, color=palette[c_i])
                    c_i += 1
                end
            end
        end
    end
    return fig
end


# function to plot all the tensor components of a given transport coefficients as a function of temperature.
function tplot_allcomp(T::Union{Array{Float64},Array{Int64}}, data::DataFrame, name::String, z::Array{Float64}=Float64[]; issym::Bool, vlines::Union{Float64,Array{Float64}}=Float64[], annotations=[], color::ColorTypes.RGBA{Float64}=ColorTypes.RGBA{Float64}(Colors.JULIA_LOGO_COLORS.blue), ylabel::LaTeXString=L"", zlabel::LaTeXString=L"", param_string::String="")
    
    num_temp = length(T)
    titles = get_headerlabels(name, issym, param_string)

    # if symmetric plot only six components
    issym ? index_map = [5,4,2,3,0,1] : index_map = collect(8:-1:0) 
    issym ? rows = 2 : rows = 3

    if isempty(z)   # data is just a function of T
        fig = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98), resolution = (1600, 500*rows))
        for i in 1:3
            for j in 1:rows
                k = (i-1)*rows+j   # running index 1 -> 6
                ax = Axis(fig[j,i], xlabel=L"$T$ $[K]$", ylabel=ylabel, xlabelsize=xlabelsize, ylabelsize=ylabelsize)
                Label(fig[j,i,Top()], titles[k], padding = (0, 0, 25, 10), fontsize=titlesize)
                lines!(ax, T, data[num_temp,end-index_map[k]:1], linewidth=1.5, color=color)
                for value in vlines
                    vlines!(ax, value, linewidth=1, linestyle = :dash)
                end
                for ann in annotations
                    text!(fig[j,i], ann.first, position = ann.second[k], color = :black, fontsize = annsize)
                end
            end
        end
        return fig
    else
        fig = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98), resolution = (1800, 500*rows))
        num_cols = length(z)
        perm_idx = sortperm(z)
        z = z[perm_idx]
        palette = cgrad(:viridis, num_cols, categorical=true)
        for i in 1:3
            for j in 1:rows
                k = (i-1)*rows+j; c_i = 1  # running index 1 -> 6 and color index
                ax = Axis(fig[j,i], xlabel=L"$T$ $[K]$", ylabel=ylabel, xlabelsize=xlabelsize, ylabelsize=ylabelsize)
                Label(fig[j,i,Top()], titles[k], padding = (0, 0, 25, 10), fontsize=titlesize)

                y_axis = Matrix{Float64}(undef, num_cols, num_temp)   # construct the matrix of data 
                for t in eachindex(T)
                    y_axis[:,t] = data[(t-1)*num_cols+1:t*num_cols,end-index_map[k]];
                end
                y_axis = y_axis[perm_idx,:]

                for row in eachrow(y_axis)
                    lines!(ax, T, row, linewidth=1.5, color=palette[c_i])
                    c_i += 1
                end
                for value in vlines
                    vlines!(ax, value, linewidth=1, linestyle = :dash)
                end
                for ann in annotations
                    text!(fig[j,i], ann.first, position = ann.second[k], color = :black, fontsize = annsize)
                end
            end
        end
        Colorbar(fig[1:rows,4], limits = (z[1],z[end]), colormap = :viridis, label = zlabel, labelsize=zlabelsize)
        rowgap!(fig.layout, 5)
        return fig
    end
end


# methods to construct a vector of titles for tplot_allcomp function
function get_headerlabels(name::String, issym::Bool, param_string::String="")
    param_string == "" ? tail = param_string : tail = ", "*param_string
    indices = Vector{String}()
    if issym
        indices = ["[1,1]","[1,2]","[2,2]","[1,3]","[3,3]","[2,3]"]
    else
        indices = ["[1,1]","[2,1]","[3,1]","[1,2]","[2,2]","[3,2]","[1,3]","[2,3]","[3,3]"]
    end
    return [latexstring("\$"*name*index*"\$ vs \$T\$"*tail) for index in indices] 
end


@doc raw"""
    plot_bandstructure(bs::BandStructure,xaxis::AbstractArray=range(-1, 1, length=100); colors=nothing, label="")
 
Plot a given parabolic band structure. 

Parameters:
bs: band structure
xaxis: x axis vector
colors: tuple of three color from `Colors`. First: conduction bands. Second: valence band. Third: Fermi level.

Examples:
band_1 = ParabBand([5.0, 5.0, 5.0, 0.0, 0.0, 0.0],1.0, 1,1);    # conduction band
band_2 = ParabBand([0.1, 0.5, 3.0, 0.0, 0.0, 0.0],0.5,-1,1);    # valence band
μ = 0.8;
model = BandStructure(2,[band_1,band_2],μ);   # build the two-band structure
fig,ax = plot_bandstructure(model,colors=(:green,:red,:blue))
"""
function plot_bandstructure(bs::BandStructure,xaxis::AbstractArray=range(-1, 1, length=100); colors=nothing, label="")
    n_bands = bs.n
    bands = bs.bands
    μ = bs.μ

    # parabolic function to represent the bands
    f(type,m,en) = (type)*m*xaxis.^2 .+ en

    titles = [L"$m_x$", L"$m_y$", L"$m_z$"]

    fig = Figure(resolution = (500*3, 500))
    axes = [Axis(fig[1,i], title=titles[i], titlesize=titlesize) for i in 1:3]
    
    if colors === nothing    # if "color" not set by the user
        colors = (jlblue, jlred, jlgreen)
    end
   
    for (i,ax) in enumerate(axes)
        for b in 1:n_bands
            m = bands[b].mstar      # eff masses
            e0 = bands[b].ϵ₀        # energy min/max
            type = bands[b].type    # bandtype
            color = type == 1 ? colors[1] : colors[2]
            if b == 1
                lines!(ax, xaxis, f(type,m[i],e0), color=color, label=label)
            else
                lines!(ax, xaxis, f(type,m[i],e0), color=color)
            end
        end
        hlines!(ax, μ, color=colors[3], linestyle = :dash)	# fermi level
        hidedecorations!(ax, grid=false);
    end

    return fig,axes
end


@doc raw"""
    plot_τ(scm::Union{ScModel,Matthiessen}, t::Union{Vector{Float64},Vector{Int64}}, type="e-T",E_argmax::Int64=50; μ=nothing, ϵ₀::Float64=-42., bandtype::Int64=0) 
 
3D plot of a given relaxation time as a function of temperature, energy and Fermi level. 

Parameters:
scm: relaxation time model
t: temperature
type: type of plot specificed in "x_axis-z_axis" format. Available types are "e-T","T-e","μ-T","T".
μ: Fermi level vector. It's **mandatory** for type="μ-T".
ϵ₀: band energy. It's **mandatory** for acoustic relaxation time.
bandtype: band type (conduction, valence). It's **mandatory** for acoustic relaxation time.

Examples:
T = collect(300:10:650);

1. Constant relaxation time
τ_form = Scattering.constant()
fig = plot_τ(τ_form, T, "T");

2. Acoustic
τ_form = Scattering.acoustic(model,T₀=180,μ_min=5,μ_max=5);
fig = plot_τ(τ_form, T, "μ-T", μ=μ, ϵ₀=ϵ₀, bandtype=type);
"""
function plot_τ(scm::Union{ScModel,Matthiessen}, t::Union{Vector{Float64},Vector{Int64}}, type="e-T",E_argmax::Int64=50; μ=nothing, ϵ₀::Float64=-42., bandtype::Int64=0) 
    
    if type ∉ ("e-T","T-e","μ-T","T")
        error("[ERROR] Check type argument.")
    end

    # Acoustic τ requires band energy type
    if ( (scm isa Acoustic) || ( scm isa Matthiessen && any(x->x isa Acoustic, scm.τ_models)) ) && (ϵ₀ == -42. || bandtype == 0)
        error("[ERROR] Select bandtype and/or band position.")
    end

    fig = Figure(resolution = (1000, 500))
    
    local τ_fun
    if isa(scm,ScModel)
        τ_fun = compute_τ
    elseif isa(scm,Matthiessen)
        τ_fun = matthiessen_rule
    end

    local x
    if type =="μ-T"
        if length(μ) == 0
            error("[ERROR] Empty μ vector")
        end
        x = μ*mu0
    elseif (type == "T") && (μ !== nothing)
        x = [μ*mu0]
    else 
        x = (tint)[1:E_argmax]
    end

    y = convert(Array{Float64},t)

    if type != "T"
        τ = Matrix{Float64}(undef,length(x),length(y))   # precompute τ
        for i in eachindex(y)
            for j in eachindex(x)
                τ[j,i] = τ_fun(scm,y[i],convert(Float64,x[j]),ϵ₀*mu0,bandtype)
            end
        end

        # PLOT
        c_i = 1
        if type == "e-T"
            ax = Axis(fig[1,1], xlabel=L"E", ylabel=L"$τ$", xlabelsize=xlabelsize, ylabelsize=ylabelsize)
            Label(fig[1,1,Top()], L"$τ$ vs $ϵ,T$", padding = (0, 0, 25, 10), fontsize=titlesize)
            palette = cgrad(:viridis,length(y),categorical=true)
            for col in eachcol(τ)
                lines!(ax, x, col, linewidth=1.5, color=palette[c_i])
                c_i += 1
            end
            Colorbar(fig[1,2],limits=(y[1],y[end]),colormap=:viridis,label=L"T [K]", labelsize=zlabelsize)

        elseif type == "T-e"
            ax = Axis(fig[1,1], xlabel=L"T [K]", ylabel=L"$τ$", xlabelsize=xlabelsize, ylabelsize=ylabelsize)
            Label(fig[1,1,Top()], L"$τ$ vs $T,ϵ$", padding = (0, 0, 25, 10), fontsize=titlesize)
            palette = cgrad(:viridis,length(x),categorical=true)
            τ = transpose(τ)
            for col in eachcol(τ)
                lines!(ax, y, col, linewidth=1.5, color=palette[c_i])
                c_i += 1
            end
            Colorbar(fig[1,2],limits=(x[1],x[end]),colormap=:viridis,label=L"ϵ", labelsize=zlabelsize)

        elseif type == "μ-T"
            ax = Axis(fig[1,1], xlabel=L"μ [eV]", ylabel=L"$τ$", xlabelsize=xlabelsize, ylabelsize=ylabelsize)
            Label(fig[1,1,Top()], L"$τ$ vs $μ,T$", padding = (0, 0, 25, 10), fontsize=titlesize)
            palette = cgrad(:viridis,length(y),categorical=true)
            for col in eachcol(τ)
                lines!(ax, x, col, linewidth=1.5, color=palette[c_i])
                c_i += 1
            end
            Colorbar(fig[1,2],limits=(y[1],y[end]),colormap=:viridis,label=L"T [K]", labelsize=zlabelsize)
        end
        current_figure()
        return fig

    # type == "T"
    else
        ax = Axis(fig[1,1], xlabel=L"T [K]", ylabel=L"$τ$", xlabelsize=xlabelsize, ylabelsize=ylabelsize)
        Label(fig[1,1,Top()], L"$τ$ vs $T$", padding = (0, 0, 25, 10), fontsize=titlesize)
        τ = Vector{Float64}(undef,length(y))
        for i in eachindex(y)
            τ[i] = τ_fun(scm,y[i],x[1],ϵ₀*mu0,bandtype)
        end
        lines!(ax, y, τ, linewidth=1.5)
        current_figure()
        return fig
    end
end


@doc raw"""
    plot3d_τ(scm, t, type="e-T"; μ, ϵ₀, bandtype) 
 
3D plot of a given relaxation time as a function of temperature, energy and Fermi level. 

Parameters:
scm: relaxation time model
t: temperature
type: type of plot specificed in "x_axis-z_axis" format. Available types are "e-T","T-e","μ-T".
μ: Fermi level vector. It's **mandatory** for type="μ-T".
ϵ₀: band energy. It's **mandatory** for acoustic relaxation time.
bandtype: band type (conduction, valence). It's **mandatory** for acoustic relaxation time.

Examples:
T = collect(300:10:650);

1. Constant relaxation time
τ_form = Scattering.constant()
fig = plot3d_τ(τ_form, T, "T-e");

2. Acoustic
τ_form = Scattering.acoustic(model,T₀=180,μ_min=5,μ_max=5);
fig = plot3d_τ(τ_form, T, "μ-T", μ=μ, ϵ₀=ϵ₀, bandtype=type);
"""
function plot3d_τ(scm::Union{ScModel,Matthiessen},t::Union{Vector{Float64},Vector{Int64}},type::String="e-T",E_argmax::Int64=50; μ::Vector{Float64}=Float64[], ϵ₀::Float64=-42., bandtype::Int64=0) 

    if type ∉ ("e-T","T-e","μ-T")
        error("[ERROR] Check type argument.")
    end

    # Acoustic τ requires band energy type
    if ( (scm isa Acoustic) || ( scm isa Matthiessen && any(x->x isa Acoustic, scm.τ_models)) ) && (ϵ₀ == -42. || bandtype == 0)
        error("[ERROR] Select bandtype or band position.")
    end

    fig = Figure(resolution = (1600,1000))

    local τ_fun
    if scm isa ScModel
        τ_fun = compute_τ
    elseif scm isa Matthiessen
        τ_fun = matthiessen_rule
    end

    local x
    if type =="μ-T"
        if length(μ) == 0
            error("[ERROR] Empty μ vector")
        end
        x = μ*mu0
    else 
        x = (tint)[1:E_argmax]
    end

    y = convert(Array{Float64},t)
    τ = Matrix{Float64}(undef,length(x),length(y))
    # precompute τ
    for i in eachindex(y)
        for j in eachindex(x)
            τ[j,i] = τ_fun(scm,y[i],convert(Float64,x[j]),ϵ₀*mu0,bandtype)
        end
    end

    # plot the surface
    if type == "e-T"
        Label(fig[1,1,Top()], L"$τ$ vs $ϵ,T$", padding = (0, 0, 25, 10), fontsize=titlesize)
        ax = Axis3(fig[1,1],
                    xlabel = L"E",
                    ylabel = L"T [K]",
                    zlabel = L"$τ$",
                    xlabelsize=xlabelsize, 
                    ylabelsize=ylabelsize,
                    zlabelsize=zlabelsize
                    );
        hm = surface!(ax, x, y, τ, shading=false)
        wirx_idxs = [i for i in collect(1:length(x)) if i%6 == 1]
        wiry_idxs = [i for i in collect(1:length(y)) if i%3 == 1]
        wireframe!(ax, x[wirx_idxs], y[wiry_idxs], τ[wirx_idxs,wiry_idxs])
        contour3d!(ax, x, y, τ; levels=10)
        Colorbar(fig[1,2], hm)

    elseif type == "T-e"
        τ = transpose(τ)
        x,y = y,x 
        Label(fig[1,1,Top()], L"$τ$ vs $T,ϵ$", padding = (0, 0, 25, 10), fontsize=titlesize)
        ax = Axis3(fig[1,1],
        xlabel = L"T [K]",
        ylabel = L"E",
        zlabel = L"$τ$",
        xlabelsize=xlabelsize, 
        ylabelsize=ylabelsize,
        zlabelsize=zlabelsize
        );
        hm = surface!(ax, x, y, τ, shading=false)
        wiry_idxs = [i for i in collect(1:length(y)) if i%6 == 1]
        wireframe!(ax, x, y[wiry_idxs], τ[:,wiry_idxs])
        contour3d!(ax, x, y, τ; levels=10)
        Colorbar(fig[1,2], hm)

    elseif type == "μ-T"
        Label(fig[1,1,Top()], L"$τ$ vs $μ,T$", padding = (0, 0, 25, 10), fontsize=titlesize)
        ax = Axis3(fig[1,1],
                    xlabel = L"μ [eV]",
                    ylabel = "T [K]",
                    ylabelrotation = 4*pi,
                    zlabel = L"$τ$",
                    xlabelsize=xlabelsize, 
                    ylabelsize=ylabelsize,
                    zlabelsize=zlabelsize
                    );
        hm = surface!(ax, x/mu0, y, τ, shading=false)
        wirx_idxs = [i for i in collect(1:length(x)) if i%2 == 1]
        wireframe!(ax, x[wirx_idxs]/mu0, y, τ[wirx_idxs,:])
        contour3d!(ax, x/mu0, y, τ; levels=10)
        Colorbar(fig[1,2], hm)
    end
    current_figure()
    return fig
end


# interface function to export plots from CLI. 
function CLI_plotresults(path::String, tensors, it_e::Union{Vector{Array{Float64}},Vector{Any},Vector{Float64},Float64}, μs::Union{Vector{Float64}, Float64}, Ts::Union{Vector{Float64},Float64,Vector{Int64},Int64}; type::String="T")
    plot_path = joinpath(path, "plot")
    if !isdir(plot_path) # if there is not a folder to save the text data -> create it
        mkdir(plot_path)
    end

    counter_i = 4
    counter_j = 1

    local x_axis,xlabel,zlabel
    if type == "T"
        x_axis = Ts
        z_axis = μs
        xlabel = L"$T$ $[K]$"
        zlabel = L"$μ$ $[eV]$"
    elseif type == "μ"
        x_axis = μs
        z_axis = Ts
        xlabel = L"$μ$ $[eV]$"
        zlabel = L"$T$ $[K]$"
    end

    tensors_names = ["El_conductivity","Seebeck","Carrier_density","Ther_conductivity"]
    titles = [L"$σ$ vs $T$", L"$S$ vs $T$",L"$n$ vs $T$",L"$\kappa_{e}$ vs $T$"]
    ylabels = [L"$\sigma$ $[(\Omega m)^{-1}]$",L"$S$ $[\mu VK^{-1}]$",L"n",L"$\kappa_{e}$ $[WK^{-1}]$"]
    
    println("--> ", counter_i, ".0 Start plotting the results.")
    for (i,tensor) in enumerate(tensors)
        if !isempty(tensor)
            if tensors_names[i] == "Seebeck"
                tensor *= 1e6
            end
            print("---> ", counter_i, ".", counter_j, " ", tensors_names[i], " ...\r")
            filepath = joinpath(plot_path, (tensors_names[i]*".png"))
            if length(μs) == 1
                fig = plot(1, x_axis, to_matrix(tensor,length(it_e),length(μs),length(Ts)), titles=titles[i], xlabels=xlabel, ylabels=ylabels[i])
                savefig(filepath, fig);
            else
                fig = plot(1, x_axis, to_matrix(tensor,length(it_e),length(μs),length(Ts)), z_axis, titles=titles[i], xlabels=xlabel, ylabels=ylabels[i], zlabel=zlabel)
                savefig(filepath, fig);
            end
            println("---> ", counter_i, ".", counter_j, " ", tensors_names[i], " done.")
            counter_j += 1
        end
    end
    println("--> ", counter_i, ".", counter_j, " Plots saved.")
    counter_i += 1
end


@doc raw"""
    savefig(fullpath::String, fig::Figure)

Export to disk a given figure to `fullpath`.
"""
function savefig(fullpath::String, fig::Figure)
    save(fullpath, fig)
end
