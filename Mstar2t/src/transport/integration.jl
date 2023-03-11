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


### SERIES INTEGRATION PART ###
"""
    k_series_int(n::Int64, m::Real, s::Real))
 
This function handles all the cases of the series integration. `n` is the index for the kinetic coefficient, `m` for the coefficient of the τ expansion.
"""
function k_series_int(n::Int64, m::Real, sr::Real)
    s = complex(sr) # Since a lot of the following math makes use of complex numbers
    f(x) = k_series(n, m, x, s) # Adds an alias for integration, as the numerical integrator doesn't like functions of multiple variables.
    var = 0.0 # The output
    if n == 0
        if m == -3/2
            sb = BigFloat(s) # We need additional numerical precision here
            var = 1/(exp(sb) + 1)
        elseif m == -1
            if sr >= 10
                var = 0.5*sqrt(pi)*exp(-s)
            else
                var = 0.5*sqrt(pi)*exp(-s)*SpecialFunctions.erfc(sqrt(10-s)) + sqrt(10-s)/exp(10)
                var += QuadGK.quadgk(f, sr < -20 ? -20 : s, 10, rtol=1.0e-10, atol=1.0e-18, maxevals=10^10)[1]
                if sr < -20
                    var += 0.5im*exp(s)*(sqrt(pi)-2custom_inc_gamma(3/2, sr+20))
                end
            end
        elseif m == -1/2
            sb = BigFloat(s) # We need additional numerical precision here
            var = log(exp(-sb) + 1) # I had to rewrite this to reduce numerical instability
        elseif m == 0
            if sr >= 10
                var = 0.75*sqrt(pi)*exp(-s)
            else
                var = 0.75*sqrt(pi)*exp(-s)*SpecialFunctions.erfc(sqrt(10-s)) + sqrt(10-s)*(23-2s)/(2*exp(10))
                var += QuadGK.quadgk(f, sr < -20 ? -20 : s, 10, rtol=1.0e-10, atol=1.0e-18, maxevals=10^10)[1]
                if sr < -20
                    var += -0.25im*exp(s)*(3sqrt(pi)-4custom_inc_gamma(5/2, sr+20))
                end
            end
        elseif m == 1/2
            if sr >= 10
                sb = BigFloat(s) # We need additional numerical precision here
                var = 2*exp(-sb)
            else
                var = ((s-22)*s + 122)/exp(10)
                var += QuadGK.quadgk(f, sr < -20 ? -20 : s, 10, rtol=1.0e-10, atol=1.0e-18, maxevals=10^10)[1]
                if sr < -20
                    var += (s*(s+42)+442)/exp(20) - 2exp(s)
                end
            end
        elseif m == 1
            if sr >= 10
                var = (15/8)*sqrt(pi)*exp(-s)
            else
                var = (15/8)*sqrt(pi)*exp(-s)*SpecialFunctions.erfc(sqrt(10-s)) + sqrt(10-s)*(4s^2-90s+515)/(4*exp(10))
                var += QuadGK.quadgk(f, sr < -20 ? -20 : s, 10, rtol=1.0e-10, atol=1.0e-18, maxevals=10^10)[1]
                if sr < -20
                    var += (1/8)*im*exp(s)*(15sqrt(pi)-8custom_inc_gamma(7/2, sr+20))
                end
            end
        elseif m == 3/2
            if sr >= 10
                var = 6*exp(-s)
            else
                var = (1366 - s*((s-33)*s + 366))/exp(10)
                var += QuadGK.quadgk(f, sr < -20 ? -20 : s, 10, rtol=1.0e-10, atol=1.0e-18, maxevals=10^10)[1]
                if sr < -20
                    var += 6exp(s) -(s*(s*(s+63) + 1326) + 9326)/exp(20)
                end
            end
        elseif m == 2
            if sr >= 10
                var = (105/16)*sqrt(pi)*exp(-s)
            else
                var = (105/16)*sqrt(pi)*exp(-s)*SpecialFunctions.erfc(sqrt(10-s)) + sqrt(10-s)*(11605-2s*(2s*(2s-67) + 1515))/(8*exp(10))
                var += QuadGK.quadgk(f, sr < -20 ? -20 : s, 10, rtol=1.0e-10, atol=1.0e-18, maxevals=10^10)[1]
                if sr < -20
                    var += -im*exp(s)*(105*sqrt(pi)/16 - custom_inc_gamma(9/2, sr+20))
                end
            end
        end
    elseif n == 1
        if m == -3/2
            # I had to modify this function to make it numerically stable.
            sb = BigFloat(s) # We need additional numerical precision here
            var = log(exp(sb) + 1) - s/(exp(-sb) + 1)
        elseif m == -1
            if sr >= 10
                var = 0.25*sqrt(pi)*exp(-s)*(2s+3)
            else
                var = 0.25*(sqrt(pi)*exp(-s)*(2s+3)*SpecialFunctions.erfc(sqrt(10-s)) + 46*sqrt(10-s)/exp(10))
                var += QuadGK.quadgk(f, sr < -20 ? -20 : s, 10, rtol=1.0e-10, atol=1.0e-18, maxevals=10^10)[1]
                if sr < -20
                    var += ((6-4s)*SpecialFunctions.dawson(sqrt(-s-20))-86sqrt(-s-20))/4exp(20)
                end
            end
        elseif m == -1/2
            sb = BigFloat(s) # We need additional numerical precision here
            if sr >= 10
                var = (-1)^(m+0.5)*SpecialFunctions.gamma(m+5/2)*SpecialFunctions.gamma(n+1)*sb^(m+n+5/2)*HypergeometricFunctions.drummond1F1(n+1, m+n+7/2, -sb)
                var /= SpecialFunctions.gamma(m+n+7/2)
                var += HypergeometricFunctions.drummond1F1(-m-3/2, -m-n-3/2, -sb)*SpecialFunctions.gamma(m+n+5/2)
            else
                var = sb*log(exp(sb) + 1) + pi^2/3 + 2*PolyLog.li2(ComplexF64(-exp(sb)))
            end
        elseif m == 0
            if sr >= 10
                var = (3/8)*sqrt(pi)*exp(-s)*(2s+5)
            else
                var = (3/8)*sqrt(pi)*exp(-s)*(2s+5)*SpecialFunctions.erfc(sqrt(10-s)) + (515-44s)*sqrt(10-s)/(4*exp(10))
                var += QuadGK.quadgk(f, sr < -20 ? -20 : s, 10, rtol=1.0e-10, atol=1.0e-18, maxevals=10^10)[1]
                if sr < -20
                    var += 3*((2s-5)*SpecialFunctions.dawson(sqrt(-s-20))+sqrt(-s-20)*(28s+605))/4exp(20)
                end
            end
        elseif m == 1/2
            sb = BigFloat(s) # We need additional numerical precision here
            if sr >= 10
                var = (-1)^(m+0.5)*SpecialFunctions.gamma(m+5/2)*SpecialFunctions.gamma(n+1)*sb^(m+n+5/2)*HypergeometricFunctions.drummond1F1(n+1, m+n+7/2, -sb)
                var /= SpecialFunctions.gamma(m+n+7/2)
                var += HypergeometricFunctions.drummond1F1(-m-3/2, -m-n-3/2, -sb)*SpecialFunctions.gamma(m+n+5/2)
            else
                var = 2*sb*PolyLog.li2(ComplexF64(-exp(sb))) - 6*PolyLog.li3(ComplexF64(-exp(sb))) - 2*pi^2*sb/3
            end
        elseif m == 1
            if sr >= 10
                var = (15/16)*sqrt(pi)*exp(-s)*(2s+7)
            else
                var = (15/16)*sqrt(pi)*exp(-s)*(2s+7)*SpecialFunctions.erfc(sqrt(10-s)) + (8s*(11s-250)+11605)*sqrt(10-s)/(8*exp(10))
                var += QuadGK.quadgk(f, sr < -20 ? -20 : s, 10, rtol=1.0e-10, atol=1.0e-18, maxevals=10^10)[1]
                if sr < -20
                    var += (15*(7-2s)*SpecialFunctions.dawson(sqrt(-s-20)) - sqrt(-s-20)*(8s*(21s+895)+76705))/8exp(20)
                end
            end
        elseif m == 3/2
            if sr >= 10
                var = 6exp(-s)*(s+4)
            else
                var = (15464+s*(s*(366-11s)-4098))/exp(10)
                var += QuadGK.quadgk(f, sr < -20 ? -20 : s, 10, rtol=1.0e-10, atol=1.0e-18, maxevals=10^10)[1]
                if sr < -20
                    var += 3*(2exp(s+20)*(s-4) + s*(s*(7s+442)+9326)+65768)/exp(20)
                end
            end
        elseif m == 2
            if sr >= 10
                var = (105/32)*sqrt(pi)*exp(-s)*(2s+9)
            else
                var = (1/32)*(105*sqrt(pi)*exp(-s)*(2s+9)*SpecialFunctions.erfc(sqrt(10-s)) + 2*(264445-4s*(4s*(11s-372)+17015))*sqrt(10-s)/exp(10))
                var += QuadGK.quadgk(f, sr < -20 ? -20 : s, 10, rtol=1.0e-10, atol=1.0e-18, maxevals=10^10)[1]
                if sr < -20
                    var += 7*(15*(2s-9)*SpecialFunctions.dawson(sqrt(-s-20))-sqrt(-s-20)*(4s*(4s*(3s+191)+16275)+464335))/16exp(20)
                end
            end
        end
    elseif n == 2
        if m == -3/2 
            sb = BigFloat(s)
            if sr >= 10
                var = (-1)^(m+0.5)*SpecialFunctions.gamma(m+5/2)*SpecialFunctions.gamma(n+1)*sb^(m+n+5/2)*HypergeometricFunctions.drummond1F1(n+1, m+n+7/2, -sb)
                var /= SpecialFunctions.gamma(m+n+7/2)
                var += HypergeometricFunctions.drummond1F1(-m-3/2, -m-n-3/2, -sb)*SpecialFunctions.gamma(m+n+5/2)
            else
                var = 2*PolyLog.li2(ComplexF64(-exp(sb))) - exp(sb)*sb^2/(exp(sb)+1) + 2sb*log(exp(sb)+1) + pi^2/3
            end
        elseif m == -1
            if sr >= 10
                var = (1/8)*sqrt(pi)*exp(-s)*(4s*(s+3)+15)
            else
                var = (1/8)*(sqrt(pi)*exp(-s)*(4s*(s+3)+15)*SpecialFunctions.erfc(sqrt(10-s))+2*sqrt(10-s)*(2s+515)/exp(10))
                var += QuadGK.quadgk(f, sr < -20 ? -20 : s, 10, rtol=1.0e-10, atol=1.0e-18, maxevals=10^10)[1]
                if sr < -20
                    var += ((-4*(s-3)*s-15)*SpecialFunctions.dawson(sqrt(-s-20))+sqrt(-s-20)*(1815-2s))/4exp(20)
                end
            end
        elseif m == -1/2
            sb = BigFloat(s)
            if sr >= 10
                var = (-1)^(m+0.5)*SpecialFunctions.gamma(m+5/2)*SpecialFunctions.gamma(n+1)*sb^(m+n+5/2)*HypergeometricFunctions.drummond1F1(n+1, m+n+7/2, -sb)
                var /= SpecialFunctions.gamma(m+n+7/2)
                var += HypergeometricFunctions.drummond1F1(-m-3/2, -m-n-3/2, -sb)*SpecialFunctions.gamma(m+n+5/2)
            else
                var = 4*sb*PolyLog.li2(ComplexF64(-exp(sb))) - 6*PolyLog.li3(ComplexF64(-exp(sb))) + s^2*log(exp(sb) + 1) - pi^2*sb/3
            end
        elseif m == 0
            if sr >= 10
                var = (3/16)*sqrt(pi)*exp(-s)*(4s*(s+5)+35)
            else
                var = (3/16)*sqrt(pi)*exp(-s)*(4s*(s+5)+35)*SpecialFunctions.erfc(sqrt(10-s))+5*sqrt(10-s)*(2321-194s)/(8*exp(10))
                var += QuadGK.quadgk(f, sr < -20 ? -20 : s, 10, rtol=1.0e-10, atol=1.0e-18, maxevals=10^10)[1]
                if sr < -20
                    var += (3*(4*(s-5)*s+35)*SpecialFunctions.dawson(sqrt(-s-20))-5sqrt(-s-20)*(706s+15341))/8exp(20)
                end
            end
        elseif m == 1/2
            if sr >= 10
                var = exp(-s)*(2s*(s+6)+24)
            else
                var = 2*(s*(61s-1366)+7732)/exp(10)
                var += QuadGK.quadgk(f, sr < -20 ? -20 : s, 10, rtol=1.0e-10, atol=1.0e-18, maxevals=10^10)[1]
                if sr < -20
                    var += 2*(s*(221s+9326)-exp(s+20)*((s-6)*s+12)+98652)/exp(20)
                end
            end
        elseif m == 1
            if sr >= 10
                var = (15/32)*sqrt(pi)*exp(-s)*(4s*(s+7)+63)
            else
                var = (1/32)*(15*sqrt(pi)*exp(-s)*(4s*(s+7)+63)*SpecialFunctions.erfc(sqrt(10-s)) + 2*sqrt(10-s)*(2s*(976s-22425)+264445)/exp(10))
                var += QuadGK.quadgk(f, sr < -20 ? -20 : s, 10, rtol=1.0e-10, atol=1.0e-18, maxevals=10^10)[1]
                if sr < -20
                    var += (sqrt(-s-20)*(7072s^2+302290s+3250345)-15*(4*(s-7)*s+63)*SpecialFunctions.dawson(sqrt(-s-20)))/16exp(20)
                end
            end
        elseif m == 3/2
            if sr >= 10
                var = 6*exp(-s)*(s*(s+8)+20)
            else
                var = (177320 - 2s*(s*(61s-2049)+23196))/exp(10)
                var += QuadGK.quadgk(f, sr < -20 ? -20 : s, 10, rtol=1.0e-10, atol=1.0e-18, maxevals=10^10)[1]
                if sr < -20
                    var += 6exp(s)*((s-8)*s + 20) - 2*(s*(s*(221s+13989)+295956)+2093260)/exp(20)
                end
            end
        elseif m == 2
            if sr >= 10
                var = (105/64)*sqrt(pi)*exp(-s)*(4s*(s+9)+99)
            else
                sc = sr < -30 ? -30 : s
                var = (1/16)*exp(-sc-10)*(exp(10)*(16*(2sc*custom_inc_gamma(11/2, 10-real(sc)) + custom_inc_gamma(13/2, 10-real(sc)))-105*sqrt(pi)*sc^2*erf(sqrt(10-sc))) - sc^2*(2*exp(sc)*sqrt(10-sc)*(2sc*(2sc*(2sc-67)+1515)-11605)-105*exp(10)*sqrt(pi)))
                var += QuadGK.quadgk(f, sr < -20 ? -20 : s, 10, rtol=1.0e-10, atol=1.0e-18, maxevals=10^10)[1]
                if sr < -20
                    var += (-2*sqrt(-20-s)*(138153795+2s*(9649415+16s*(28221+442s))) + 210*(99+4*(-9+s)*s)*SpecialFunctions.dawson(sqrt(-20-sr)))/(64*exp(20))
                end
            end
        end
    end
    return real(var)
end

# Integrand
function k_series(n, m, x, s)
    return x^n*(x-s)^(1.5 + m)*exp(x)/((1+exp(x))^2)
end

# The Julia implementation of the incomplete gamma function is insufficient for our needs, so here is a custom one
function custom_inc_gamma(a, x)
    if x >= 0
        return SpecialFunctions.gamma_inc(a, x, 0)[2]*SpecialFunctions.gamma(a)
    else
        return SpecialFunctions.gamma(a) + QuadGK.quadgk(y->y^(a-1)*exp(-y), complex(x), 0+0im, rtol=1.0e-10, atol=1.0e-18, maxevals=10^10)[1]
    end
end


### NUMERICAL INTEGRATION PART ###

tint, wint = FastGaussQuadrature.gausslaguerre(Nint)    # Gaussian quadrature nodes computation 

# f(x) in the Gauss-Laguerre quadrature (Ref.: https://en.wikipedia.org/wiki/Gauss%E2%80%93Laguerre_quadrature)
function kn_arb_int(t::Float64, s::Float64, n::Int64)
    if t < 3e2 # Protects against numerical overflow
        return (t+s)^n*exp(t)*t^1.5*(exp(t+s)/(1+exp(t+s))^2)
    else
        return exp(-s)*(t+s)^n*t^1.5
    end
end