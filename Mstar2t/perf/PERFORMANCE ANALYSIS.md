## PERFORMANCE ANALYSIS

```bash
julia --project=.. parallel_perf_test.jl
```


### 1. v1
Number of calculation for each tensor: 2016

| Time [ms] | Time per it [μs] | allocations | allocations (MiB) |
| --------- | ---------------- | ----------- | ----------------- |
| 89.627    | 44               | 708564      | 27.96             |
| 186.043   | 92               | 1745397     | 60.68             |
| 363.245   | 180              | 5243344     | 131.88            |

### 2. v2

From REPL

Number of calculation for each tensor: 2016

| Time [ms] | Time per it [μs] | allocations | allocations (MiB) |
| --------- | ---------------- | ----------- | ----------------- |
| 87.913    | 43.6             | 682466      | 26.64             |
| 162.269   | 80.5             | 1168796     | 44.46             |
| 209.827   | 104.1            | 1669113     | 72.33             |

### 3. MKL

Number of calculation for each tensor: 2016

| Time [ms] | Time per it [μs] | allocations | allocations (MiB) |
| --------- | ---------------- | ----------- | ----------------- |
| 86.121    | 42.7             | 706551      | 28.45             |
| 174.550   | 86.6             | 1732847     | 62.25             |
| 336.379   | 166.9            | 4744556     | 123.05            |

### 4. MKL + PARALLEL 

Number of calculation for each tensor: 2016

| Time [ms] | Time per it [μs] | allocations | allocations (MiB) |
| --------- | ---------------- | ----------- | ----------------- |
| 17.680    | 8.77             | 680536      | 27.16             |
| 33.080    | 16.4             | 1164899     | 46.41             |
| 45.619    | 22.6             | 1665917     | 71.22             |


Test script

```julia
using Mstar2t
using BenchmarkTools
using Mstar2t: Scattering

function test()
    # BAND STRUCTURE DEFINITION
    m_1 = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0];
    ϵ₀_1 = 1.0;
    type_1 = 1;
    deg_1 = 1;
    band_1 = ParabBand(m_1,ϵ₀_1,type_1,deg_1);   # create the conduction band
    
    m_2 = [0.1, 1.0, 5.0, 0.0, 0.0, 0.0];
    ϵ₀_2 = 0.0;
    type_2 = -1;
    deg_2 = 1;
    band_2 = ParabBand(m_2,ϵ₀_2,type_2,deg_2);   # create the valence band

    μ = collect(-.1:0.01:.1);
    model = BandStructure(2,[band_1,band_2],μ);   # build the two-band structure

    T = collect(50.:10:750);

    τ_form = Scattering.constant();

    num_calc = length(μ)*length(T)
    println("Number of calculation for each tensor: ", num_calc)

    # TENSORS COMPUTATION
    @btime electrical_conductivity($model,$T,$τ_form);
    @btime seebeck_coefficient($model,$T,$τ_form);
    @btime carrier_concentration($model,$T,$τ_form);
end;

test();
```

