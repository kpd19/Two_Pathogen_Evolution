using Pkg
using DifferentialEquations
using CSV
using DataFrames
using LinearAlgebra
using Tables
using Dates
using Random

post_df = DataFrame(CSV.File("files/thinned_posteriors.csv"))

const gamma = 0.2

const μ = 0.55
const δ = 0.083
const k = Int(20)

function onestrain_SEIR(du,u,p,t)
    ν, δ, k, C, S0, μ = p

    du[1] = -ν*u[1]*u[2 + k]*((u[1]/S0)^(C^2)) # S
    du[2] = ν*u[1]*u[2 + k]*((u[1]/S0)^(C^2)) - k*δ*u[2]

    for i in 2:k
        du[2 + i - 1] = (k*δ*u[2 + i - 2]) - k*δ*u[2 + i - 1]                      # E1 second
    end

    du[2 + k] = k*δ*u[2 + k - 1] - μ*u[2 + k]                                       # P1
    du[3 + k] = k*δ*u[2 + k - 1]                                                    # psum1
    nothing
end

function run_ode(S0, p0,ν,C)
    e = zeros(k)

    u0 = hcat(S0, transpose(e),p0, 0.0)
    p = (ν, δ, k, C, S0, μ)

    num_eq = length(u0)

    tspan = (0.0,70.0)
    prob = ODEProblem{true}(onestrain_SEIR,u0,tspan,p)
    sol = solve(prob, Tsit5(), abstol = 1e-8, reltol = 1e-8, isoutofdomain=(u,p,t) -> any(x -> x < 0, u), save_everystep = false, verbose = false)

    return(sol,num_eq)

end

function run_simulation(ν, C, morph, tree, draw)

    initS = 1:1:1000
    initZ = 0.1

    frac_I = zeros(length(initS))

    for i in 1:length(initS)
        output = run_ode(initS[i],initZ,ν,C)

        simulation = output[1]
        num_equations = output[2]
            
        end_t = length(simulation.t)::Int64
    
        S_end = simulation[1,end_t]::Float64
        I_end = simulation[num_equations,end_t]::Float64
    
        frac_I[i] = I_end/initS[i]

    end

    dens_df = DataFrame(host_dens = initS, frac_I = frac_I)
    dens_df.morph .= morph

    dens_df.tree .= tree
    dens_df.draw .= draw
    dens_df.C .= C
    dens_df.nu_bar .= ν
    
    return dens_df

end

global nu_DO = post_df.nu_DO
global nu_GR = post_df.nu_GR

global C_DO = post_df.C_DO
global C_GR = post_df.C_GR

global morph_id = post_df.morph
global draw_id = post_df.draw

dir = "output/"

if isdir(dir) == false
    mkdir(dir) 
end

global all_dens = DataFrame()

for t in ["DO","GR"]

    if t == "DO"
        C_arr = C_DO
        nu_arr = nu_DO
    else 
        C_arr = C_GR
        nu_arr = nu_GR
    end

    for j in 1:length(C_arr)

        output_sim = run_simulation(nu_arr[j], C_arr[j], morph_id[j],t,draw_id[j])
    
        global all_dens = vcat(all_dens,output_sim)

    end

end
    

CSV.write(dir*"posteriors_df.csv", all_dens, header = true)
