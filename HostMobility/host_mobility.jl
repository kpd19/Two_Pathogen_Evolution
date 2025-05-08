using Pkg
using DifferentialEquations
using CSV
using DataFrames
using Distributions
using LinearAlgebra
using Tables
using Dates
using Random

# running array jobs on the midway, idx is the index of the array
idx = parse(Int64, ENV["SLURM_ARRAY_TASK_ID"]) 

dist_df = DataFrame(CSV.File("files/coord_distances_meta13.csv"))
param_df = DataFrame(CSV.File("files/stoch/best_parameter_set.csv"))

meta = 3

dir = "output_sims/host_mobility/"

a = now()
println("Time started: "*string(a)*" hours")

const d = 0.2
const eps_val = 10 #0.1
const gen = 200 ### MAKE SURE THIS IS 200!!!!!
const index = 30 # 0.3 index

const bi = 1
const br = 1

const r = 0.2
const gamma = 0.2

const μ1 = 0.55
const μ2 = 0.55
const δ1 = 0.083
const δ2 = 0.083
const k1 = Int(20)
const k2 = Int(20)

const νSDO = 0.1
const νSGR = 0.1
const νMDO = 0.1
const νMGR = 0.1

const C_SNPV_DO = 3.7
const C_SNPV_GR = 3.7
const C_MNPV_DO = 2.5
const C_MNPV_GR = 3.3

# C_SNPV_DO = 3.50 # model 5
# C_SNPV_GR = 3.18 # model 5
# C_MNPV_DO = 2.48 # model 6
# C_MNPV_GR = 3.32 # model 6

# C_SNPV_DO = 3.38 # model 3 no tree species
# C_SNPV_GR = 3.38 # model 3 no tree species
# C_MNPV_DO = 2.65 # model 3 no tree species
# C_MNPV_GR = 2.65 # model 3 no tree species

# C_SNPV_DO = 3.1 # model 3 no tree species or morph
# C_SNPV_GR = 3.1 # model 3 no tree species or morph
# C_MNPV_DO = 3.1 # model 3 no tree species or morph
# C_MNPV_GR = 3.1 # model 3 no tree species or morph

global number_of_pops = length(Tables.matrix(unique(select(dist_df, "num"))))::Int
#number_of_pops = 1

global pop_nums = (1:number_of_pops)

function twostrain_SEIR(du,u,p,t)
    μ1, μ2, δ1, δ2, k1, k2, C1, C2, rho = p

    du[1] = -(u[4 + k1 + k2]*u[2]*u[1]) - (u[5 + k1 + k2]*u[3]*u[1])                     # S
    du[2] = -(C1^2*u[2]^2*u[4 + k1 + k2]) - (u[5 + k1 + k2]*rho*C1*C2*u[2]*u[3])            # nu_bar 1
    du[3] = -(C2^2*u[3]^2*u[5 + k1 + k2]) - (u[4 + k1 + k2]*rho*C1*C2*u[2]*u[3])            # nu_bar2

    du[4] = u[4 + k1 + k2]*u[2]*u[1] - k1*δ1*u[4]                                        # E1 first

    for i in 2:k1
        du[4 + i - 1] = (k1*δ1*u[4 + i - 2]) - k1*δ1*u[4 + i - 1]                      # E1 second
    end

    du[4 + k1] = u[5 + k1 + k2]*u[3]*u[1] - k2*δ2*u[4 + k1]                              # E2 first

    for i in 2:k2
        du[4 + k1 + i - 1] = k2*δ2*u[4 + k1 + i - 2] - k2*δ2*u[4 + k1 + i - 1]       # E2 second
    end

    du[4 + k1 + k2] = k1*δ1*u[4 + k1 - 1] - μ1*u[4 + k1 + k2]                        # P1
    du[5 + k1 + k2] = k2*δ2*u[4 + k1 + k2 - 1] - μ2*u[5 + k1 + k2]                   # P2
    du[6 + k1 + k2] = k1*δ1*u[4 + k1 - 1]                                            # psum1
    du[7 + k1 + k2] = k2*δ2*u[4 + k1 + k2 - 1]                                       # psum2
    nothing
end

function run_ode(S0, p0_1, p0_2, ν1, ν2, C1, C2, rho)
    rho = rho

    e1 = zeros(k1)
    e2 = zeros(k2)

    u0 = hcat(S0, ν1, ν2,transpose(e1),transpose(e2),p0_1, p0_2, 0,0)
    p = (μ1, μ2, δ1, δ2, k1, k2, C1, C2, rho)

    num_eq = length(u0)

    tspan = (0.0,70.0)
    prob = ODEProblem{true}(twostrain_SEIR,u0,tspan,p)
    sol = solve(prob, Tsit5(), abstol = 1e-8, reltol = 1e-8, isoutofdomain=(u,p,t) -> any(x -> x < 0, u), save_everystep = false, verbose = false)

    return(sol,num_eq)

end

function inter_annual_evolution(S_end, I1_end, I2_end, νi, νr, S_old, Z1old, Z2old, Ci, Cr, si, sr, ϕ1, ϕ2, ρ, stoch)

    frac_i1 = (I1_end)/(S_old)
    frac_i2 = (I2_end)/(S_old)
    Nt1 = stoch*S_end*(r + r*(si*νi + sr*νr))

    n1t1 = (1/(1 + si*νi + sr*νr))*(νi + si*(νi^2)*((bi^2)*(Ci^2) + 1) + sr*νr*νi*(ρ*bi*br*Ci*Cr + 1))
    n2t1 = (1/(1 + si*νi + sr*νr))*(νr + sr*(νr^2)*((br^2)*(Cr^2) + 1) + si*νr*νi*(ρ*bi*br*Ci*Cr + 1))

    Z1t1 = ϕ1*I1_end + gamma*Z1old
    Z2t1 = ϕ2*I2_end + gamma*Z2old

    ow_list = [frac_i1, frac_i2, Nt1, Z1t1, Z2t1, n1t1, n2t1]

    return ow_list
end

function place_tree_sp(n_do, n)
    if n_do >= 1
        tree_nums = sample(1:n,n_do, replace = false)
        tree_arr = repeat(["GR"],outer = 37)
        tree_arr[tree_nums] .= "DO"
    else
        tree_arr = repeat(["GR"],outer = 37)
    end
    return tree_arr
end

function get_prop(ϵ, distance_matrix)
    exp_mat = exp.(-ϵ*distance_matrix)
    exp_mat[diagind(exp_mat)] .= 0.00
    exp_row_sum = sum(exp_mat,dims = 2)
    prop_mat = exp_mat./exp_row_sum
    return(prop_mat)
end

function move_pops(pop_mat, t, prop_disp, d)
    B = Array{Float64}(undef, 1, number_of_pops)
    for i in 1:number_of_pops
        B[1,i] = pop_mat[i,t]*(1 - d) + sum(d.*pop_mat[:,t].*prop_disp[:,i])
    end

    return(B)
end

function move_pops_nu(S_mat, nu_mat, t, prop_disp, d)
    B = Array{Float64}(undef, 1, number_of_pops)
    for i in 1:number_of_pops
        B[1,i] = (S_mat[i,t]*nu_mat[i,t]*(1 - d) + sum(d.*S_mat[:,t].*prop_disp[:,i].*nu_mat[:,t]))/(S_mat[i,t]*(1-d) + sum(d.*S_mat[:,t].*prop_disp[:,i]))
    end
    return(B)
end

dist_mat = Array{Float64}(undef, number_of_pops, number_of_pops)
for i in pop_nums
    temp = select(filter(:ref=>==(i),dist_df), [:distance])
    dist_mat[:,i] = transpose(convert(Vector{Float64}, temp.distance))
end

if idx == 1
    p1 = get_prop(eps_val,dist_mat)::Matrix{Float64}
    CSV.write("prop_disp_"*string(meta)*".csv", Tables.table(p1),writeheader=false)
end

function dbinom(n,k,p)
    b = Binomial(k,p)
    return pdf(b,n)
end

function run_simulation(pdoug1::Int64, pdoug2::Int64, phiS::Float64, phiM::Float64, rep::Int64, tree_vals, sS_DO::Float64,sS_GR::Float64,sM_DO::Float64,sM_GR::Float64,rho::Float64,sigma::Float64)

    if number_of_pops >2
        prop_disp = get_prop(eps_val,dist_mat)::Matrix{Float64}
    elseif number_of_pops == 2
        prop_disp = [0 1; 1 0]::Matrix{Float64}
    else
        prop_disp = Array{Float64}(undef, 1, number_of_pops)
        prop_disp[1,1] = 0
    end

    S_pop = zeros(number_of_pops, gen+1)::Matrix{Float64}
    Z1_pop = zeros(number_of_pops, gen+1)::Matrix{Float64}
    Z2_pop = zeros(number_of_pops, gen+1)::Matrix{Float64}
    frac_pop1 = zeros(number_of_pops, gen+1)::Matrix{Float64}
    frac_pop2 = zeros(number_of_pops, gen+1)::Matrix{Float64}
    nu1_track = zeros(number_of_pops, gen+1)::Matrix{Float64}
    nu2_track = zeros(number_of_pops, gen+1)::Matrix{Float64}

    initS = rand(Uniform(0.1,2),number_of_pops)
    initMNPV = rand(LogNormal(-1,0.5),number_of_pops)
    initSNPV = rand(LogNormal(-1,0.5),number_of_pops)

    for p in 1:number_of_pops

        S_pop[p,1] = initS[p]
        Z1_pop[p,1] = initSNPV[p]
        Z2_pop[p,1] = initMNPV[p]

        if tree_vals[p] == "DO"
            ν_SNPV = νSDO # ν_SNPV_DO
            ν_MNPV = νMDO
        else
            ν_SNPV = νSGR
            ν_MNPV = νMGR
        end

        nu1_track[p,1] = ν_SNPV
        nu2_track[p,1] = ν_MNPV

    end

    for t in 2:(gen + 1)

        S_move = move_pops(S_pop,t-1, prop_disp, d)
        Z1_move = move_pops(Z1_pop,t-1, prop_disp, d)
        Z2_move = move_pops(Z2_pop,t-1, prop_disp, d)
        nu1_move = move_pops_nu(S_pop, nu1_track,t-1, prop_disp, d)
        nu2_move = move_pops_nu(S_pop, nu2_track,t-1, prop_disp, d)

        # S_pop_AD[:,t] .= transpose(S_move)
        # Z1_pop_AD[:,t] .= transpose(Z1_move)
        # Z2_pop_AD[:,t] .= transpose(Z2_move)
        # nu1_AD[:,t] .= transpose(nu1_move)
        # nu2_AD[:,t] .= transpose(nu2_move)

        en = rand(Normal(0,sigma),1)
        stoch = exp(en[1])

        for p in 1:number_of_pops

            if tree_vals[p] == "DO"
                C_SNPV = C_SNPV_DO
                C_MNPV = C_MNPV_DO

                sS = sS_DO
                sM = sM_DO
            else
                C_SNPV = C_SNPV_GR
                C_MNPV = C_MNPV_GR

                sS = sS_GR
                sM = sM_GR
            end

            S_init = S_move[p]::Float64
            Z1_init = Z1_move[p]::Float64 #+ 1e-8
            Z2_init = Z2_move[p]::Float64 #+ 1e-8
            nu1_init = nu1_move[p]::Float64
            nu2_init = nu2_move[p]::Float64

            output = run_ode(S_init,Z1_init,Z2_init,nu1_init,nu2_init,C_SNPV,C_MNPV,rho)

            simulation = output[1]
            num_equations = output[2]
            
            end_t = length(simulation.t)::Int64

            #num_equations = size(output)[2]

            S_end = simulation[1,end_t]::Float64
            nu1_end = simulation[2,end_t]::Float64
            nu2_end = simulation[3,end_t]::Float64
            I1_end = simulation[num_equations - 1,end_t]::Float64
            I2_end = simulation[num_equations,end_t]::Float64

            #S_end, I1_end, I2_end, νi, νr, Ci, Cr, bi, br, sr, si, r, ϕ1, ϕ2, γ, ρ
            over_winter = inter_annual_evolution(S_end, I1_end, I2_end, nu1_end, nu2_end, S_init, Z1_init, Z2_init, C_SNPV, C_MNPV, sS, sM, phiS, phiM, rho, stoch)::Vector{Float64}
            #over_winter = inter_annual(S_end, I1_end, I2_end, lambda, phi1, phi2)

            frac_pop1[p,t] = over_winter[1]
            frac_pop2[p,t] = over_winter[2]
            S_pop[p,t] = over_winter[3]
            Z1_pop[p,t] = over_winter[4]
            Z2_pop[p,t] = over_winter[5]
            nu1_track[p,t] = over_winter[6]
            nu2_track[p,t] = over_winter[7]

        end
    end

    frac1_df = DataFrame(frac_pop1, :auto)
    frac1_df.tree_sp = tree_vals
    frac1_df.pop = 1:74

    frac1_df = stack(frac1_df,1:201)

    rename!(frac1_df,:value => :FracI1)

    frac2_df = DataFrame(frac_pop2, :auto)
    frac2_df.pop = 1:74

    frac2_df = stack(frac2_df,1:201)

    rename!(frac2_df,:value => :FracI2)

    all = outerjoin(frac1_df,frac2_df, on = [:pop,:variable])          

    S_df = DataFrame(S_pop, :auto)
    S_df.pop = 1:74

    S_df = stack(S_df,1:201)

    rename!(S_df,:value => :S)

    all = outerjoin(all,S_df, on = [:pop,:variable])  
    
    Z1_df = DataFrame(Z1_pop, :auto)
    Z1_df.pop = 1:74

    Z1_df = stack(Z1_df,1:201)

    rename!(Z1_df,:value => :Z1)

    all = outerjoin(all,Z1_df, on = [:pop,:variable])  

    Z2_df = DataFrame(Z2_pop, :auto)
    Z2_df.pop = 1:74

    Z2_df = stack(Z2_df,1:201)

    rename!(Z2_df,:value => :Z2)

    all = outerjoin(all,Z2_df, on = [:pop,:variable])  

    nu1_df = DataFrame(nu1_track, :auto)
    nu1_df.pop = 1:74

    nu1_df = stack(nu1_df,1:201)

    rename!(nu1_df,:value => :nu1)

    all = outerjoin(all,nu1_df, on = [:pop,:variable])  

    nu2_df = DataFrame(nu2_track, :auto)
    nu2_df.pop = 1:74

    nu2_df = stack(nu2_df,1:201)

    rename!(nu2_df,:value => :nu2)

    all = outerjoin(all,nu2_df, on = [:pop,:variable])  

    all.pd1 .= pdoug1
    all.pd2 .= pdoug2
    all.rep .= rep

    return all

end

global sSD = param_df.sS_DO[1]
global sSG = param_df.sS_GR[1]

global sMD = param_df.sM_DO[1]
global sMG = param_df.sM_GR[1]

global phi_S = global phi_M = param_df.phi[1]
global rho = param_df.rho[1]

global sigma = 0.5 # param_df.sigma #[1e-6, 1e-2, 0.1, 0.5]

if isdir(dir) == false
    mkdir(dir) 
end

global ll_data = DataFrame()

for t1 in [2,35]

    for t2 in [2,35]

        for rp in 1:10

            tree_vals1 = place_tree_sp(t1,37)
            tree_vals2 = place_tree_sp(t2,37)

            tree_vals = [tree_vals1;tree_vals2]

            output_sim = run_simulation(t1,t2,phi_S,phi_S,rp,tree_vals,Float64(sSD),Float64(sSG),Float64(sMD),Float64(sMG),rho,sigma)
            global ll_data = vcat(ll_data,output_sim)

 
        end
    
    end
end

sim_info = "meta"*"_"*string(meta)*"_eps10.csv"
CSV.write(dir*sim_info, ll_data, header = true)

b = now()
time_elapsed = round((b-a).value/(60000*60),digits = 2)
println("Time elapsed: "*string(time_elapsed)*" hours")
