using Pkg
using DifferentialEquations
using CSV
using DataFrames
using Distributions
using LinearAlgebra
using Tables
using Dates
using MPI
#using Random

MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)

global num_reps = size

global pdoug = 18
global num_sims = 6

model  = "no trees"
dir = "linesearch_op/op_nt1"

if rank == 0
    a = now()
    println("Time started: "*string(a)*" hours")
    println("Size is $size")
    
end

dist_df = DataFrame(CSV.File("files/coord_distances_R3.csv"))
ll_df = DataFrame(CSV.File("files/morphotype_dist_data.csv"))
ll_df[:,:n_trees] = convert.(Int,round.(ll_df.n_trees, digits = 0))

const d = 0.2
const eps_val = 0.1
const gen = 200
const index = 30

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

if model == "no trees"
    const C_SNPV_DO = 3.4 # model 3 no tree species
    const C_SNPV_GR = 3.4 # model 3 no tree species
    const C_MNPV_DO = 2.9 # model 3 no tree species
    const C_MNPV_GR = 2.9 # model 3 no tree species
elseif model == "no morphotype"
    const C_SNPV_DO = 3.1 # model 3 no tree species or morph
    const C_SNPV_GR = 3.1 # model 3 no tree species or morph
    const C_MNPV_DO = 3.1 # model 3 no tree species or morph
    const C_MNPV_GR = 3.1 # model 3 no tree species or morph
else 
    print("Error: No model specified")
end

global number_of_pops = length(Tables.matrix(unique(select(dist_df, "num"))))::Int

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

function dbinom(n,k,p)
    b = Binomial(k,p)
    return pdf(b,n)
end

function run_simulation(pdoug::Int64, phiS::Float64, phiM::Float64, rep::Int64, tree_vals, si::Float64,sr::Float64,rho::Float64,sigma::Float64)

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

    initS = rand(Uniform(0.4,1.8),1)[1]::Float64
    initMNPV = rand(LogNormal(-1,0.5),1)[1]::Float64
    initSNPV = rand(LogNormal(-1,0.5),1)[1]::Float64

    for p in 1:number_of_pops

        S_pop[p,1] = initS
        Z1_pop[p,1] = initSNPV
        Z2_pop[p,1] = initMNPV

        if tree_vals[p] == "DO"
            ν_SNPV = νSDO 
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

        en = rand(Normal(0,sigma),1)
        stoch = exp(en[1])

        for p in 1:number_of_pops

            if tree_vals[p] == "DO"
                C_SNPV = C_SNPV_DO
                C_MNPV = C_MNPV_DO
            else
                C_SNPV = C_SNPV_GR
                C_MNPV = C_MNPV_GR
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

            S_end = simulation[1,end_t]::Float64
            nu1_end = simulation[2,end_t]::Float64
            nu2_end = simulation[3,end_t]::Float64
            I1_end = simulation[num_equations - 1,end_t]::Float64
            I2_end = simulation[num_equations,end_t]::Float64

            over_winter = inter_annual_evolution(S_end, I1_end, I2_end, nu1_end, nu2_end, S_init, Z1_init, Z2_init, C_SNPV, C_MNPV, si, sr, phiS, phiM, rho, stoch)::Vector{Float64}

            frac_pop1[p,t] = over_winter[1]
            frac_pop2[p,t] = over_winter[2]
            S_pop[p,t] = over_winter[3]
            Z1_pop[p,t] = over_winter[4]
            Z2_pop[p,t] = over_winter[5]
            nu1_track[p,t] = over_winter[6]
            nu2_track[p,t] = over_winter[7]

        end
    end

    total_frac = (frac_pop1 + frac_pop2)::Matrix{Float64}
    total_frac_means = mean(total_frac,dims = 1)::Matrix{Float64}

    frac2_means = mean(frac_pop2,dims = 1)::Matrix{Float64}

    mnpv_frac = (frac2_means ./ total_frac_means)::Matrix{Float64}

    f_index = findall(x -> x >= index/100, total_frac_means[1,:])::Vector{Int64}

    m_frac_mean = mean(mnpv_frac[:,f_index])::Float64

    mean_last_S = mean(S_pop[:,(gen-50):gen])::Float64
    max_last_S = maximum(S_pop[:,(gen-50):gen])::Float64
    min_last_S = minimum(S_pop[:,(gen-50):gen])::Float64

    mean_last_Z1 = mean(Z1_pop[:,(gen-50):gen])::Float64
    mean_last_Z2 = mean(Z2_pop[:,(gen-50):gen])::Float64

    if mean_last_S <= 12 && mean_last_S >= 1e-2 && max_last_S <= 200 && min_last_S >= 1e-10
        extinct_S = 0::Int
    else
        extinct_S = 1::Int
    end

    if mean_last_Z1 <= 1e-10 && extinct_S == 0
        extinct_SNPV = 1::Int
    else
        extinct_SNPV = 0::Int
    end

    if mean_last_Z2 <= 1e-10 && extinct_S == 0
        extinct_MNPV = 1::Int
    else
        extinct_MNPV = 0::Int
    end

    if extinct_S == 1 #|| extinct_SNPV == 1 || extinct_MNPV == 1
        qual = 0::Int
    else
        qual = 1::Int
    end

    small_df = ll_df::DataFrame

    probs_df = DataFrame()
    for i in 1:length(small_df.id)
        if qual == 1
            if isnan(m_frac_mean) == true
                p1 = NaN
            else
                p1 = dbinom(small_df.MNPV[i],small_df.total[i],m_frac_mean)::Float64
            end
        else
            p1 = 0 #1e-150

        end
 
        temp_df = DataFrame(LL = p1, id = small_df.id[i])

        probs_df = vcat(probs_df,temp_df)
    end

    return probs_df

end

function get_lhoods(ndoug::Int64,sS::Float64,sM::Float64,phi::Float64,rho::Float64,nreps::Int64,sigma::Float64)
    summary_data = DataFrame()
    for rep in 1:nreps

        tree_vals = place_tree_sp(ndoug,number_of_pops)

        output_sim = run_simulation(ndoug,phi,phi,rep,tree_vals,sS,sM,rho,sigma)

        summary_data = vcat(summary_data,output_sim)
    end
    gd = groupby(summary_data, [:id])
    ll_df_new = combine(gd, :LL => mean)

    return ll_df_new

end

function calc_global_ll(ll_arr)
    sum_reps = sum(ll_arr,dims =1)./num_reps
    ll_vals = log.(sum_reps)
    sumll = sum(ll_vals)
    return sumll
end

number_of_params = 4

if number_of_params == 4
    global par_names = ["sS","sM","phi","rho"]
    global par_uppers = [500,500,60,0.9]
    global par_lowers = [0.01,0.01,10,0.05]
    global global_pars = length(par_names)

elseif number_of_params == 3
    global par_names = ["sS","sM","phi","rho"]
    global par_uppers = [500,60,0.9]
    global par_lowers = [0.01,10,0.05]
    global global_pars = length(par_names)
end

global half_njump = 7

global par_uppers_log = log.(par_uppers)
global par_lowers_log = log.(par_lowers)

global iter = 3
global searches = 40

searches*iter*(half_njump + 0.5*half_njump + 0.5)*global_pars*num_sims*2.1/(60*60)

for search in 1:searches

    if rank == 0
        print("Starting line search #"*string(search)*"\n")
    end

    if rank == 0
        global old_lhood = -1e10
        global current_pars =  Array{Float64}(undef, global_pars)

        for j in 1:global_pars
            global current_pars[j] = rand(Uniform(par_lowers[j],par_uppers[j]),1)[1]
        end

        if number_of_params == 4
            track_par = DataFrame(sS = Float64[], sM = Float64[], phi = Float64[], rho = Float64[])
            track_ll = Array{Float64,1}()
    
            track_par_tried = DataFrame(sS = Float64[], sM = Float64[], phi = Float64[], rho = Float64[])
            track_ll_tried = Array{Float64,1}()

        elseif number_of_params == 3
            track_par = DataFrame(s = Float64[], phi = Float64[], rho = Float64[])
            track_ll = Array{Float64,1}()
    
            track_par_tried = DataFrame(s = Float64[], phi = Float64[], rho = Float64[])
            track_ll_tried = Array{Float64,1}()

        end

    end

    for i in 1:iter

        if rank == 0

            jumps = Array{Float64}(undef, 1, global_pars)
            
            for j in 1:global_pars
                adder = rand(Uniform(0,1),1)[1]
                jumper = half_njump + half_njump*adder
                jumps[j] = (par_uppers_log[j] - par_lowers_log[j])/jumper
            end
            
            par_samp = sample(1:global_pars,global_pars,replace = false)
        else
            jumps = Array{Float64}(undef, 1, global_pars)
            par_samp = Array{Float64}(undef,1,global_pars)
            global current_pars = Array{Float64}(undef, global_pars)
        end

        MPI.Barrier(comm)
        par_order = MPI.bcast(par_samp,comm; root = 0)
        MPI.Bcast!(jumps,comm; root = 0)
        
        for k in par_order

            MPI.Barrier(comm)
            par_ptrs = MPI.bcast(current_pars,comm; root = 0)

            par_arr = (par_lowers_log[k]+0.05*(jumps[k])):jumps[k]:(par_uppers_log[k]+0.05*(jumps[k]))

            if rank == 0
                print("Line in log \n")
                print(par_arr)               
                print("\n")
            end

            for j in 1:length(par_arr)

                sites_arr = zeros(1,length(ll_df.Site))
                par_hold = par_ptrs[k]
                global par_ptrs[k] = exp(par_arr[j])

                if rank == 0

                    push!(track_par_tried,par_ptrs)
                end

                if number_of_params == 4
                    likelihood_df = get_lhoods(pdoug,par_ptrs[1],par_ptrs[2],par_ptrs[3],par_ptrs[4],num_sims,0.5)
                elseif number_of_params == 3
                    likelihood_df = get_lhoods(pdoug,par_ptrs[1],par_ptrs[1],par_ptrs[2],par_ptrs[3],num_sims,0.5)
                end


                for l in 1:length(likelihood_df.id)
                    sites_arr[likelihood_df.id[l]] = likelihood_df.LL_mean[l]
                end

                # gather the likelihoods across pdoug and rep sets for pset

                sites_sum = MPI.Gather(sites_arr, comm; root = 0)

                if rank == 0
                    # calculate gobal likelihoods for pset

                    sites_sum = hcat(sites_sum)
                    sites_sum2 = transpose((reshape(sites_sum, 128,size))) 
                    
                    new_lhood = calc_global_ll(sites_sum2) # rand(Uniform(0,1),1)

                    if new_lhood > old_lhood
                        global old_lhood = new_lhood
                        global current_pars[k] = par_ptrs[k]
                    else
                        global current_pars[k] = par_hold
                    end

                    push!(track_ll_tried,new_lhood)

                    push!(track_par,current_pars)
                    push!(track_ll,old_lhood)
                end

            end
        
        end

    end

    if rank == 0
        track_par.ll = track_ll
        track_par.round .= search
        track_par_tried.ll = track_ll_tried
        track_par_tried.round .= search

        dir = "linesearch_op/nts_noext_fix2/"

        if isdir(dir) == false
            mkdir(dir)
        end

        CSV.write(dir*"best_track"*string(search)*".csv", track_par, header = true)
        CSV.write(dir*"tried_track"*string(search)*".csv", track_par_tried, header = true)

        b = now()
        time_elapsed = round((b-a).value/(60000*60),digits = 2)
        print("\n")
        println("Time elapsed during simulation: "*string(time_elapsed)*" hours")
    end


    MPI.Barrier(comm)
end

MPI.Barrier(comm)
