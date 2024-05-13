using Gurobi
using JuMP
using Plots
using LinearAlgebra
using CSV
using DataFrames
using FileIO
using JLD2
using StatsPlots
using Statistics
using Random

include("functions.jl")

function main()
    total_time_days = 30; # should be a multiple of 30
    timestep_hours = 0.25; # dataset timestep is also 0.25 hours
    number_of_timesteps = floor(Int,total_time_days*24/timestep_hours);
    # Type of data: 0 = perfect, 1 = 24h average, 1.1 = 1 month average, 2 = shuffled days, 2.1 = shuffled days 24h average
    simulation_data_number = 0;
    model_data_number = 1;
    # Type of model
    model_number = 1; # 0 - detailed forecast MILP (DFM), 1 - average forecast model solved using solver, 2 - optimal benchmark MILP (OBM), 4 - average forecast greeedy (AFG) model
    threshold_constant_days = 1;
    load_priority_order = [2;4;3;1];
    start_day = 1;
    # Add cases
    case_array = [2] # 0 = no energy management, 1 = fixed thresholds, 2 = optimal thresholds
    # Add real recharge monthly frequency
    real_recharge_logic_array = [1];
    # Add real recharge amount in percent
    recharge_percent_array = [70,80,90,100];    
    experiment_parameters = Dict("total_time_days"=> total_time_days, "timestep_hours"=>timestep_hours, "simulation_data_number"=>simulation_data_number,"model_data_number"=>model_data_number,"model_number"=>model_number,"threshold_constant_days"=>threshold_constant_days,"load_priority_order"=>load_priority_order,"start_day"=>start_day,"case_array"=>case_array,"real_recharge_logic_array"=>real_recharge_logic_array,"recharge_percent_array"=>recharge_percent_array)    
    if model_number == 0 || model_number == 1 || model_number == 3 || model_number == 4
        time_array, simulation_results, model_results, parameters = run_threshold_model_simulation(experiment_parameters);
    elseif model_number == 2
        run_hem_model_simulation(experiment_parameters);        
    end
end

# function to run home energy manager based models and simulations
function run_hem_model_simulation(experiment_parameters)
    total_time_days = experiment_parameters["total_time_days"]
    timestep_hours = experiment_parameters["timestep_hours"]
    model_number = experiment_parameters["model_number"]
    simulation_data_number = experiment_parameters["simulation_data_number"]
    model_data_number = experiment_parameters["model_data_number"]
    start_day = experiment_parameters["start_day"]
    case_array = experiment_parameters["case_array"]
    real_recharge_logic_array = experiment_parameters["real_recharge_logic_array"];
    recharge_percent_array = experiment_parameters["recharge_percent_array"];
    load_priority_order = experiment_parameters["load_priority_order"]   
    number_of_timesteps = floor(Int,total_time_days*24/timestep_hours) 
    # To store metrics
    number_of_results = length(real_recharge_logic_array)*length(recharge_percent_array)*length(case_array)
    metrics_results = zeros(number_of_results,12)
    metrics_results_index = 0
    model_metrics_results = zeros(number_of_results,11)
    model_metrics_results_index = 0    
    disconnection_counter = 0
    directory = pwd()

    for case in case_array
    for recharge_percent in recharge_percent_array
    for real_recharge_logic in real_recharge_logic_array
    
    if model_number == 2
    # Initialize load parameters
    demand_state, demand_power_W, monthly_demand_energy_kWh, number_of_loads, load_priority_factor, daily_energy_Wh, daily_energy_array_Wh = initialize_load_parameters(total_time_days, timestep_hours, load_priority_order, start_day, directory);             
    # Real recharge parameters
    daily_real_recharge, daily_real_cost, cost_perkWh, real_recharge_schedule, real_cost_schedule, monthly_recharge_amount = initialize_real_recharge_parameters(total_time_days, timestep_hours, monthly_demand_energy_kWh, recharge_percent, real_recharge_logic);
    # Generate modified load data
    demand_power_W_modified, demand_state_modified, daily_energy_array_Wh_modified, daily_energy_Wh_modified = generate_modified_data(demand_power_W,demand_state,number_of_timesteps,-1,timestep_hours,total_time_days,model_data_number,number_of_loads,directory)    
    # Modify load data for model if necessary
    if model_data_number == 0
        model_demand_state = demand_state
        model_demand_power_W = demand_power_W
    elseif model_data_number == 1 || model_data_number == 2 || model_data_number == 2.1 || model_data_number == 1.1
        model_demand_state = demand_state_modified
        model_demand_power_W = demand_power_W_modified
    end
    # Calculate actuation states from model    
    model_actuation_state, model_psf, model_esf, model_load_tsf, solve_time_model = compute_hem_schedule(total_time_days, timestep_hours, model_demand_state, model_demand_power_W, cost_perkWh, number_of_loads, load_priority_factor, daily_real_recharge, directory)
    # Modify load data for simulation if necessary
    if simulation_data_number == 1
        demand_power_W = demand_power_W_modified
        demand_state = demand_state_modified 
    end    
    # Starting timesteps
    day_start_timesteps = zeros(Int,total_time_days,1);
    for i in 1:total_time_days
        day_start_timesteps[i] = floor(Int,(24/timestep_hours)*(i-1) + 1);
    end
    real_wallet = zeros(number_of_timesteps)
    actuation_state = zeros(number_of_loads,number_of_timesteps)
    # Start simulation
    for t in 1:number_of_timesteps
        if t in day_start_timesteps
            real_wallet[t] = real_wallet[t] + daily_real_recharge[floor(Int,1 + t/(24/timestep_hours))]
        end
        if real_wallet[t] < 0
            actuation_state[:,t] = zeros(number_of_loads,1)
        else
            for k in 1:number_of_loads
                if demand_state[k,t] == 0
                    actuation_state[k,t] = 0
                else
                    actuation_state[k,t] = model_actuation_state[k,t]
                end
            end
        end
        if t < number_of_timesteps
            real_wallet[t+1] = real_wallet[t] - sum(actuation_state[k,t]*demand_power_W[k,t] for k in 1:number_of_loads)*timestep_hours*(cost_perkWh/1000)
        end
    end

    disconnection_counter = 0
    for t = 1:number_of_timesteps
        # Count disconnections
        if t == 1
            if real_wallet[t] < 0 
                disconnection_counter = disconnection_counter + 1
            end
        else
            if real_wallet[t-1] > 0 && real_wallet[t] <0
                disconnection_counter = disconnection_counter + 1
            end
        end    
    end  

    priority_sf = 100*sum(load_priority_factor[k]*(sum(value.(actuation_state[k,t]) for t in 1:number_of_timesteps)/(sum(demand_state[k,t] for t in 1:number_of_timesteps))) for k in 1:number_of_loads)
    energy_sf = sum(sum(value.(actuation_state[k,t])*demand_power_W[k,t] for k in 1:number_of_loads) for t in 1:number_of_timesteps)*100/sum(sum(demand_state[k,t]*demand_power_W[k,t] for k in 1:number_of_loads) for t in 1:number_of_timesteps)
    # SF per load
    time_service_factor_perload = zeros(number_of_loads,1);
    for k in 1:number_of_loads
        time_service_factor_perload[k] = (sum(value.(actuation_state[k,t]) for t in 1:number_of_timesteps))*100/(sum(demand_state[k,t] for t in 1:number_of_timesteps))     
    end    
    
    # Add to metrics matrix
    metrics_results_index = metrics_results_index + 1
    metrics_results[metrics_results_index,1] = case
    metrics_results[metrics_results_index,2] = real_recharge_logic
    metrics_results[metrics_results_index,3] = recharge_percent
    metrics_results[metrics_results_index,4] = priority_sf
    metrics_results[metrics_results_index,5:8] = time_service_factor_perload
    metrics_results[metrics_results_index,9] = -1
    metrics_results[metrics_results_index,10] = energy_sf
    metrics_results[metrics_results_index,11] = disconnection_counter   
    metrics_results[metrics_results_index,12] = solve_time_model       

    # Add to model metrics matrix
    model_metrics_results_index = model_metrics_results_index + 1
    model_metrics_results[model_metrics_results_index,1] = case
    model_metrics_results[model_metrics_results_index,2] = real_recharge_logic
    model_metrics_results[model_metrics_results_index,3] = recharge_percent
    model_metrics_results[model_metrics_results_index,4] = model_psf
    model_metrics_results[model_metrics_results_index,5:8] = model_load_tsf
    model_metrics_results[model_metrics_results_index,9] = -1
    model_metrics_results[model_metrics_results_index,10] = model_esf
    model_metrics_results[model_metrics_results_index,11] = -1       

    end

    end
    end
    end

    
    column_names = ["case", "real_recharge_logic", "recharge_percent", "priority_sf", "load1_tsf", "load2_tsf", "load3_tsf", "load4_tsf", "time_sf", "energy_sf", "disconnections", "solve time"]
    column_names_model = ["case", "real_recharge_logic", "recharge_percent", "priority_sf", "load1_tsf", "load2_tsf", "load3_tsf", "load4_tsf", "time_sf", "energy_sf", "disconnections"]
    metrics_results_df = DataFrame(metrics_results, column_names)
    model_metrics_results_df = DataFrame(model_metrics_results, column_names_model)
    # Write metrics to csv file
    # CSV.write(directory * "metrics.csv", metrics_results_df);
    # CSV.write(directory * "model_metrics.csv", model_metrics_results_df);

end

# function to run threshold-based models and simulations
function run_threshold_model_simulation(experiment_parameters)
    ###### Read experiment parameters ######
    simulation_time_days = experiment_parameters["total_time_days"]
    simulation_timestep_hours = experiment_parameters["timestep_hours"]
    number_of_timesteps = floor(Int,simulation_time_days*24/simulation_timestep_hours);
    time_array = (0:(number_of_timesteps-1))*simulation_timestep_hours;    
    simulation_data_number = experiment_parameters["simulation_data_number"];
    model_data_number = experiment_parameters["model_data_number"];
    model_number = experiment_parameters["model_number"];
    threshold_constant_days = experiment_parameters["threshold_constant_days"];
    load_priority_order = experiment_parameters["load_priority_order"];
    start_day = experiment_parameters["start_day"];
    case_array = experiment_parameters["case_array"]
    real_recharge_logic_array = experiment_parameters["real_recharge_logic_array"];
    recharge_percent_array = experiment_parameters["recharge_percent_array"];
    # To store metrics
    number_of_results = length(real_recharge_logic_array)*length(recharge_percent_array)*length(case_array)
    metrics_results = zeros(number_of_results,12)
    metrics_results_index = 0
    model_metrics_results = zeros(number_of_results,11)
    model_metrics_results_index = 0
    # Data storage directory
    directory = pwd()
    # Creating variables to return results
    simulation_results = 0
    model_results = 0
    parameters = 0
    
    ###### Set up simulation for each combination of case, recharge_percent, real_recharge_logic ######
    for case in case_array
    for recharge_percent in recharge_percent_array
    for real_recharge_logic in real_recharge_logic_array
    
    # Initialize load parameters
    demand_state, demand_power_W, monthly_demand_energy_kWh, number_of_loads, load_priority_factor, daily_energy_Wh, daily_energy_array_Wh = initialize_load_parameters(simulation_time_days, simulation_timestep_hours, load_priority_order, start_day, directory);
    # Recharge data
    daily_real_recharge, daily_real_cost, cost_perkWh, real_recharge_schedule, real_cost_schedule, monthly_recharge_amount = initialize_real_recharge_parameters(simulation_time_days, simulation_timestep_hours, monthly_demand_energy_kWh, recharge_percent, real_recharge_logic);
    daily_virtual_recharge = compute_daily_virtual_update(simulation_time_days, real_recharge_schedule)
    daily_virtual_cost = compute_daily_virtual_update(simulation_time_days, real_cost_schedule)
    # Time data
    simulation_time_months = floor(Int, simulation_time_days/30)
    month_start_timesteps = zeros(Int, simulation_time_months,1);
    for i in 1:simulation_time_months
        month_start_timesteps[i] = floor(Int,30*(24/simulation_timestep_hours)*(i-1) + 1);
    end 
    # Horizon
    horizon_time_days = 30
    horizon_start_days = zeros(Int,ceil(Int,simulation_time_days/horizon_time_days),1);
    horizon_start_timesteps = zeros(Int, ceil(Int,simulation_time_days/horizon_time_days),1);
    # Update
    update_time_days = 30
    if model_number == 0
        update_time_days = 30
    elseif model_number == 1
        update_time_days = 30
    elseif model_number == 3
        update_time_days = 30
    elseif model_number == 4
        update_time_days = 30
    end  
    update_start_days = zeros(Int,ceil(Int,simulation_time_days/update_time_days),1);
    update_start_timesteps = zeros(Int, ceil(Int,simulation_time_days/update_time_days),1);
    # Starting timesteps
    day_start_timesteps = zeros(Int,simulation_time_days,1);
    for i in 1:simulation_time_days
        day_start_timesteps[i] = floor(Int,(24/simulation_timestep_hours)*(i-1) + 1);
        if mod(i,horizon_time_days) == 1
            horizon_start_days[floor(Int,div(i,horizon_time_days))+1] = i;
            horizon_start_timesteps[floor(Int,div(i,horizon_time_days))+1] = floor(Int,(24/simulation_timestep_hours)*(i-1) + 1);
        end
        if update_time_days == 1
            update_start_days[i] = i;
            update_start_timesteps[i] = floor(Int,(24/simulation_timestep_hours)*(i-1) + 1);
        elseif mod(i,update_time_days) == 1
            update_start_days[floor(Int,div(i,update_time_days))+1] = i;
            update_start_timesteps[floor(Int,div(i,update_time_days))+1] = floor(Int,(24/simulation_timestep_hours)*(i-1) + 1);
        end
    end

    # Creating simulation variables to store results
    actuation_state = zeros(number_of_loads,number_of_timesteps);
    real_wallet = zeros(number_of_timesteps);
    real_enable = zeros(number_of_loads,number_of_timesteps);
    virtual_wallet = zeros(number_of_timesteps);
    virtual_enable = zeros(number_of_loads,number_of_timesteps);
    virtual_wallet_thresholds = zeros(number_of_loads,number_of_timesteps);
    computed_daily_virtual_recharge = zeros(update_time_days);
    # Auxiliary variable
    virtual_recharge = 0
    # For model_number = 0
    model_real_wallet = zeros(number_of_timesteps);
    model_virtual_wallet = zeros(number_of_timesteps);
    model_virtual_enable = zeros(number_of_loads,number_of_timesteps);
    model_real_enable = zeros(number_of_loads,number_of_timesteps);
    model_actuation_state = zeros(number_of_loads,number_of_timesteps);
    model_demand_state = zeros(number_of_loads,number_of_timesteps);
    model_virtual_wallet_thresholds = zeros(number_of_loads, simulation_time_days)
    # For model_number = 1
    daily_load_sf = zeros(number_of_loads,simulation_time_days)
    daily_psf = zeros(simulation_time_days,1)
    model_daily_load_sf = zeros(number_of_loads,simulation_time_days)
    model_daily_psf = zeros(simulation_time_days,1)
    # Counting disconnections
    disconnection_counter = 0;
    # Solution time
    solve_time_model = 0;
    # Wallet initial conditions
    real_wallet[1] = 0;
    virtual_wallet[1] = 0; 
    # Generate modified load data
    demand_power_W_modified, demand_state_modified, daily_energy_array_Wh_modified, daily_energy_Wh_modified = generate_modified_data(demand_power_W,demand_state,number_of_timesteps,horizon_time_days,simulation_timestep_hours,simulation_time_days,model_data_number,number_of_loads,directory)    
    # Modify load data for simulation if necessary
    if simulation_data_number == 1
        demand_power_W = demand_power_W_modified
        demand_state = demand_state_modified 
    end

    # Start simulation
    for t = 1:number_of_timesteps-1
        # case = 2 -> optimized thresholds
        if case == 2
            # Solve optimization problem
            if t in update_start_timesteps
                time_data = Dict("number_of_days" => horizon_time_days, "timestep_hours" => simulation_timestep_hours, "update_time_days" => update_time_days);
                Nw = horizon_time_days*floor(Int,24/simulation_timestep_hours) # horizon duration in timesteps, reduced when current day + horizon_time_days exceeds simulation time
                Nw_days = horizon_time_days # horizon duration in days, reduced when current day + horizon_time_days exceeds simulation time
                if t + Nw - 1 > number_of_timesteps
                    Nw = number_of_timesteps + 1 - t;
                    Nw_days = floor(Int,Nw/(24/simulation_timestep_hours));
                    time_data["number_of_days"] = floor(Int,Nw/(24/simulation_timestep_hours));
                end
                if model_data_number == 1 || model_data_number == 2 || model_data_number == 2.1 || model_data_number == 1.1
                    demand_data = Dict("number_of_loads" => number_of_loads, "demand_power_W" => demand_power_W_modified[:,t:t+Nw-1], "demand_state" => demand_state_modified[:,t:t+Nw-1], "load_priority_factor" => load_priority_factor);
                else
                    demand_data = Dict("number_of_loads" => number_of_loads, "demand_power_W" => demand_power_W[:,t:t+Nw-1], "demand_state" => demand_state[:,t:t+Nw-1], "load_priority_factor" => load_priority_factor);
                end
                real_wallet_balance = real_wallet[t];
                virtual_wallet_balance = virtual_wallet[t]
                day_index = floor(Int,(t-1)/(24/simulation_timestep_hours))+1;
                real_wallet_data = Dict("real_wallet_balance" => real_wallet_balance, "daily_real_recharge" => daily_real_recharge[day_index:day_index+Nw_days-1], "daily_real_cost" => daily_real_cost[day_index:day_index+Nw_days-1], "cost_perkWh" => cost_perkWh, "recharge_percent" => recharge_percent); 
                virtual_wallet_data = Dict("virtual_wallet_balance" => virtual_wallet_balance, "daily_virtual_recharge" => daily_virtual_recharge[day_index:day_index+Nw_days-1], "daily_virtual_cost" => daily_virtual_cost[day_index:day_index+Nw_days-1], "threshold_constant_days" => threshold_constant_days);
                Nu = floor(Int,update_time_days*24/simulation_timestep_hours)-1 # update time in timesteps
                # Invoke required model
                if model_number == 0
                    virtual_wallet_thresholds[:,t:t+Nu], model_real_wallet[t:t+Nu], model_virtual_wallet[t:t+Nu], model_virtual_enable[:,t:t+Nu], model_real_enable[:,t:t+Nu], model_actuation_state[:,t:t+Nu], model_demand_state[:,t:t+Nu], model_virtual_wallet_thresholds[:,floor(Int,t/(24/simulation_timestep_hours))+1:floor(Int,t/(24/simulation_timestep_hours)) + update_time_days], solve_time_model = compute_thresholds(time_data, demand_data, real_wallet_data, virtual_wallet_data);
                elseif model_number == 1
                    virtual_wallet_thresholds[:,t:t+Nu], computed_daily_virtual_recharge, model_daily_load_sf[:,floor(Int,t/(24/simulation_timestep_hours))+1:floor(Int,t/(24/simulation_timestep_hours)) + update_time_days], model_real_wallet[t:t+Nu], model_virtual_wallet[t:t+Nu], model_virtual_enable[:,t:t+Nu], model_real_enable[:,t:t+Nu], model_actuation_state[:,t:t+Nu], model_demand_state[:,t:t+Nu], model_virtual_wallet_thresholds[:,floor(Int,t/(24/simulation_timestep_hours))+1:floor(Int,t/(24/simulation_timestep_hours)) + update_time_days], solve_time_model = compute_virtualrecharge_thresholds(time_data, demand_data, real_wallet_data, virtual_wallet_data, directory);                    
                elseif model_number == 3
                    virtual_wallet_thresholds[:,t:t+Nu], model_real_wallet[t:t+Nu], model_virtual_wallet[t:t+Nu], model_virtual_enable[:,t:t+Nu], model_real_enable[:,t:t+Nu], model_actuation_state[:,t:t+Nu], model_demand_state[:,t:t+Nu], model_virtual_wallet_thresholds[:,floor(Int,t/(24/simulation_timestep_hours))+1:floor(Int,t/(24/simulation_timestep_hours)) + update_time_days], computed_daily_virtual_recharge = compute_thresholds_extra_constraints(time_data, demand_data, real_wallet_data, virtual_wallet_data);
                elseif model_number == 4
                    virtual_wallet_thresholds[:,t:t+Nu], computed_daily_virtual_recharge, model_daily_load_sf[:,floor(Int,t/(24/simulation_timestep_hours))+1:floor(Int,t/(24/simulation_timestep_hours)) + update_time_days], model_real_wallet[t:t+Nu], model_virtual_wallet[t:t+Nu], model_virtual_enable[:,t:t+Nu], model_real_enable[:,t:t+Nu], model_actuation_state[:,t:t+Nu], model_demand_state[:,t:t+Nu], model_virtual_wallet_thresholds[:,floor(Int,t/(24/simulation_timestep_hours))+1:floor(Int,t/(24/simulation_timestep_hours)) + update_time_days], solve_time_model = compute_enable_durations(time_data, demand_data, real_wallet_data, virtual_wallet_data, directory)
                end
            end
        # case = 1 -> fixed thresholds
        elseif case == 1
            if t in month_start_timesteps
                descending_priority_indices = sortperm(load_priority_factor[:,1], rev=true)
                i = 0
                for k in descending_priority_indices
                    i = i+1
                    virtual_wallet_thresholds[k,t:t+floor(Int,(30*24/simulation_timestep_hours))-1] = 0.05*monthly_recharge_amount*(i/number_of_loads)*ones(1,floor(Int,(30*24/simulation_timestep_hours)))
                end
            end
        end
        # Update wallet
        # Daily recharge and daily costs
        if t in day_start_timesteps
            if ((model_number == 1) || (model_number == 3) || (model_number == 4))  && (case == 2)
                day_index = floor(Int,1 + (t)*simulation_timestep_hours/24)
                update_window_index = mod(day_index,update_time_days)
                if update_window_index == 0
                    update_window_index = update_time_days
                end
                virtual_recharge = computed_daily_virtual_recharge[update_window_index];
            else
                virtual_recharge = daily_virtual_recharge[floor(Int,1 + (t)*simulation_timestep_hours/24)];
            end    
            real_wallet[t] = real_wallet[t] + daily_real_recharge[floor(Int,1 + (t)*simulation_timestep_hours/24)] - daily_real_cost[floor(Int,1 + (t)*simulation_timestep_hours/24)];
            virtual_wallet[t] = virtual_wallet[t] + virtual_recharge - daily_virtual_cost[floor(Int,1 + (t)*simulation_timestep_hours/24)];
        end        
        actuation_state[:,t], virtual_enable[:,t], real_enable[:,t] = compute_actuation_state(number_of_loads,virtual_wallet[t], virtual_wallet_thresholds[:,t], real_wallet[t], demand_state[:,t], case);   
        # Determine wallet values for next time step
        virtual_wallet[t+1] = virtual_wallet[t] - sum(actuation_state[k,t]*demand_power_W[k,t]*simulation_timestep_hours*(cost_perkWh/1000) for k in 1:number_of_loads);
        real_wallet[t+1] = real_wallet[t] - sum(actuation_state[k,t]*demand_power_W[k,t]*simulation_timestep_hours*(cost_perkWh/1000) for k in 1:number_of_loads);
        # Prevent disconnections
        if (case == 1) || (case == 2)
            if real_wallet[t+1] < 0
                actuation_state[:,t] = zeros(number_of_loads,1)
                virtual_wallet[t+1] = virtual_wallet[t] - sum(actuation_state[k,t]*demand_power_W[k,t]*simulation_timestep_hours*(cost_perkWh/1000) for k in 1:number_of_loads);
                real_wallet[t+1] = real_wallet[t] - sum(actuation_state[k,t]*demand_power_W[k,t]*simulation_timestep_hours*(cost_perkWh/1000) for k in 1:number_of_loads); 
            end
        end
    end
    # Compute actuation_state for last timestep since the for loop does not
    actuation_state[:,number_of_timesteps], virtual_enable[:,number_of_timesteps], real_enable[:,number_of_timesteps] = compute_actuation_state(number_of_loads,virtual_wallet[number_of_timesteps], virtual_wallet_thresholds[:,number_of_timesteps], real_wallet[number_of_timesteps], demand_state[:,number_of_timesteps],case);
    for t = 1:number_of_timesteps
        # Count disconnections
        if t == 1
            if real_wallet[t] < 0 
                disconnection_counter = disconnection_counter + 1
            end
        else
            if real_wallet[t-1] > 0 && real_wallet[t] <0
                disconnection_counter = disconnection_counter + 1
                print("t = ")
                println(t)                
            end
        end    
    end

    for d = 1:simulation_time_days
        for k = 1:number_of_loads
            if sum(demand_state[k,t] for t in floor(Int,(24/simulation_timestep_hours)*(d-1)+1):floor(Int,(24/simulation_timestep_hours)*(d))) == 0
                daily_load_sf[k,d] = -1
            else
                daily_load_sf[k,d] = 100 * sum(actuation_state[k,t] for t in floor(Int,(24/simulation_timestep_hours)*(d-1)+1):floor(Int,(24/simulation_timestep_hours)*(d))) / sum(demand_state[k,t] for t in floor(Int,(24/simulation_timestep_hours)*(d-1)+1):floor(Int,(24/simulation_timestep_hours)*(d)))
            end
        end
        for k = 1:number_of_loads
            if daily_load_sf[k,d] == -1
                daily_psf[d,1] = daily_psf[d] + load_priority_factor[k]*100
            else
                daily_psf[d,1] = daily_psf[d] + load_priority_factor[k]*daily_load_sf[k,d]
            end
            if model_daily_load_sf[k,d] == -1
                model_daily_psf[d,1] = model_daily_psf[d] + load_priority_factor[k]*100
            else
                model_daily_psf[d,1] = model_daily_psf[d] + load_priority_factor[k]*model_daily_load_sf[k,d]
            end
        end
    end
    column_names = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30"]    
    daily_load_sf_results_df = DataFrame([daily_load_sf; daily_psf'; model_daily_load_sf; model_daily_psf'],column_names)

    # Store results
    simulation_results = Dict("actuation_state" => actuation_state, "demand_state" => demand_state, "virtual_enable" => virtual_enable, "real_wallet" => real_wallet, "virtual_wallet" => virtual_wallet, "virtual_wallet_thresholds" => virtual_wallet_thresholds);
    model_results = Dict("actuation_state" => model_actuation_state, "demand_state" => model_demand_state, "virtual_enable" => model_virtual_enable, "real_wallet" => model_real_wallet, "virtual_wallet" => model_virtual_wallet, "virtual_wallet_thresholds" => model_virtual_wallet_thresholds );
    parameters = Dict("number_of_loads" => number_of_loads, "number_of_timesteps" => number_of_timesteps, "demand_power_W" => demand_power_W, "load_priority_factor" => load_priority_factor, "case" => case, "recharge_percent" => recharge_percent, "real_recharge_logic" => real_recharge_logic, "simulation_timestep_hours" => simulation_timestep_hours, "model_data_number" => model_data_number, "simulation_data_number" => simulation_data_number)

    # Add to metrics matrix
    metrics_results_index = metrics_results_index + 1
    metrics_results[metrics_results_index,1] = case
    metrics_results[metrics_results_index,2] = real_recharge_logic
    metrics_results[metrics_results_index,3] = recharge_percent
    # Compute metrics and add to metrics matrix
    priority_sf, load_tsf, time_sf, energy_sf, model_priority_sf, model_load_tsf, model_time_sf, model_energy_sf = process_results(time_array, simulation_results, model_results, parameters, directory);
    metrics_results[metrics_results_index,4] = priority_sf
    metrics_results[metrics_results_index,5:8] = load_tsf
    metrics_results[metrics_results_index,9] = time_sf
    metrics_results[metrics_results_index,10] = energy_sf
    metrics_results[metrics_results_index,11] = disconnection_counter
    metrics_results[metrics_results_index,12] = solve_time_model
    # Add to model metrics matrix
    model_metrics_results_index = model_metrics_results_index + 1
    model_metrics_results[model_metrics_results_index,1] = case
    model_metrics_results[model_metrics_results_index,2] = real_recharge_logic
    model_metrics_results[model_metrics_results_index,3] = recharge_percent
    model_metrics_results[model_metrics_results_index,4] = model_priority_sf
    model_metrics_results[model_metrics_results_index,5:8] = model_load_tsf
    model_metrics_results[model_metrics_results_index,9] = model_time_sf
    model_metrics_results[model_metrics_results_index,10] = model_energy_sf
    model_metrics_results[model_metrics_results_index,11] = -1
    

    end
    end
    end

    column_names = ["case", "real_recharge_logic", "recharge_percent", "priority_sf", "load1_tsf", "load2_tsf", "load3_tsf", "load4_tsf", "time_sf", "energy_sf", "disconnections", "solve time"]
    column_names_model = ["case", "real_recharge_logic", "recharge_percent", "priority_sf", "load1_tsf", "load2_tsf", "load3_tsf", "load4_tsf", "time_sf", "energy_sf", "disconnections"]
    metrics_results_df = DataFrame(metrics_results, column_names)
    model_metrics_results_df = DataFrame(model_metrics_results, column_names_model)
    # Write metrics to csv file
    # CSV.write(directory * "metrics.csv", metrics_results_df);
    # CSV.write(directory * "model_metrics.csv", model_metrics_results_df);

    return time_array, simulation_results, model_results, parameters;
end

main();