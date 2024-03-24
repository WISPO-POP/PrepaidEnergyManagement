function initialize_load_parameters(simulation_time_days, timestep_hours, load_priority_order, start_day, directory)
    number_of_loads = 4;
    number_of_timesteps = floor(Int,simulation_time_days*24/timestep_hours);
    # Read load data file
    load_data = CSV.read(pwd()*"\\..\\datasets\\monthly-experiment-data.csv",DataFrame);

    # Extract required data 
    lb = (start_day - 1)*floor(Int,24/timestep_hours) + 1;; # lower bound
    ub = lb + number_of_timesteps - 1; # upper bound
    time_array = load_data.local_15min[lb:ub];
    demand_power_W = zeros(number_of_loads,number_of_timesteps)
    demand_power_W[1,:] = 1000*load_data.air1[lb:ub]; # converting power values from kW to W
    demand_power_W[2,:] = 1000*load_data.clotheswasher_dryg1[lb:ub];
    demand_power_W[3,:] = 1000*load_data.microwave1[lb:ub];
    demand_power_W[4,:] = 1000*load_data.refrigerator1[lb:ub];
    demand_state = zeros(number_of_loads,number_of_timesteps); # demand signal (1 for load ON, 0 for load OFF)
    for t = 1:number_of_timesteps
        # De-noising step; for air compressor values below 1% of maximum are assumed to be noise
        # and for washer, microwave, refrigerator values below 10% of maximum are assumed to be noise
        # These % noise thresholds are assumptions and were decided by looking at data plots
        # Then multiply power demand by demand state
        demand_state[1,t] = ceil((demand_power_W[1,t] - 0.01*maximum(demand_power_W[1,:]))/(abs(demand_power_W[1,t] - 0.01*maximum(demand_power_W[1,:])) + 0.001));
        demand_power_W[1,t] = demand_power_W[1,t]*demand_state[1,t];
        demand_state[2,t] = ceil((demand_power_W[2,t] - 0.1*maximum(demand_power_W[2,:]))/(abs(demand_power_W[2,t] - 0.1*maximum(demand_power_W[2,:])) + 0.001));
        demand_power_W[2,t] = demand_power_W[2,t]*demand_state[2,t];
        demand_state[3,t] = ceil((demand_power_W[3,t] - 0.1*maximum(demand_power_W[3,:]))/(abs(demand_power_W[3,t] - 0.1*maximum(demand_power_W[3,:])) + 0.001));
        demand_power_W[3,t] = demand_power_W[3,t]*demand_state[3,t];
        demand_state[4,t] = ceil((demand_power_W[4,t] - 0.1*maximum(demand_power_W[4,:]))/(abs(demand_power_W[4,t] - 0.1*maximum(demand_power_W[4,:])) + 0.001));
        demand_power_W[4,t] = demand_power_W[4,t]*demand_state[4,t];
    end

    # Plots for demand state and power demand for all loads
    time_array = (0:number_of_timesteps-1)*timestep_hours/24;
    air1_power_plot = plot(time_array, demand_power_W[1,:], title = "Load A: air compressor", titlefontsize = 10, yticks = [0, 1500], ytickfontsize = 10, xtickfontsize = 10, legend = false)
    air1_demandstate_plot = plot(time_array, demand_state[1,:], title = "Load A demand state", titlefontsize = 10, yticks = [0, 1500], ytickfontsize = 10, xtickfontsize = 10, legend = false)    
    clotheswasher_dryg1_power_plot = plot(time_array, demand_power_W[2,:], title = "Load B: washing machine", yticks = [0, 400], ytickfontsize = 10, xtickfontsize = 10, titlefontsize = 10, legend = false)
    clotheswasher_demandstate_plot = plot(time_array, demand_state[2,:], title = "Load B: demand state", yticks = [0, 400], ytickfontsize = 10, xtickfontsize = 10, titlefontsize = 10, legend = false)
    microwave1_power_plot = plot(time_array,demand_power_W[3,:], title = "Load C: microwave", yticks = [0, 1200], ytickfontsize = 10, xtickfontsize = 10, titlefontsize = 10, legend = false)
    microwave_demandstate_plot = plot(time_array,demand_state[3,:], title = "Load C: demand state", yticks = [0, 1200], ytickfontsize = 10, xtickfontsize = 10, titlefontsize = 10, legend = false)
    refrigerator1_power_plot = plot(time_array,demand_power_W[4,:], title = "Load D: refrigerator", ytickfontsize = 10, xtickfontsize = 10, xlabel = "days", xguidefontsize = 12, titlefontsize = 10, legend = false)
    refrigerator_demandstate_plot = plot(time_array,demand_state[4,:], title = "Load D: demand state", ytickfontsize = 10, xtickfontsize = 10, xlabel = "days", xguidefontsize = 12, titlefontsize = 10, legend = false)
    load_power_plot1 = plot(air1_power_plot,air1_demandstate_plot,layout = (2,1))
    load_power_plot2 = plot(clotheswasher_dryg1_power_plot,clotheswasher_demandstate_plot,layout = (2,1))
    load_power_plot3 = plot(microwave1_power_plot,microwave_demandstate_plot,layout = (2,1))
    load_power_plot4 = plot(refrigerator1_power_plot,refrigerator_demandstate_plot,layout = (2,1))

    # Compute load priority factors
    factor = 0;
    load_priority_factor = zeros(number_of_loads,1);
    for k = 1:number_of_loads
        factor = factor + 1/load_priority_order[k];
    end
    factor = 1/factor;
    for k = 1:number_of_loads
        load_priority_factor[k] = factor/load_priority_order[k];
    end

    # Compute monthly energy use
    month_start_days = zeros(Int, simulation_time_days*30,1);
    simulation_time_months = floor(Int,simulation_time_days/30);
    month_start_timesteps = zeros(Int, simulation_time_months,1);
    for i in 1:simulation_time_months
        month_start_days[i] = floor(Int,30*(i-1) + 1);
        month_start_timesteps[i] = floor(Int,30*(24/timestep_hours)*(i-1) + 1);
    end  
    
    demand_energy_Wh = zeros(number_of_loads,simulation_time_months)
    for i in 1:simulation_time_months
        for k in 1:number_of_loads
            demand_energy_Wh[k,i] = sum(demand_power_W[k,t]*demand_state[k,t] for t in month_start_timesteps[i]:floor(Int,month_start_timesteps[i]+30*(24/timestep_hours)-1))*timestep_hours;
        end
    end     
    monthly_energy_Wh = sum(demand_energy_Wh[k,:] for k in 1:number_of_loads)
    monthly_energy_kWh = monthly_energy_Wh/1000;
    total_energy_kWh = sum(monthly_energy_kWh);    
    # Compute daily energy use
    day_start_timesteps = zeros(Int, simulation_time_days,1);
    for i in 1:simulation_time_days
        day_start_timesteps[i] = floor(Int,(24/timestep_hours)*(i-1) + 1);
    end
    demand_daily_energy_Wh = zeros(number_of_loads,simulation_time_days)
    for i in 1:simulation_time_days
        for k in 1:number_of_loads
            demand_daily_energy_Wh[k,i] = sum(demand_power_W[k,t]*demand_state[k,t] for t in day_start_timesteps[i]:floor(Int,day_start_timesteps[i]+(24/timestep_hours)-1))*timestep_hours;
        end
    end
    daily_energy_array_Wh = [demand_daily_energy_Wh[1,:] demand_daily_energy_Wh[2,:] demand_daily_energy_Wh[3,:] demand_daily_energy_Wh[4,:]]
    daily_energy_Wh = zeros(simulation_time_days,1)
    for i in 1:simulation_time_days
        daily_energy_Wh[i,1] = sum(demand_daily_energy_Wh[k,i] for k in 1:number_of_loads)
    end
    #groupedbar([air1_daily_energy_Wh clotheswasher_dryg1_daily_energy_Wh microwave1_daily_energy_Wh refrigerator1_daily_energy_Wh], bar_position = :dodge, bar_width=0.7, legend=:bottomright, label = ["compressor" "washer" "microwave" "refrigerator"], xlabel = "Days", ylabel = "energy (Wh)")
    air1_daily_energy_Wh_plot = plot(1:simulation_time_days,demand_daily_energy_Wh[1,:])
    clotheswasher_dryg1_daily_energy_Wh_plot = plot(1:simulation_time_days,demand_daily_energy_Wh[2,:])
    microwave1_daily_energy_Wh_plot = plot(1:simulation_time_days,demand_daily_energy_Wh[3,:])
    refrigerator1_daily_energy_Wh_plot = plot(1:simulation_time_days,demand_daily_energy_Wh[4,:])
    daily_energy_Wh_plot = plot(air1_daily_energy_Wh_plot,clotheswasher_dryg1_daily_energy_Wh_plot,microwave1_daily_energy_Wh_plot,refrigerator1_daily_energy_Wh_plot,layout=(4,1))

    return demand_state, demand_power_W, monthly_energy_kWh, number_of_loads, load_priority_factor, daily_energy_Wh, daily_energy_array_Wh
end

function initialize_real_recharge_parameters(simulation_time_days, timestep_hours, monthly_energy_kWh, recharge_percent, real_recharge_logic)
    simulation_time_months = floor(Int, simulation_time_days/30);
    # Electricity rates
    # according to my MG&E bill
    grid_connection_customer_service = 0.56; # $/day
    distribution_service = 0.050; # $/kWh
    electricity_service = 0.098; # $/kWh
    state_tax = 5; # % on total bill
    county_tax = 0.5; # % on total bill
    # monthly_fixed_cost = (grid_connection_customer_service*30*(100 + state_tax + county_tax)/100); # $
    monthly_fixed_cost = 0; # We are keeping fixed cost zero for now
    elec_rate = (distribution_service + electricity_service)*(100 + state_tax + county_tax)/100; # $/kWh
    monthly_variable_cost = elec_rate * monthly_energy_kWh; # $
    disconnection_cost = 25; # $
    reconnection_cost = 25; # $
    # recharge amount for per month 
    # (this can be split into multiple small recharges according to real_recharge_logic)
    monthly_recharge_amount = (recharge_percent/100)*(monthly_fixed_cost .+ monthly_variable_cost);

    # Set real wallet recharge schedules
    if real_recharge_logic == 1
        monthly_real_recharge_schedule = [1 1];    
    elseif real_recharge_logic == 2
        monthly_real_recharge_schedule = [1 0.5; 15 0.5];
    elseif real_recharge_logic == 21
        monthly_real_recharge_schedule = [1 0.75; 22 0.25];
    elseif real_recharge_logic == 3
        monthly_real_recharge_schedule = [1 0.33; 11 0.33; 21 0.34];
    elseif real_recharge_logic == 4
        monthly_real_recharge_schedule = [1 0.25; 8 0.25; 15 0.25; 22 0.25];
    elseif real_recharge_logic == 5
        monthly_real_recharge_schedule = [1 0.2; 7 0.2; 13 0.2; 19 0.2; 25 0.2];
    elseif real_recharge_logic == 6
        monthly_real_recharge_schedule = [1 0.16; 6 0.16; 11 0.17; 16 0.17; 21 0.17; 26 0.17];
    elseif real_recharge_logic == 7
        monthly_real_recharge_schedule = [1 0.14; 5 0.14; 9 0.14; 13 0.14; 17 0.14; 21 0.15; 25 0.15];
    else
        monthly_real_recharge_schedule = [1 0];
    end
    real_recharge_schedule = zeros(simulation_time_months*size(monthly_real_recharge_schedule)[1],size(monthly_real_recharge_schedule)[2]);
    for i = 1:simulation_time_months
        row_indices = (i-1)*size(monthly_real_recharge_schedule)[1]+1:i*size(monthly_real_recharge_schedule)[1];
        real_recharge_schedule[row_indices,:] = [30*(i-1).+monthly_real_recharge_schedule[:,1]  monthly_recharge_amount[i]*monthly_real_recharge_schedule[:,2]];
    end

    # Fixed cost deduction schedule
    month_end_days = zeros(Int, simulation_time_months, 1);
    for i in 1:simulation_time_months
        month_end_days[i] = floor(Int, 30*i);
    end
    real_cost_schedule = [month_end_days monthly_fixed_cost*ones(simulation_time_months,1)];

    # Computing daily recharge and cost schedules
    daily_real_recharge = zeros(simulation_time_months*30);
    daily_real_cost = zeros(simulation_time_months*30);
    for iDay = 1:simulation_time_days
        if iDay in real_recharge_schedule[:,1]
            real_recharge_schedule_index = findall(x->x== iDay, real_recharge_schedule[:,1]);
            daily_real_recharge[iDay] = real_recharge_schedule[real_recharge_schedule_index[1],2];
        else
            daily_real_recharge[iDay] = 0;
        end
        if iDay in real_cost_schedule[:,1]
            real_cost_schedule_index = findall(x->x== iDay, real_cost_schedule[:,1]);
            daily_real_cost[iDay] = real_cost_schedule[real_cost_schedule_index[1], 2];
        else
            daily_real_cost[iDay] = 0;
        end
    end    

    return daily_real_recharge, daily_real_cost, elec_rate, real_recharge_schedule, real_cost_schedule, monthly_recharge_amount

end

function compute_daily_virtual_update(number_of_days, recharge_cost_schedule)
    daily_update = zeros(number_of_days)
    if size(recharge_cost_schedule)[1] == 1
        # if there is only one recharge/deduction in the duration and it is not on day 1, then the daily update till the day of the recharge/deduction is zero
        if recharge_cost_schedule[1,1] != 1
            daily_update[1:floor(Int,recharge_cost_schedule[1,1][1])-1] .= 0
        end
        # the daily update on and after that day till the end is the recharge/deduction divided by the number of days till the end
        daily_update[floor(Int,recharge_cost_schedule[1,1][1]):number_of_days] .= recharge_cost_schedule[1,2][1]/(number_of_days - recharge_cost_schedule[1,1][1] + 1)
    else
        # if there are multiple recharges/deductions in the duration and the first one is not on day 1, then the daily update till that day is zero
        if recharge_cost_schedule[1,1] != 1
            daily_update[1:floor(Int,recharge_cost_schedule[1,1][1])-1] .= 0     
        end
        # the daily update on and after the first recharge/deduction day is the latest recharge/deduction divided by the number of days till the next one
        for i in 1:floor(Int,size(recharge_cost_schedule)[1]) - 1
            daily_update[floor(Int,recharge_cost_schedule[i,1][1]):floor(Int,recharge_cost_schedule[i+1,1][1])-1] .= (recharge_cost_schedule[i,2][1])/(recharge_cost_schedule[i+1,1][1] - recharge_cost_schedule[i,1][1])
        end
        # the daily update on and after the last recharge/deduction day is the last recharge/deduction divided by the number of days till the end
        i = floor(Int,size(recharge_cost_schedule)[1])
        daily_update[floor(Int,recharge_cost_schedule[i,1][1]):number_of_days] .= recharge_cost_schedule[i,2][1]/(number_of_days - recharge_cost_schedule[i,1][1] + 1)
    end
    return daily_update
end

function compute_thresholds(time_data, demand_data, real_wallet_data, virtual_wallet_data)
    model = Model(Gurobi.Optimizer)
    number_of_timesteps = floor(Int,time_data["number_of_days"]*24/time_data["timestep_hours"]);
    timestep_hours = time_data["timestep_hours"];
    day_start_timesteps = zeros(Int,time_data["number_of_days"],1);
    load_priority_factor = demand_data["load_priority_factor"];
    for i in 1:time_data["number_of_days"]
        day_start_timesteps[i] = floor(Int,(24/timestep_hours)*(i-1) + 1);
    end
    # Add Variables
    @variable(model, 0 <= actuation_state[k in 1:demand_data["number_of_loads"], t in 1:number_of_timesteps] <= 1, Int);
    @variable(model, real_wallet[t in 1:(number_of_timesteps+1)]);
    @variable(model, 0 <= real_enable[k in 1:demand_data["number_of_loads"], t in 1:(number_of_timesteps+1)] <= 1, Int);
    @variable(model, 0 <= virtual_enable[k in 1:demand_data["number_of_loads"], t in 1:number_of_timesteps] <= 1, Int); 
    @variable(model, virtual_wallet[t in 1:number_of_timesteps]);
    @variable(model, virtual_wallet_thresholds[k in 1:demand_data["number_of_loads"], i in 1:time_data["number_of_days"]]);    

    # Add Constraints
    m = -10000;
    M = 10000;
    epsilon = 1e-6;
    # Real Wallet update constraint
    real_recharge_schedule_index = 0
    @constraint(model, real_wallet[1] == real_wallet_data["real_wallet_balance"] + real_wallet_data["daily_real_recharge"][1] - real_wallet_data["daily_real_cost"][1]);
    for t in 2:number_of_timesteps+1
        tm1 = t-1;
        recharge = 0;
        cost = 0;
        if t in day_start_timesteps
            recharge = real_wallet_data["daily_real_recharge"][floor(Int,1 + (t - 1)*timestep_hours/24)];
            cost = real_wallet_data["daily_real_cost"][floor(Int,1 + (t - 1)*timestep_hours/24)];
        end
        @constraint(model, real_wallet[t] == real_wallet[tm1] + (timestep_hours*real_wallet_data["cost_perkWh"]/1000)*(- sum(demand_data["demand_power_W"][k,tm1]*actuation_state[k,tm1] for k in 1:demand_data["number_of_loads"])) + recharge - cost);
    end
    # Real wallet enable signal constraints
    for t in 1:number_of_timesteps+1
        for k in 1:demand_data["number_of_loads"]
            @constraint(model, m*real_enable[k,t] <= -real_wallet[t]);
            @constraint(model, (M+epsilon)*(1-real_enable[k,t]) >= epsilon-real_wallet[t]);
        end
    end
    # Virtual wallet update constraints
    @constraint(model, virtual_wallet[1] == virtual_wallet_data["virtual_wallet_balance"] + virtual_wallet_data["daily_virtual_recharge"][1] - virtual_wallet_data["daily_virtual_cost"][1]);
    for iDays = 1:time_data["number_of_days"]
        for t in round(Int,(iDays-1)*24/timestep_hours + 2):round(Int,(iDays)*24/timestep_hours)
            tm1 = t-1;
            @constraint(model, virtual_wallet[t] == virtual_wallet[tm1] + (timestep_hours*real_wallet_data["cost_perkWh"]/1000)*(- sum(demand_data["demand_power_W"][k,tm1]*actuation_state[k,tm1] for k in 1:demand_data["number_of_loads"])));
        end
        if iDays >= 2
            tnp1 = round(Int, 1 + (iDays-1)*24/timestep_hours);
            tnp = round(Int, (iDays-1)*24/timestep_hours);
            @constraint(model, virtual_wallet[tnp1] == virtual_wallet[tnp] + virtual_wallet_data["daily_virtual_recharge"][iDays] - virtual_wallet_data["daily_virtual_cost"][iDays] + (timestep_hours*real_wallet_data["cost_perkWh"]/1000)*(- sum(demand_data["demand_power_W"][k,tnp]*actuation_state[k,tnp] for k in 1:demand_data["number_of_loads"])));
        end
    end
    # Virtual wallet enable constraints
    for t in 1:number_of_timesteps
        if mod(t,24/timestep_hours) == 0
            i = floor(Int,div(t,24/timestep_hours));
        else
            i = floor(Int,div(t,24/timestep_hours)) + 1;
        end
        for k in 1:demand_data["number_of_loads"]
            @constraint(model, virtual_wallet[t] - virtual_wallet_thresholds[k,i] - (M + epsilon)*virtual_enable[k,t] <= -epsilon);
            @constraint(model, virtual_wallet[t] - virtual_wallet_thresholds[k,i] -m*(1 - virtual_enable[k,t]) >= 0);
        end
    end
    # Actuation state constraints
    for t in 1:number_of_timesteps    
        for k in 1:demand_data["number_of_loads"]
            @constraint(model, actuation_state[k,t] <= demand_data["demand_state"][k,t]*virtual_enable[k,t]);
            @constraint(model, actuation_state[k,t] <= real_enable[k,t]);
            tp1 = t + 1;
            @constraint(model, actuation_state[k,t] <= real_enable[k,tp1])
            @constraint(model, demand_data["demand_state"][k,t]*virtual_enable[k,t] + real_enable[k,t] + real_enable[k,tp1] <= 2 + actuation_state[k,t]);
        end
    end  

    # Constraints for when thresholds have to be the same across more than 1 day
    if virtual_wallet_data["threshold_constant_days"] > 1
        for k in 1:demand_data["number_of_loads"]
            for j = 1:virtual_wallet_data["threshold_constant_days"]:time_data["number_of_days"]
                for i = j:j+(virtual_wallet_data["threshold_constant_days"]-2)
                    if i < time_data["number_of_days"]
                        @constraint(model, virtual_wallet_thresholds[k,i] == virtual_wallet_thresholds[k,i+1])
                    end
                end
            end
        end
    end 

    # Objective function
    @objective(model, Max, sum(load_priority_factor[k]*(sum(actuation_state[k,t] for t in 1:number_of_timesteps)/(sum(demand_data["demand_state"][k,t] for t in 1:number_of_timesteps)+0.0001)) for k in 1:demand_data["number_of_loads"]))

    # Solve model
    #set_optimizer_attribute(model, "LogToConsole", 0)                        
    optimize!(model)
    solve_time_model = solve_time(model)

    # Extract thresholds for days in update window
    update_number_of_timesteps = floor(Int,time_data["update_time_days"]*24/timestep_hours);
    virtual_wallet_thresholds_updatewindow = zeros(demand_data["number_of_loads"],update_number_of_timesteps);
    for t in 1:update_number_of_timesteps
        if mod(t,24/timestep_hours) == 0
            day_index = floor(Int,div(t,24/timestep_hours));
        else
            day_index = floor(Int,div(t,24/timestep_hours))+1;
        end
        virtual_wallet_thresholds_updatewindow[:,t] = value.(virtual_wallet_thresholds[:,day_index])
    end     

    return virtual_wallet_thresholds_updatewindow, value.(real_wallet[1:update_number_of_timesteps]), value.(virtual_wallet[1:update_number_of_timesteps]), value.(virtual_enable[:,1:update_number_of_timesteps]), value.(real_enable[:,1:update_number_of_timesteps]), value.(actuation_state[:,1:update_number_of_timesteps]), demand_data["demand_state"][:,1:update_number_of_timesteps], value.(virtual_wallet_thresholds[:,1:time_data["update_time_days"]]), solve_time_model

end

function compute_virtualrecharge_thresholds(time_data, demand_data, real_wallet_data, virtual_wallet_data, directory)
    model = Model(Gurobi.Optimizer)
    number_of_loads = demand_data["number_of_loads"];
    number_of_days = time_data["number_of_days"];
    number_of_timesteps = floor(Int,time_data["number_of_days"]*24/time_data["timestep_hours"]);
    timestep_hours = time_data["timestep_hours"];
    update_number_of_timesteps = floor(Int,time_data["update_time_days"]*24/timestep_hours);    
    day_start_timesteps = zeros(Int,time_data["number_of_days"],1);
    for i in 1:time_data["number_of_days"]
        day_start_timesteps[i] = floor(Int,(24/timestep_hours)*(i-1) + 1);
    end    
    load_priority_factor = demand_data["load_priority_factor"];
    demand_avg_power_W = zeros(number_of_loads,number_of_days);
    enable_duration_max = zeros(number_of_loads,number_of_days); 
    # max duration for which loads are to be enabled 
    # = 24h if average power demand is non-zero on that day, else = 0
    for d = 1:number_of_days
        timestep = day_start_timesteps[d];
        demand_avg_power_W[:,d] = demand_data["demand_power_W"][:,timestep];
        for k = 1:number_of_loads
            if demand_avg_power_W[k,d] > 0
                enable_duration_max[k,d] = 24;
            end
        end
    end
    daily_virtual_recharge_updatewindow = zeros(1,time_data["update_time_days"])
    daily_real_recharge = real_wallet_data["daily_real_recharge"]
    recharge_days = findall(x -> x!=0, daily_real_recharge)
    
    # Add Variables
    @variable(model, 0 <= enable_duration[k in 1:number_of_loads, d in 1:number_of_days]);
    @variable(model, daily_virtual_recharge[d in 1:number_of_days]);

    # Add Constraints
    for d = 1:number_of_days
        # upper bound on enable duration
        for k = 1:number_of_loads
            @constraint(model, enable_duration[k,d] <= enable_duration_max[k,d])
        end
        # upper bound on energy consumption and lower bound on virtual recharge
        if d == 1
            @constraint(model, sum(enable_duration[k,d]*demand_avg_power_W[k,d]*real_wallet_data["cost_perkWh"]/1000 for k in 1:number_of_loads) <= daily_virtual_recharge[d] + virtual_wallet_data["virtual_wallet_balance"])
            @constraint(model, daily_virtual_recharge[d] + virtual_wallet_data["virtual_wallet_balance"] >= 0)
        else            
            @constraint(model, sum(enable_duration[k,d]*demand_avg_power_W[k,d]*real_wallet_data["cost_perkWh"]/1000 for k in 1:number_of_loads) <= daily_virtual_recharge[d])
            @constraint(model, daily_virtual_recharge[d] >= 0)
        end
    end
    # upper bound on virtual recharge
    tolerance = 1e-4
    @constraint(model, virtual_wallet_data["virtual_wallet_balance"] + sum(daily_virtual_recharge[d] for d in 1:number_of_days) <= real_wallet_data["real_wallet_balance"] + sum(daily_real_recharge) - tolerance)
    # the following for loop can be omitted if recharge frequency is 1
    # the constraints in the loop assume that there will be a real wallet recharge on day 1
    for i = 1:length(recharge_days)-1
        @constraint(model, virtual_wallet_data["virtual_wallet_balance"] + sum(daily_virtual_recharge[d] for d in 1:recharge_days[i+1]-1) <= real_wallet_data["real_wallet_balance"] + sum(daily_real_recharge[recharge_days[j]] for j in 1:i) - tolerance)
    end

    # Objective function
    @objective(model, Max, sum(load_priority_factor[k]*(sum(enable_duration[k,d] for d in 1:number_of_days)/(sum(enable_duration_max[k,d] for d in 1:number_of_days)+0.0001)) for k in 1:number_of_loads))
    # Solve model
    optimize!(model)
    solve_time_model = solve_time(model)

    #model_psf = sum(load_priority_factor[k]*(sum(enable_duration[k,d] for d in 1:number_of_days)/(sum(enable_duration_max[k,d] for d in 1:number_of_days)+0.0001)) for k in 1:number_of_loads)

    days_array = 1:number_of_days
    groupedbar(reshape(value.(daily_virtual_recharge),length(value.(daily_virtual_recharge)),1),bar_position = :dodge, bar_width=0.7, title = "recharge", titlefontsize = 8, legend=false)
    virtual_wallet_thresholds = zeros(number_of_loads,number_of_days)
    for d = 1:number_of_days
        # value larger than virtual wallet balance at any point during the day
        virtual_wallet_thresholds[:,d] = 2*(virtual_wallet_data["virtual_wallet_balance"] + value.(daily_virtual_recharge[d]))*ones(number_of_loads,1)
    end
    # For debugging
    energy_used = 0;
    balance_added = 0;

    for d = 1:number_of_days
        # find indices of elements with value zero in enable_duration
        zero_element_indices = findall(x -> x==0, vec(value.(enable_duration[:,d])))
        # find indices of elements in enable_duration such that the values are arranged in ascending order
        ascending_enable_duration_indices = sortperm(value.(enable_duration[:,d]));
        # upper bound for index variable since length of ascending_enable_duration_indices will change in for loops
        j_ub = length(ascending_enable_duration_indices)
        # delete indices corresponding to elements with value zero
        for i =  1:length(zero_element_indices)
            for j = 1:j_ub
                if j > length(ascending_enable_duration_indices)
                    break
                elseif ascending_enable_duration_indices[j] == zero_element_indices[i]
                    ascending_enable_duration_indices = deleteat!(ascending_enable_duration_indices,j)
                end
            end
        end

        # determine thresholds for each day
        i = 0;
        for k in ascending_enable_duration_indices
            i = i + 1;
            if i == 1
                virtual_wallet_thresholds[k,d] = virtual_wallet_data["virtual_wallet_balance"] + value.(daily_virtual_recharge[d]) - value.(enable_duration[k,d])*(real_wallet_data["cost_perkWh"]/1000)*sum(demand_avg_power_W[k,d] for k in ascending_enable_duration_indices[i:length(ascending_enable_duration_indices)])
            else
                virtual_wallet_thresholds[k,d] = virtual_wallet_thresholds[ascending_enable_duration_indices[i-1],d] - (value.(enable_duration[k,d]) - value.(enable_duration[ascending_enable_duration_indices[i-1],d]))*(real_wallet_data["cost_perkWh"]/1000)*sum(demand_avg_power_W[k,d] for k in ascending_enable_duration_indices[i:length(ascending_enable_duration_indices)])
            end
        end        
    end
    # extract virtual recharge and thresholds for update window
    virtual_wallet_thresholds_updatewindow = zeros(number_of_loads,update_number_of_timesteps)
    daily_virtual_recharge_updatewindow = value.(daily_virtual_recharge[1:time_data["update_time_days"]]);

    for t in 1:update_number_of_timesteps
        if mod(t,24/timestep_hours) == 0
            day_index = floor(Int,div(t,24/timestep_hours));
        else
            day_index = floor(Int,div(t,24/timestep_hours))+1;
        end
        virtual_wallet_thresholds_updatewindow[:,t] = value.(virtual_wallet_thresholds[:,day_index])
    end

    daily_load_sf = zeros(number_of_loads,time_data["update_time_days"])
    for d = 1:time_data["update_time_days"]
        for k = 1:number_of_loads
            if enable_duration_max[k,d] == 0
                daily_load_sf[k,d] = -1
            else
                daily_load_sf[k,d] = 100*value.(enable_duration[k,d])/enable_duration_max[k,d]
            end
        end
    end
     
    return virtual_wallet_thresholds_updatewindow, daily_virtual_recharge_updatewindow, daily_load_sf, zeros(1:update_number_of_timesteps), zeros(1:update_number_of_timesteps), zeros(1:number_of_loads,1:update_number_of_timesteps), zeros(1:number_of_loads,1:update_number_of_timesteps), zeros(1:number_of_loads,1:update_number_of_timesteps), demand_data["demand_state"][:,1:update_number_of_timesteps], virtual_wallet_thresholds[:,1:time_data["update_time_days"]], solve_time_model

end

function compute_actuation_state(number_of_loads,virtual_wallet, virtual_wallet_thresholds, real_wallet, demand_state, case)
    virtual_enable = zeros(number_of_loads,1);
    real_enable = zeros(number_of_loads,1);
    actuation_state = zeros(number_of_loads,1);
    for k in 1:number_of_loads
        if virtual_wallet >= virtual_wallet_thresholds[k,1]
            virtual_enable[k,1] = 1;
        else
            virtual_enable[k,1] = 0;
        end
    end
    if real_wallet > 0
        real_enable[:,1] = ones(number_of_loads,1);
    else
        real_enable[:,1] = zeros(number_of_loads,1);
    end
    if case == 0
        actuation_state[:,1] = real_enable[:,1].*demand_state[:,1];
    elseif (case == 1) || (case == 2)
        actuation_state[:,1] = virtual_enable[:,1].*real_enable[:,1].*demand_state[:,1];
    end
    return actuation_state, virtual_enable, real_enable
end

function calculate_distance(vector_a, vector_b)
    distance = sqrt(sum((vector_a .- vector_b).^2)/length(vector_a))
    return distance
end

function process_results(time_array, simulation_results, model_results, parameters, directory)
    case = parameters["case"]
    number_of_loads = parameters["number_of_loads"];
    number_of_timesteps = parameters["number_of_timesteps"];
    demand_power_W = parameters["demand_power_W"];
    # Plot results
    real_wallet_plot = plot(time_array, simulation_results["real_wallet"], label = "simulation");
    plot!(time_array, model_results["real_wallet"], label = "model");
    virtual_wallet_plot = plot(time_array, simulation_results["virtual_wallet"], label = "simulation");
    plot!(time_array, simulation_results["virtual_wallet_thresholds"]', label = ["air1" "clotheswasher_dryg1" "microwave1" "refrigerator1"], legend=:bottomleft)
    plot!(time_array, model_results["virtual_wallet"], label = "model");
    actuation_state_plot = plot(time_array,simulation_results["actuation_state"]',layout = (number_of_loads,1));
    demand_state_plot = plot(time_array, simulation_results["demand_state"]', layout = (number_of_loads,1));
    
    # Compute metrics
    # Cumulative SF
    time_service_factor = (sum(sum(simulation_results["actuation_state"][k,t] for k in 1:number_of_loads) for t in 1:number_of_timesteps))*100/(sum(sum(simulation_results["demand_state"][k,t] for k in 1:number_of_loads) for t in 1:number_of_timesteps)) 
    # Average SF
    avg_time_service_factor = 100*(1/number_of_loads)*sum((sum(simulation_results["actuation_state"][k,t] for t in 1:number_of_timesteps)/sum(simulation_results["demand_state"][k,t] for t in 1:number_of_timesteps)) for k in 1:number_of_loads)
    # Energy SF
    energy_service_factor = (sum(sum(demand_power_W[k,t]*simulation_results["actuation_state"][k,t] for k in 1:number_of_loads) for t in 1:number_of_timesteps))*100/(sum(sum(demand_power_W[k,t]*simulation_results["demand_state"][k,t] for k in 1:number_of_loads) for t in 1:number_of_timesteps))     
    print("Time Service Factor = ")
    println(time_service_factor)
    print("Energy Service Factor = ")
    println(energy_service_factor)
    print("Average Time Service Factor = ")
    println(avg_time_service_factor)
    # Compute metrics using values from optimization model
    model_time_service_factor = (sum(sum(model_results["actuation_state"][k,t] for k in 1:number_of_loads) for t in 1:number_of_timesteps))*100/(sum(sum(model_results["demand_state"][k,t] for k in 1:number_of_loads) for t in 1:number_of_timesteps)) 
    model_energy_service_factor = (sum(sum(demand_power_W[k,t]*model_results["actuation_state"][k,t] for k in 1:number_of_loads) for t in 1:number_of_timesteps))*100/(sum(sum(demand_power_W[k,t]*model_results["demand_state"][k,t] for k in 1:number_of_loads) for t in 1:number_of_timesteps))     
    print("Model Time Service Factor = ")
    println(model_time_service_factor)
    print("Model Energy Service Factor = ")
    println(model_energy_service_factor)

    # Priority SF
    priorityweighted_time_service_factor = 100*sum(parameters["load_priority_factor"][k]*(sum(simulation_results["actuation_state"][k,t] for t in 1:number_of_timesteps)/sum(simulation_results["demand_state"][k,t] for t in 1:number_of_timesteps)) for k in 1:number_of_loads)
    print("Priority Weighted Time Service Factor = ")
    println(priorityweighted_time_service_factor)
    # Priority SF using values from optimization model
    model_priorityweighted__time_service_factor = 100*sum(parameters["load_priority_factor"][k]*(sum(model_results["actuation_state"][k,t] for t in 1:number_of_timesteps)/sum(model_results["demand_state"][k,t] for t in 1:number_of_timesteps)) for k in 1:number_of_loads)
    print("Priority Weighted Model Time Service Factor = ")
    println(model_priorityweighted__time_service_factor)

    # SF per load
    time_service_factor_perload = zeros(number_of_loads,1);
    model_time_service_factor_perload = zeros(number_of_loads,1);
    for k in 1:number_of_loads
        time_service_factor_perload[k] = (sum(simulation_results["actuation_state"][k,t] for t in 1:number_of_timesteps))*100/(sum(simulation_results["demand_state"][k,t] for t in 1:number_of_timesteps))     
        model_time_service_factor_perload[k] = (sum(model_results["actuation_state"][k,t] for t in 1:number_of_timesteps))*100/(sum(model_results["demand_state"][k,t] for t in 1:number_of_timesteps))     
    end
    print("Time Service Factor Per Load = ")
    println(time_service_factor_perload)
    print("Model Time Service Factor Per Load = ")
    println(model_time_service_factor_perload)

    return priorityweighted_time_service_factor, time_service_factor_perload', time_service_factor, energy_service_factor, model_priorityweighted__time_service_factor, model_time_service_factor_perload', model_time_service_factor, model_energy_service_factor

end    

function generate_modified_data(demand_power_W,demand_state,number_of_timesteps,horizon_time_days,simulation_timestep_hours,simulation_time_days, model_data_number, number_of_loads, directory)
    demand_power_W_modified = zeros(number_of_loads,number_of_timesteps);
    demand_state_modified = zeros(number_of_loads,number_of_timesteps);        
    if model_data_number == 1
        ## Average Power
        averaging_interval_timesteps = floor(Int,24/simulation_timestep_hours);
        for i = 1:averaging_interval_timesteps:number_of_timesteps
            ub = i+averaging_interval_timesteps-1;
            if ub > number_of_timesteps
                ub = number_of_timesteps
            end
            for k = 1:4
                demand_power_W_modified[k,i:ub] .= mean(demand_power_W[k,i:ub].*demand_state[k,i:ub])
                if mean(demand_power_W[k,i:ub].*demand_state[k,i:ub]) > 0
                    demand_state_modified[k,i:ub] .= 1;
                else
                    demand_state_modified[k,i:ub] .= 0;
                end
            end
        end
    elseif model_data_number == 2
        ## Random
        # Generated using c = shuffle(Vector(1:30));
        c = [9  23  5  27  14  30  17  15  2  13 16  7  26  28  12  4  25  3  6  8 20  1  11  10  19  29  24  21  18  22];
        demand_power_W_modified = zeros(4,number_of_timesteps);
        demand_state_modified = zeros(4,number_of_timesteps);
        for i = 1:length(c)
            demand_power_W_modified[:,floor(Int,(i-1)*(24/simulation_timestep_hours))+1:floor(Int,i*(24/simulation_timestep_hours))] = demand_power_W[:,floor(Int,(c[i]-1)*(24/simulation_timestep_hours))+1:floor(Int,c[i]*(24/simulation_timestep_hours))]
            demand_state_modified[:,floor(Int,(i-1)*(24/simulation_timestep_hours))+1:floor(Int,i*(24/simulation_timestep_hours))] = demand_state[:,floor(Int,(c[i]-1)*(24/simulation_timestep_hours))+1:floor(Int,c[i]*(24/simulation_timestep_hours))]
        end        
    elseif model_data_number == 2.1
        ## Random and averaged
        # Generated using c = shuffle(Vector(1:30));
        c = [9  23  5  27  14  30  17  15  2  13 16  7  26  28  12  4  25  3  6  8 20  1  11  10  19  29  24  21  18  22];
        demand_power_W_shuffled = zeros(4,number_of_timesteps);
        demand_state_shuffled = zeros(4,number_of_timesteps);
        for i = 1:length(c)
            demand_power_W_shuffled[:,floor(Int,(i-1)*(24/simulation_timestep_hours))+1:floor(Int,i*(24/simulation_timestep_hours))] = demand_power_W[:,floor(Int,(c[i]-1)*(24/simulation_timestep_hours))+1:floor(Int,c[i]*(24/simulation_timestep_hours))]
            demand_state_shuffled[:,floor(Int,(i-1)*(24/simulation_timestep_hours))+1:floor(Int,i*(24/simulation_timestep_hours))] = demand_state[:,floor(Int,(c[i]-1)*(24/simulation_timestep_hours))+1:floor(Int,c[i]*(24/simulation_timestep_hours))]
        end
        averaging_interval_timesteps = floor(Int,24/simulation_timestep_hours);
        for i = 1:averaging_interval_timesteps:number_of_timesteps
            ub = i+averaging_interval_timesteps-1;
            if ub > number_of_timesteps
                ub = number_of_timesteps
            end
            for k = 1:4
                demand_power_W_modified[k,i:ub] .= mean(demand_power_W_shuffled[k,i:ub].*demand_state_shuffled[k,i:ub])
                if mean(demand_power_W_shuffled[k,i:ub].*demand_state_shuffled[k,i:ub]) > 0
                    demand_state_modified[k,i:ub] .= 1;
                else
                    demand_state_modified[k,i:ub] .= 0;
                end
            end
        end            
    elseif model_data_number == 1.1
        ## Average Power
        averaging_interval_timesteps = floor(Int,24*30/simulation_timestep_hours);
        for i = 1:averaging_interval_timesteps:number_of_timesteps
            ub = i+averaging_interval_timesteps-1;
            if ub > number_of_timesteps
                ub = number_of_timesteps
            end
            for k = 1:4
                demand_power_W_modified[k,i:ub] .= mean(demand_power_W[k,i:ub].*demand_state[k,i:ub])
                if mean(demand_power_W[k,i:ub].*demand_state[k,i:ub]) > 0
                    demand_state_modified[k,i:ub] .= 1;
                else
                    demand_state_modified[k,i:ub] .= 0;
                end
            end
        end
    end
    
    lb = 1; ub = floor(Int,24/simulation_timestep_hours)*30;
    power_modified_plot = plot(lb:ub,demand_power_W_modified[:,lb:ub]',layout=(4,1))
    state_modified_plot = plot(lb:ub,demand_state_modified[:,lb:ub]',layout=(4,1))
    power_nonmodified_plot = plot(lb:ub,demand_power_W[:,lb:ub]',layout=(4,1))
    state_nonmodified_plot = plot(lb:ub,demand_state[:,lb:ub]',layout=(4,1))

    day_start_timesteps = zeros(Int, simulation_time_days,1);
    for i in 1:simulation_time_days
        day_start_timesteps[i] = floor(Int,(24/simulation_timestep_hours)*(i-1) + 1);
    end
    
    demand_energy_Wh_modified = zeros(number_of_loads,simulation_time_days)
    for i in 1:simulation_time_days
        for k in 1:number_of_loads
            demand_energy_Wh_modified[k,i] = sum(demand_power_W_modified[k,t]*demand_state_modified[k,t] for t in day_start_timesteps[i]:floor(Int,day_start_timesteps[i]+(24/simulation_timestep_hours)-1))*simulation_timestep_hours;
        end
    end
    daily_energy_array_Wh_modified = [demand_energy_Wh_modified[1,:] demand_energy_Wh_modified[2,:] demand_energy_Wh_modified[3,:] demand_energy_Wh_modified[4,:]]
    daily_energy_Wh_modified = zeros(simulation_time_days,1)
    for i in 1:simulation_time_days
        daily_energy_Wh_modified[i,1] = sum(demand_energy_Wh_modified[k,i] for k in 1:number_of_loads)
    end

    return demand_power_W_modified, demand_state_modified, daily_energy_array_Wh_modified, daily_energy_Wh_modified
end

function compute_hem_schedule(total_time_days, timestep_hours, demand_state, demand_power_W, cost_perkWh, number_of_loads, load_priority_factor, daily_real_recharge, directory)
    number_of_timesteps = floor(Int,total_time_days*24/timestep_hours)   
    model = Model(Gurobi.Optimizer)
    @variable(model, 0 <= actuation_state[k in 1:number_of_loads, t in 1:number_of_timesteps] <= 1, Int);
    for k in 1:number_of_loads
        for t in 1:number_of_timesteps
            @constraint(model, actuation_state[k,t] <= demand_state[k,t])
        end
    end
    recharge_days = findall(x -> x!=0, daily_real_recharge)
    recharge_timesteps = zeros(Int,length(recharge_days))
    for i = 1:length(recharge_timesteps)
        recharge_timesteps[i] = floor(Int,(recharge_days[i]-1)*24/timestep_hours + 1)
    end
    # upper bound on energy usage
    tolerance = 1e-4    
    @constraint(model,(cost_perkWh/1000)*timestep_hours*sum(sum(actuation_state[k,t]*demand_power_W[k,t] for k in 1:number_of_loads) for t in 1:number_of_timesteps) <= sum(daily_real_recharge[recharge_days[i]] for i in 1:length(recharge_days)) - tolerance)
    # the following for loop can be omitted if recharge frequency is 1
    # the constraints in the loop assume that there will be a real wallet recharge on day 1    
    for i = 1:length(recharge_days)-1
        @constraint(model, (cost_perkWh/1000)*timestep_hours*sum(sum(actuation_state[k,t]*demand_power_W[k,t] for k in 1:number_of_loads) for t in 1:recharge_timesteps[i+1]-1) <= sum(daily_real_recharge[recharge_days[j]] for j in 1:i) - tolerance)
    end
    
    @objective(model, Max, sum(load_priority_factor[k]*(sum(actuation_state[k,t] for t in 1:number_of_timesteps)/(sum(demand_state[k,t] for t in 1:number_of_timesteps)+0.0001)) for k in 1:number_of_loads))
    optimize!(model)
    solve_time_model = solve_time(model)

    load1_plot = plot(value.(actuation_state[1,:])); plot!(0.5*demand_state[1,:])
    load2_plot = plot(value.(actuation_state[2,:])); plot!(0.5*demand_state[2,:])
    load3_plot = plot(value.(actuation_state[3,:])); plot!(0.5*demand_state[3,:])
    load4_plot = plot(value.(actuation_state[4,:])); plot!(0.5*demand_state[4,:])
    demand_actuation_plot = plot(load1_plot,load2_plot,load3_plot,load4_plot,layout = (4,1))    

    priority_sf = 100*sum(load_priority_factor[k]*(sum(value.(actuation_state[k,t]) for t in 1:number_of_timesteps)/(sum(demand_state[k,t] for t in 1:number_of_timesteps))) for k in 1:number_of_loads)
    energy_sf = sum(sum(value.(actuation_state[k,t])*demand_power_W[k,t] for k in 1:number_of_loads) for t in 1:number_of_timesteps)*100/sum(sum(demand_state[k,t]*demand_power_W[k,t] for k in 1:number_of_loads) for t in 1:number_of_timesteps)
    # SF per load
    time_service_factor_perload = zeros(number_of_loads,1);
    for k in 1:number_of_loads
        time_service_factor_perload[k] = (sum(value.(actuation_state[k,t]) for t in 1:number_of_timesteps))*100/(sum(demand_state[k,t] for t in 1:number_of_timesteps))     
    end
    return value.(actuation_state), priority_sf, energy_sf, time_service_factor_perload, solve_time_model
end


function compute_thresholds_extra_constraints(time_data, demand_data, real_wallet_data, virtual_wallet_data)
    model = Model(Gurobi.Optimizer)
    number_of_timesteps = floor(Int,time_data["number_of_days"]*24/time_data["timestep_hours"]);
    timestep_hours = time_data["timestep_hours"];
    day_start_timesteps = zeros(Int,time_data["number_of_days"],1);
    load_priority_factor = demand_data["load_priority_factor"];
    for i in 1:time_data["number_of_days"]
        day_start_timesteps[i] = floor(Int,(24/timestep_hours)*(i-1) + 1);
    end
    # Add Variables
    @variable(model, 0 <= actuation_state[k in 1:demand_data["number_of_loads"], t in 1:number_of_timesteps] <= 1, Int);
    @variable(model, real_wallet[t in 1:(number_of_timesteps+1)]);
    @variable(model, 0 <= real_enable[k in 1:demand_data["number_of_loads"], t in 1:(number_of_timesteps+1)] <= 1, Int);
    @variable(model, 0 <= virtual_enable[k in 1:demand_data["number_of_loads"], t in 1:number_of_timesteps] <= 1, Int); 
    @variable(model, 0 <= virtual_wallet[t in 1:number_of_timesteps]);
    @variable(model, virtual_wallet_thresholds[k in 1:demand_data["number_of_loads"], i in 1:time_data["number_of_days"]]);    
    @variable(model, daily_virtual_recharge[d in 1:time_data["number_of_days"]])

    # Add Constraints
    m = -10000;
    M = 10000;
    epsilon = 1e-6;
    # Real Wallet update constraint
    real_recharge_schedule_index = 0
    @constraint(model, real_wallet[1] == real_wallet_data["real_wallet_balance"] + real_wallet_data["daily_real_recharge"][1] - real_wallet_data["daily_real_cost"][1]);
    for t in 2:number_of_timesteps+1
        tm1 = t-1;
        recharge = 0;
        cost = 0;
        if t in day_start_timesteps
            recharge = real_wallet_data["daily_real_recharge"][floor(Int,1 + (t - 1)*timestep_hours/24)];
            cost = real_wallet_data["daily_real_cost"][floor(Int,1 + (t - 1)*timestep_hours/24)];
        end
        @constraint(model, real_wallet[t] == real_wallet[tm1] + (timestep_hours*real_wallet_data["cost_perkWh"]/1000)*(- sum(demand_data["demand_power_W"][k,tm1]*actuation_state[k,tm1] for k in 1:demand_data["number_of_loads"])) + recharge - cost);
    end
    # Real wallet enable signal constraints
    for t in 1:number_of_timesteps+1
        for k in 1:demand_data["number_of_loads"]
            @constraint(model, m*real_enable[k,t] <= -real_wallet[t]);
            @constraint(model, (M+epsilon)*(1-real_enable[k,t]) >= epsilon-real_wallet[t]);
        end
    end
    # Virtual wallet update constraints
    @constraint(model, virtual_wallet[1] == virtual_wallet_data["virtual_wallet_balance"] + daily_virtual_recharge[1]);
    for iDays = 1:time_data["number_of_days"]
        for t in round(Int,(iDays-1)*24/timestep_hours + 2):round(Int,(iDays)*24/timestep_hours)
            tm1 = t-1;
            @constraint(model, virtual_wallet[t] == virtual_wallet[tm1] + (timestep_hours*real_wallet_data["cost_perkWh"]/1000)*(- sum(demand_data["demand_power_W"][k,tm1]*actuation_state[k,tm1] for k in 1:demand_data["number_of_loads"])));
        end
        if iDays >= 2
            tnp1 = round(Int, 1 + (iDays-1)*24/timestep_hours);
            tnp = round(Int, (iDays-1)*24/timestep_hours);
            @constraint(model, virtual_wallet[tnp1] == virtual_wallet[tnp] + daily_virtual_recharge[iDays] + (timestep_hours*real_wallet_data["cost_perkWh"]/1000)*(- sum(demand_data["demand_power_W"][k,tnp]*actuation_state[k,tnp] for k in 1:demand_data["number_of_loads"])));
        end
    end
    recharge_days = findall(x -> x!=0, real_wallet_data["daily_real_recharge"])    
    # upper bound on virtual wallet recharge
    @constraint(model, virtual_wallet_data["virtual_wallet_balance"] + sum(daily_virtual_recharge[d] for d in 1:time_data["number_of_days"]) <= real_wallet_data["real_wallet_balance"] + sum(real_wallet_data["daily_real_recharge"]))
    # the following for loop can be omitted if recharge frequency is 1
    # the constraints in the loop assume that there will be a real wallet recharge on day 1    
    for i = 1:length(recharge_days)-1
        @constraint(model, virtual_wallet_data["virtual_wallet_balance"] + sum(daily_virtual_recharge[d] for d in 1:recharge_days[i+1]-1) <= real_wallet_data["real_wallet_balance"] + sum(real_wallet_data["daily_real_recharge"][recharge_days[j]] for j in 1:i))
    end



    # Virtual wallet enable constraints
    for t in 1:number_of_timesteps
        if mod(t,24/timestep_hours) == 0
            i = floor(Int,div(t,24/timestep_hours));
        else
            i = floor(Int,div(t,24/timestep_hours)) + 1;
        end
        for k in 1:demand_data["number_of_loads"]
            @constraint(model, virtual_wallet[t] - virtual_wallet_thresholds[k,i] - (M + epsilon)*virtual_enable[k,t] <= -epsilon);
            @constraint(model, virtual_wallet[t] - virtual_wallet_thresholds[k,i] -m*(1 - virtual_enable[k,t]) >= 0);
        end
    end
    # Actuation state constraints
    for t in 1:number_of_timesteps    
        for k in 1:demand_data["number_of_loads"]
            @constraint(model, actuation_state[k,t] <= demand_data["demand_state"][k,t]*virtual_enable[k,t]);
            @constraint(model, actuation_state[k,t] <= real_enable[k,t]);
            tp1 = t + 1;
            @constraint(model, actuation_state[k,t] <= real_enable[k,tp1])
            @constraint(model, demand_data["demand_state"][k,t]*virtual_enable[k,t] + real_enable[k,t] + real_enable[k,tp1] <= 2 + actuation_state[k,t]);
        end
    end  

    # Constraints for when thresholds have to be the same across more than 1 day
    if virtual_wallet_data["threshold_constant_days"] > 1
        for k in 1:demand_data["number_of_loads"]
            for j = 1:virtual_wallet_data["threshold_constant_days"]:time_data["number_of_days"]
                for i = j:j+(virtual_wallet_data["threshold_constant_days"]-2)
                    if i < time_data["number_of_days"]
                        @constraint(model, virtual_wallet_thresholds[k,i] == virtual_wallet_thresholds[k,i+1])
                    end
                end
            end
        end
    end 

    # Objective function
    @objective(model, Max, sum(load_priority_factor[k]*(sum(actuation_state[k,t] for t in 1:number_of_timesteps)/(sum(demand_data["demand_state"][k,t] for t in 1:number_of_timesteps)+0.0001)) for k in 1:demand_data["number_of_loads"]))

    # Solve model
    #set_optimizer_attribute(model, "LogToConsole", 0)                        
    optimize!(model)

    # Extract thresholds for days in update window
    update_number_of_timesteps = floor(Int,time_data["update_time_days"]*24/timestep_hours);
    virtual_wallet_thresholds_updatewindow = zeros(demand_data["number_of_loads"],update_number_of_timesteps);
    for t in 1:update_number_of_timesteps
        if mod(t,24/timestep_hours) == 0
            day_index = floor(Int,div(t,24/timestep_hours));
        else
            day_index = floor(Int,div(t,24/timestep_hours))+1;
        end
        virtual_wallet_thresholds_updatewindow[:,t] = value.(virtual_wallet_thresholds[:,day_index])
    end     
    daily_virtual_recharge_updatewindow = value.(daily_virtual_recharge[1:time_data["update_time_days"]]);
    return virtual_wallet_thresholds_updatewindow, value.(real_wallet[1:update_number_of_timesteps]), value.(virtual_wallet[1:update_number_of_timesteps]), value.(virtual_enable[:,1:update_number_of_timesteps]), value.(real_enable[:,1:update_number_of_timesteps]), value.(actuation_state[:,1:update_number_of_timesteps]), demand_data["demand_state"][:,1:update_number_of_timesteps], value.(virtual_wallet_thresholds[:,1:time_data["update_time_days"]]), daily_virtual_recharge_updatewindow

end

