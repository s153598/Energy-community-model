using CSV, DataFrames
using JuMP, Gurobi
import XLSX
using Plots
using StatsPlots

function ModelDataImport()

    # Sets of the model
    global Sets = DataFrame(XLSX.readtable("Data/Input_Data_Tariff_Model.xlsx", "Sets")...)

    # Demand profile of household types
    global Demand_profiles = DataFrame(XLSX.readtable("Data/El_Profiles.xlsx", "El_Demand")...)
    global Demand_profiles2 = DataFrame(XLSX.readtable("Data/El_Profiles.xlsx", "El_Demand2")...)
    global Heat_profiles = DataFrame(XLSX.readtable("Data/Heat_Profiles.xlsx", "Heat_Demand")...)
    global Heat_profiles2 = DataFrame(XLSX.readtable("Data/Heat_Profiles.xlsx", "Heat_Demand2")...)

    # Hourly electricity profiles
    #global El_price = CSV.read("Data/Electricity_prices.csv",DataFrame)
    global El_price = DataFrame(XLSX.readtable("Data/ELprice.xlsx","Sheet1")...)

    # Network tariffs including distribution tariffs, PSO and energy tax
    global Network_tariffs  = CSV.File("Data/Network_tariffs.csv") |> Dict

    # Technical and economic parameters of PV array
    global PV_par = CSV.File("Data/PV_par.csv") |> Dict

    # Technical and economic parameters of battery
    global Battery_par = CSV.File("Data/Battery_par.csv") |> Dict
    global Grid_par = CSV.File("Data/Grid_par.csv") |> Dict
    global Scalars = CSV.File("Data/Scalars.csv") |> Dict

    # Capacity factor of PV Array
    global PV_CF = DataFrame(XLSX.readtable("Data/Input_Data_Tariff_Model.xlsx", "SolarCF")...)

    # Parameter defining the technologies in each technology type
    global Household_types  = DataFrame(XLSX.readtable("Data/Input_Data_Tariff_Model.xlsx", "Household_types")...)

    # EV parameters
    global EV_DF = DataFrame(XLSX.readtable("Data/EV_Profiles.xlsx", "EV_Avail")...)
    global EV_par = CSV.File("Data/EV_par.csv") |> Dict
    global EV_cons = DataFrame(XLSX.readtable("Data/EV_Profiles.xlsx", "EV_Cons")...)

    # Technical and economic parameters of P2H
    global P2H_par = CSV.File("Data/HP_par.csv") |> Dict

    #CO2 intentities
    global EL_CO2 = DataFrame(XLSX.readtable("Data/Intensities.xlsx", "Intensities (kg per kWh)")...)

    #Boiler
    global Boiler_par = CSV.File("Data/Boiler_par.csv") |> Dict

    #Heat Storage
    global Heat_Storage_par = CSV.File("Data/Heat_Storage_par.csv") |> Dict

    #Energy community
    global Ec_par = CSV.File("Data/Ec_par.csv") |> Dict

    # Assigning sets
    global T = Sets[:,"T"]
    #global Y = collect(skipmissing(Sets[:,"Y"]))
    global S = collect(skipmissing(Sets[:,"S"]))
    global H = collect(skipmissing(Sets[:,"Household"]))

end

function InitializeModel(scheme)
    M = Model(Gurobi.Optimizer)

    #Battery variables.
    @variable(M,C_BT[s=S],lower_bound=0, upper_bound=50 ,base_name="C_BT[s=S]:") #Capacity
    @variable(M,b_st[t=T,s=S],lower_bound=0,     base_name="b_st[t=T,s=S]:") #Variable holding the state of charging
    @variable(M,b_dh[t=T,s=S],lower_bound=0,     base_name="b_dh[t=T,s=S]:") #Showing the amount discharged
    @variable(M,b_ch[t=T,s=S],lower_bound=0,base_name="b_ch[t=T,s=S]:") #Showing the charged amount
    @variable(M,b_dh_load[t=T,s=S],lower_bound=0,base_name="b_dh_load[t=T,s=S]:")
    @variable(M,b_dh_ex[t=T,s=S],lower_bound=0,base_name="b_dh_ex[t=T,s=S]:")
    @variable(M,b_dh_ev[t=T,s=S],lower_bound=0,base_name="b_dh_ev[t=T,s=S]:")
    @variable(M,b_dh_p2h[t=T,s=S],lower_bound=0,base_name="b_dh_p2h[t=T,s=S]:")
    @variable(M,b_dh_eb[t=T,s=S],lower_bound=0,base_name="b_dh_eb[t=T,s=S]:")
    @variable(M,b_dh_ec[t=T,s=S],lower_bound=0,base_name="b_dh_ec[t=T,s=S]:")



    # EV related variables
    @variable(M,C_EV[h=H,s=S],lower_bound=0,base_name="C_EV[h=H,s=S]:")
    @variable(M,ev_st[t=T,h=H,s=S],lower_bound=0,base_name="ev_st[t=T,h=H,s=S]:")
    @variable(M,ev_dh[t=T,h=H,s=S],lower_bound=0,base_name="ev_dh[t=T,h=H,s=S]:")
    @variable(M,ev_ch[t=T,h=H,s=S],lower_bound=0,base_name="ev_ch[t=T,h=H,s=S]:")
    @variable(M,ev_dh_load[t=T,h=H,s=S],lower_bound=0,base_name="ev_dh_load[t=T,h=H,s=S]:")
    @variable(M,ev_dh_ex[t=T,h=H,s=S],lower_bound=0,base_name="ev_dh_ex[t=T,h=H,s=S]:")
    @variable(M,ev_dh_p2h[t=T,h=H,s=S],lower_bound=0,base_name="ev_dh_p2h[t=T,h=H,s=S]:")
    @variable(M,ev_dh_eb[t=T,h=H,s=S],lower_bound=0,base_name="ev_dh_eb[t=T,h=H,s=S]:")
    @variable(M,ev_dh_ec[t=T,h=H,s=S],lower_bound=0,base_name="ev_dh_ec[t=T,h=H,s=S]:")

    # PV related constraints
    @variable(M,C_PV[s=S],lower_bound=0,upper_bound=70, base_name="C_PV[s=S]:")
    @variable(M,p_PV[t=T,s=S],lower_bound=0,base_name="p_PV[t=T,s=S]:")
    @variable(M,p_PV_load[t=T,s=S],lower_bound=0,base_name="p_PV_load[t=T,s=S]:")
    @variable(M,p_PV_bat[t=T,s=S],lower_bound=0,base_name="p_PV_bat[t=T,s=S]:")
    @variable(M,p_PV_ev[t=T,s=S],lower_bound=0,base_name="p_PV_ev[t=T,s=S]:")
    @variable(M,p_PV_ex[t=T,s=S],lower_bound=0,base_name="p_PV_ex[t=T,s=S]:")
    @variable(M,p_PV_p2h[t=T,s=S],lower_bound=0,base_name="p_PV_p2h[t=T,s=S]:")
    @variable(M,p_PV_eb[t=T,s=S],lower_bound=0,base_name="p_PV_eb[t=T,s=S]:")
    @variable(M,p_PV_ec[t=T,s=S],lower_bound=0,base_name="p_PV_ec[t=T,s=S]:")

    #P2H related constraints
    @variable(M,C_P2H[s=S],lower_bound=0,upper_bound = 60, base_name="C_P2H[s=S]:")
    @variable(M,p2h_p[t=T,s=S],lower_bound=0,base_name="p2h_p[t=T,s=S]:")
    @variable(M,p2h_load[t=T,s=S],lower_bound=0,base_name="p2h_load[t=T,s=S]:")
    @variable(M,p2h_im[t=T,s=S],lower_bound=0,base_name="p2h_im[t=T,s=S]:")
    @variable(M,p2h_hs[t=T,s=S],lower_bound=0,base_name="p2h_hs[t=T,s=S]:")


    #Electric Boiler
    @variable(M,C_EB[s=S],lower_bound=0,base_name="C_EB[s=S]:")
    @variable(M,eb_p[t=T,s=S],lower_bound=0,base_name="eb_p[t=T,s=S]:")
    @variable(M,eb_load[t=T,s=S],lower_bound=0,base_name="eb_load[t=T,s=S]:")
    @variable(M,eb_im[t=T,s=S],lower_bound=0,base_name="eb_im[t=T,s=S]:")
    @variable(M,eb_hs[t=T,s=S],lower_bound=0,base_name="eb_hs[t=T,s=S]:")


    #Heat Storage Constaints
    @variable(M,C_HS[s=S],lower_bound=0, upper_bound=180 ,base_name="C_HS[s=S]:") #Capacity
    @variable(M,hs_st[t=T,s=S],lower_bound=0,     base_name="hs_st[t=T,s=S]:") #Variable holding the state of charging
    @variable(M,hs_dh[t=T,s=S],lower_bound=0,     base_name="hs_dh[t=T,s=S]:") #Showing the amount discharged
    @variable(M,hs_ch[t=T,s=S],lower_bound=0,base_name="hs_ch[t=T,s=S]:") #Showing the charged amount
    @variable(M,hs_dh_load[t=T,s=S],lower_bound=0,base_name="hs_dh_load[t=T,s=S]:") #Showing the charged amount

    # Grid related constraints
    @variable(M,g_ex[t=T,s=S],lower_bound=0,base_name="g_ex[t=T,s=S]:")
    @variable(M,g_im[t=T,s=S],lower_bound=0,base_name="g_im[t=T,s=S]:")
    @variable(M,g_im_load[t=T,s=S],lower_bound=0,base_name="g_im_load[t=T,s=S]:")
    @variable(M,g_im_bat[t=T,s=S],lower_bound=0,base_name="g_im_bat[t=T,s=S]:")
    @variable(M,g_im_ev[t=T,s=S],lower_bound=0,base_name="g_im_ev[t=T,s=S]:")
    @variable(M,g_im_p2h[t=T,s=S],lower_bound=0,base_name="g_im_p2h[t=T,s=S]:")
    @variable(M,g_im_eb[t=T,s=S],lower_bound=0,base_name="g_im_eb[t=T,s=S]:")

    #EC variables
    @variable(M,ec_ex[t=T,c=S,s=S],lower_bound=0,base_name="ec_ex[t=T,c=S,s=S]:")
    @variable(M,ec_im[t=T,c=S,s=S],lower_bound=0,base_name="ec_im[t=T,c=S,s=S]:")
    @variable(M,ec_im_load[t=T,s=S],lower_bound=0,base_name="ec_im_load[t=T,s=S]:")
    @variable(M,ec_im_bat[t=T,s=S],lower_bound=0,base_name="ec_im_bat[t=T,s=S]:")
    @variable(M,ec_im_ev[t=T,s=S],lower_bound=0,base_name="ec_im_ev[t=T,s=S]:")
    @variable(M,ec_im_p2h[t=T,s=S],lower_bound=0,base_name="ec_im_p2h[t=T,s=S]:")
    @variable(M,ec_im_eb[t=T,s=S],lower_bound=0,base_name="ec_im_eb[t=T,s=S]:")
    if scheme == "new"
    @variable(M,ec_im_var[t=T,s=S],Bin,base_name="ec_im_var[t=T,s=S]:")
    elseif scheme == "base"
    end
    return(M)
end

# -------------------------------------------------------------------------------------------------------------------
#                                        FIXING THE VARIABLES FUNCTION
# -------------------------------------------------------------------------------------------------------------------

# Fixing the capacity variables
function FixingCap(M, ref_EV_cap)
    for h in H
        #if Household_types[Household_types[!, "Type"] .== h, "PV"][1] == 0
            #for s in S
            #    fix(M[:C_PV][s], 0; force=true)
        #    end
    #    else
        #    for s in S
    #        fix(M[:C_PV][s], ref_PV_cap; force=true)
        #    end
    #    end
        #if Household_types[Household_types[!, "Type"] .== h, "BT"][1] == 0
        #    for s in S
            #    fix(M[:C_BT][s], 0; force=true)
        #    end
        #else
        #    for s in S
    #            fix(M[:C_BT][s], ref_BT_cap; force=true)
        #    end
    #    end
        if Household_types[Household_types[!, "Type"] .== h, "EV"][1] == 0
           for s in S
                fix(M[:C_EV][h,s], 0; force=true)
            end
        else
            for s in S
               fix(M[:C_EV][h,s], ref_EV_cap; force=true)
            end
        end
    end
end

# -------------------------------------------------------------------------------------------------------------------
#                                   CALCULATING THE PARAMETERS TO USE IN THE MODEL
# -------------------------------------------------------------------------------------------------------------------
function CalculatingParameters(max_charge)

    # Initializing empty 4-dimentional array based on sets t and S
    # Demand = Array{Float64}(undef, length(T), length(Y), length(H), length(S))
    Demand = Dict{Tuple{Int64,String,Int64},Float64}()
    # Assigning the imported data frame with specific demand profile to demand dictionary
    for s = 1
        for h in H
            for t in T
                Demand[(t,h,s)] = Demand_profiles[t,h]
            end
        end
    end
    for s= 2
        for h in H
            for t in T
                Demand[(t,h,s)] = Demand_profiles2[t,h]
            end
        end
    end
    # Initializing empty 4-dimentional array based on sets t and S
    # Demand = Array{Float64}(undef, length(T), length(Y), length(H), length(S))
    H_Demand = Dict{Tuple{Int64,String,Int64},Float64}()
    # Assigning the imported data frame with specific demand profile to demand dictionary
    for s = 1
            for h in H
                for t in T
                    H_Demand[(t,h,s)] = Heat_profiles[t,h]
                end
            end
    end
    for s = 2
            for h in H
                for t in T
                    H_Demand[(t,h,s)] = Heat_profiles2[t,h]
                end
            end
    end
    # Initializing empty 3-dimentional array based on sets t and S
    EV_avail = Dict{Tuple{Int64,String,Int64},Float64}()
    # Assigning the availability data imported from csv file
    for s in S
            for h in H
                for t in T
                    EV_avail[(t,h,s)] = EV_DF[t,h]
                end
            end
    end

    # Initializing empty 3-dimentional array based on sets t and S
    EV_demand = Dict{Tuple{Int64,String,Int64},Float64}()
    # Assigning data to the array based on the availability array as well as the second input to the function.
    for s in S
            for h in H
                for i=1:(size(EV_DF,1)-1)
                    if Household_types[Household_types[!, "Type"] .== h, "EV"][1] != 0
                        if EV_DF[i,h]==1 && EV_DF[i+1,h]==0
                            if  EV_cons[i,h]>=max_charge/1.2
                                EV_demand[(i,h,s)]= max_charge/1.2
                            else
                                EV_demand[(i,h,s)]=EV_cons[i,h]
                            end
                        else
                            EV_demand[(i,h,s)] = 0
                        end
                    else
                        EV_demand[(i,h,s)] = 0
                    end
                end
                EV_demand[(size(EV_DF,1),h,s)] = 0
            end
    end

    # Initializing empty 3-dimentional array based on sets t and S
    EV_SOC_goal = Dict{Tuple{Int64,String,Int64},Float64}()
    #EV_SOC_goal[:,:,:] .= 0
    # Assigning data to the array based on the availability array as well as the first input to the function.
    #the 1.2 makes sure that the demand is always set to be below the state of charge, to ensure that there always is a buffer.
    for s in S
            for h in H
                for i=1:(size(EV_DF,1)-1)
                    if Household_types[Household_types[!, "Type"] .== h, "EV"][1] != 0
                        if EV_DF[i,h]==1 && EV_DF[i+1,h]==0

                            if  EV_cons[i,h]>=max_charge/1.2
                                EV_SOC_goal[(i,h,s)] = max_charge
                            else
                                EV_SOC_goal[(i,h,s)]=EV_cons[i,h]*1.2
                            end
                        else
                            EV_SOC_goal[(i,h,s)] = 0
                        end
                    else
                        EV_SOC_goal[(i,h,s)] = 0
                    end
                end
                EV_SOC_goal[(size(EV_DF,1),h,s)] = 0
            end
    end
    return Demand, H_Demand, EV_avail, EV_demand, EV_SOC_goal

end


# -------------------------------------------------------------------------------------------------------------------
#                                     CONSTRAINTS
# -------------------------------------------------------------------------------------------------------------------
function  DefineConstraints(M, scheme)
    if scheme == "new"

    @objective(M, Min, sum((PV_par["Capital_cost"]/PV_par["Lifetime"])*M[:C_PV][s]
    + (Battery_par["Capital_cost"]/Battery_par["Lifetime"])*M[:C_BT][s]
    + (P2H_par["Capital_cost"]/P2H_par["Lifetime"])*M[:C_P2H][s]
    + (Boiler_par["Capital_cost"]/Boiler_par["Lifetime"] )*M[:C_EB][s]
    + (Heat_Storage_par["Capital_cost"]/Heat_Storage_par["Lifetime"])*M[:C_HS][s]
    + sum(PV_par["OM_cost"] * M[:p_PV][t,s]
    +  Battery_par["OM_cost"] * (M[:b_dh][t,s] + M[:b_ch][t,s])
    + P2H_par["OM_cost"] * M[:p2h_p][t,s]
    + Boiler_par["OM_cost"] * M[:eb_p][t,s]
    + Heat_Storage_par["OM_cost"] * (M[:hs_dh][t,s] + M[:hs_ch][t,s]) for t in T )
     #+ Battery_par["OP_cost"] * sum(M[:b_dh][t,s] + M[:b_ch][t,s] for t in T )
    + sum(M[:g_im][t,s] * (El_price[t,"Import"] + Network_tariffs["Var_dist"] + Network_tariffs["TSO"] + Network_tariffs["Tax"]+ EL_CO2[t,"EL average"]*0.020  )
    - M[:g_ex][t,s] * El_price[t,"Export"] for t in T )
    #+ sum(M[:ec_im][t,c,s] * El_price[t,"Import"]
    #- M[:ec_ex][t,s,c] * El_price[t,"Export"] for t in T for c in S )
    + sum(Network_tariffs["Fixed_dist"] for h in H) for s in S)
    + Ec_par["Capital_cost"]/40)

    elseif scheme == "base"
    @objective(M, Min, sum((PV_par["Capital_cost"]/PV_par["Lifetime"])*M[:C_PV][s]
    + (Battery_par["Capital_cost"]/Battery_par["Lifetime"])*M[:C_BT][s]
    + (P2H_par["Capital_cost"]/P2H_par["Lifetime"])*M[:C_P2H][s]
    + (Boiler_par["Capital_cost"]/Boiler_par["Lifetime"] )*M[:C_EB][s]
    + (Heat_Storage_par["Capital_cost"]/Heat_Storage_par["Lifetime"])*M[:C_HS][s]
    + sum(PV_par["OM_cost"] * M[:p_PV][t,s]
    +  Battery_par["OM_cost"] * (M[:b_dh][t,s] + M[:b_ch][t,s])
    + P2H_par["OM_cost"] * M[:p2h_p][t,s]
    + Boiler_par["OM_cost"] * M[:eb_p][t,s]
    + Heat_Storage_par["OM_cost"] * (M[:hs_dh][t,s] + M[:hs_ch][t,s]) for t in T )
     #+ Battery_par["OP_cost"] * sum(M[:b_dh][t,s] + M[:b_ch][t,s] for t in T )
    + sum(M[:g_im][t,s] * ((El_price[t,"Import"]) + Network_tariffs["Var_dist"] + Network_tariffs["TSO"] + Network_tariffs["Tax"]+ EL_CO2[t,"EL average"]*0.020  )
    - M[:g_ex][t,s] * (El_price[t,"Export"]) for t in T )
    #+ sum(M[:ec_im][t,c,s] * El_price[t,"Import"] + Network_tariffs["Var_dist"] + Network_tariffs["PSO"] + Network_tariffs["Tax"]
    #- M[:ec_ex][t,s,c] * El_price[t,"Export"] for t in T for c in S )
    + sum(Network_tariffs["Fixed_dist"] for h in H) for s in S))




    end

    # Balancing constraint taking into account only load flows
    @constraint(M, Balance[t in T, s in S], M[:g_im_load][t,s] + M[:b_dh_load][t,s] + sum(M[:ev_dh_load][t,h,s] for h in H ) + M[:p_PV_load][t,s] + M[:ec_im_load][t,s] - sum(Demand[(t,h,s)] for h in H ) == 0)

    # Heat balancing constraint taking into account only load flows
    @constraint(M, H_Balance[t in T, s in S], M[:p2h_load][t,s] + M[:eb_load][t,s] + M[:hs_dh_load][t,s] - sum(H_Demand[(t,h,s)] for h in H ) == 0)

    ################ PV

    # Definition of the PV array production
    @constraint(M, PV_prod_def[t in T,  s in S], M[:C_PV][s]*PV_CF[t,s] == M[:p_PV][t,s])
    # Balance of the PV energy
    @constraint(M, PV_prod_bal[t in T, s in S], M[:p_PV][t,s] == M[:p_PV_ev][t,s] + M[:p_PV_bat][t,s] + M[:p_PV_ex][t,s] + M[:p_PV_load][t,s] + M[:p_PV_p2h][t,s] + M[:p_PV_eb][t,s] + M[:p_PV_ec][t,s])

    @constraint(M, PV_prod_sens[t in T,  s in S ; s == 2], M[:C_PV][s] <= 40)

    ############# Battery

    # Limit on the maximum charge state of charge of the battery
    @constraint(M, SOC_lim_up[t in T,  s in S], M[:b_st][t,s] <= Battery_par["Max_charge"]*M[:C_BT][s])
    # Limit on the minimum battery state of charge
    @constraint(M, SOC_lim_down[t in T,  s in S], M[:b_st][t,s] >= Battery_par["Min_charge"]*M[:C_BT][s])

    # State Of Charge regular balance when the hours set is not 1
    @constraint(M, SOC[t in T,  s in S; t>1], M[:b_st][t,s] == M[:b_st][t-1,s] - M[:b_dh][t,s]/Battery_par["Discharging_eff"] + M[:b_ch][t,s]*Battery_par["Charging_eff"])
    # State Of Charge balance for the first hour and NOT first year
    #@constraint(M, SOC_LastT[t in T,  s in S; t == 1 && y!=1], M[:b_st][t,s] == M[:b_st][last(T),y-1,s] - M[:b_dh][t,s]/Battery_par["Discharging_eff"] + M[:b_ch][t,s]*Battery_par["Charging_eff"])
    # State Of Charge balance for the first hour and first year
    @constraint(M, SOC_First[t in T,  s in S; t==1], M[:b_st][t,s] == M[:C_BT][s] - M[:b_dh][t,s]/Battery_par["Discharging_eff"] + M[:b_ch][t,s]*Battery_par["Charging_eff"] )

    # Limit on the maximum hourly charging
    @constraint(M, Charge_limit[t in T,  s in S], M[:b_ch][t,s] <= Battery_par["Charging_lim"]*M[:C_BT][s])
    # Limit on the maximum hourly discharging
    @constraint(M, Discharge_limit[t in T,  s in S], M[:b_dh][t,s] <= Battery_par["Discharging_lim"]*M[:C_BT][s])

    # Balance of the charging energy
    @constraint(M, bat_ch_def[t in T,  s in S], M[:b_ch][t,s] == M[:p_PV_bat][t,s] + M[:g_im_bat][t,s] + M[:ec_im_bat][t,s])

    # Balance of the discharging energy
    @constraint(M, bat_dh_def[t in T,  s in S], M[:b_dh][t,s] == M[:b_dh_ex][t,s] + M[:b_dh_load][t,s] + M[:b_dh_ev][t,s] + M[:b_dh_p2h][t,s] +  M[:b_dh_eb][t,s] + M[:b_dh_ec][t,s])

    ############# EV

    # SOC (EV) regular balance when the hours set is not 1
    @constraint(M, SOC_EV[t in T, h in H, s in S; t>1], M[:ev_st][t,h,s] == M[:ev_st][t-1,h,s] - M[:ev_dh][t,h,s]/EV_par["Discharging_eff"] + M[:ev_ch][t,h,s]*EV_par["Charging_eff"] - EV_demand[(t,h,s)])
    # SOC (EV) balance for the first hour
    @constraint(M, SOC_First_EV[t in T, h in H, s in S; t==1], M[:ev_st][t,h,s] == M[:C_EV][h,s] - M[:ev_dh][t,h,s]/EV_par["Discharging_eff"] + M[:ev_ch][t,h,s]*EV_par["Charging_eff"]  - EV_demand[(t,h,s)])

    # Limit on the maximum charge state of charge of the battery
    @constraint(M, SOC_lim_up_EV[t in T, h in H, s in S], M[:ev_st][t,h,s] <= EV_par["Max_charge"]*M[:C_EV][h,s])
    # Limit on the minimum battery state of charge
    @constraint(M, SOC_lim_down_EV[t in T, h in H, s in S], M[:ev_st][t,h,s] >= EV_par["Min_charge"]*M[:C_EV][h,s])

    # Definition of the EV charging limit taking into accout availability
    @constraint(M, EV_charge_lim[t in T, h in H, s in S], M[:ev_ch][t,h,s] <=  EV_avail[(t,h,s)] * EV_par["Charging_lim"] * M[:C_EV][h,s])
    # Definition of the EV discharging limit taking into accout availability
    @constraint(M, EV_discharge_lim[t in T, h in H, s in S], M[:ev_dh][t,h,s] <=  EV_avail[(t,h,s)] * EV_par["Discharging_lim"] * M[:C_EV][h,s])
    # Definition of the EV discharging limit taking into accout availability
    @constraint(M, EV_SOC_goal_cons[t in T, h in H, s in S], M[:ev_st][t,h,s] >=  EV_SOC_goal[(t,h,s)])

    # Balance of the charging energy in EV
    @constraint(M, ev_ch_def[t in T, s in S], sum(M[:ev_ch][t,h,s] for h in H) == M[:p_PV_ev][t,s] + M[:g_im_ev][t,s] + M[:b_dh_ev][t,s] + M[:ec_im_ev][t,s] )
    # Balance of the discharging energy in EV
    @constraint(M, ev_dh_def[t in T, h in H, s in S], M[:ev_dh][t,h,s] == M[:ev_dh_ex][t,h,s] + M[:ev_dh_load][t,h,s] + M[:ev_dh_ec][t,h,s]+M[:ev_dh_p2h][t,h,s]+M[:ev_dh_eb][t,h,s] )
    # dh ec 0
    #@constraint(M, ev_dhno[t in T, h in H, s in S], M[:ev_dh_ec][t,h,s]== 0)

    ################### Grid

    # Limit on the amount of hourly exported electricity
    @constraint(M, grid_ex_lim[t in T,  s in S], M[:g_ex][t,s] <= Grid_par["Ex_lim"])
    # Limit on the amount of hourly imported electricity
    @constraint(M, grid_im_lim[t in T,  s in S], M[:g_im][t,s] <= Grid_par["Im_lim"])


    # Balance of the imported energy
    @constraint(M, grid_im_def[t in T ,s in S], M[:g_im][t,s] == M[:g_im_load][t,s] + M[:g_im_bat][t,s] + M[:g_im_ev][t,s] + M[:g_im_p2h][t,s] + M[:g_im_eb][t,s])
    # Balance of the exported energy
    @constraint(M, grid_ex_def[t in T,  s in S], M[:g_ex][t,s] == M[:p_PV_ex][t,s] + M[:b_dh_ex][t,s] + sum(M[:ev_dh_ex][t,h,s] for h in H))


    #################### EC

    if scheme == "new"
    @constraint(M, ec_bal[t in T ,s in S, c in S], M[:ec_im][t,c,s] == M[:ec_ex][t,s,c])

    @constraint(M, ec_bal2[t in T ,s in S, c in S; c==s], M[:ec_im][t,c,s] + M[:ec_ex][t,c,s]  == 0)

    @constraint(M, ec_rest[t in T ,s in S, c in S], sum(M[:ec_im][t,c,s] for c in S)  <= M[:ec_im_var][t,s]*5000)

    @constraint(M, ec_rest1[t in T ,s in S, c in S], sum(M[:ec_ex][t,c,s] for c in S) <= (1-M[:ec_im_var][t,s])*5000)

    elseif scheme == "base"
    @constraint(M, ec_bal[t in T ,s in S, c in S], M[:ec_im][t,c,s] == 0)

    @constraint(M, ec_bal2[t in T ,s in S, c in S; c==s], M[:ec_ex][t,c,s]  == 0)
    end

    # Limit on the amount of hourly exported electricity
    @constraint(M, Ec_ex_lim[t in T,  s in S], sum(M[:ec_ex][t,c,s] for c in S) <=  Ec_par["Ex_lim"])
    # Limit on the amount of hourly imported electricity
    @constraint(M, Ec_im_lim[t in T,  s in S], sum(M[:ec_im][t,c,s] for c in S) <=  Ec_par["Im_lim"])


    # Balance of the imported energy
    @constraint(M, Ec_im_def[t in T , s in S], sum(M[:ec_im][t,c,s] for c in S) == M[:ec_im_load][t,s] + M[:ec_im_bat][t,s] + M[:ec_im_ev][t,s] + M[:ec_im_p2h][t,s] + M[:ec_im_eb][t,s])
    # Balance of the exported energy
    @constraint(M, Ec_ex_def[t in T,  s in S], sum(M[:ec_ex][t,c,s] for c in S) == M[:p_PV_ec][t,s] + M[:b_dh_ec][t,s] + sum(M[:ev_dh_ec][t,h,s] for h in H))



    ################### P2H
    @constraint(M, p2h_cap[t in T,  s in S], M[:p2h_p][t,s] <= M[:C_P2H][s])

    @constraint(M, p2h_im_def[t in T,  s in S], M[:p2h_im][t,s] == M[:p_PV_p2h][t,s] + M[:g_im_p2h][t,s] + M[:b_dh_p2h][t,s]  +  sum(M[:ev_dh_p2h][t,h,s] for h in H) + M[:ec_im_p2h][t,s])

    @constraint(M, p2h_ex_def[t in T,  s in S], M[:p2h_p][t,s]/P2H_par["COP"] == M[:p2h_im][t,s])

    @constraint(M, p2h_prod_bal[t in T,  s in S], M[:p2h_p][t,s] == M[:p2h_load][t,s] + M[:p2h_hs][t,s])

    ###################### EB
    @constraint(M, eb_cap[t in T,  s in S], M[:eb_p][t,s] <= M[:C_EB][s])

    @constraint(M, eb_im_def[t in T,  s in S], M[:eb_im][t,s] == M[:p_PV_eb][t,s] + M[:g_im_eb][t,s] + M[:b_dh_eb][t,s] + sum(M[:ev_dh_eb][t,h,s] for h in H) + M[:ec_im_eb][t,s])

    @constraint(M, eb_ex_def[t in T,  s in S], M[:eb_p][t,s] == M[:eb_im][t,s])

    @constraint(M, eb_prod_bal[t in T,  s in S], M[:eb_p][t,s] == M[:eb_load][t,s] + M[:eb_hs][t,s])

    ###################### HS
    # SOC (hs) regular balance when the hours set is not 1
    @constraint(M, SOC_HS[t in T,  s in S; t>1], M[:hs_st][t,s] == M[:hs_st][t-1,s] - M[:hs_dh][t,s]/Heat_Storage_par["Discharging_eff"] + M[:hs_ch][t,s]*Heat_Storage_par["Charging_eff"])
    # SOC (hs) balance for the first hour and NOT first year
    #@constraint(M, SOC_LastT_HS[t in T,  s in S; t == 1 && y!=1], M[:hs_st][t,s] == M[:hs_st][last(T),y-1,s] - M[:hs_dh][t,s]/Heat_Storage_par["Discharging_eff"] + M[:hs_ch][t,s]*Heat_Storage_par["Charging_eff"])
    # SOC (hs) balance for the first hour and first year
    @constraint(M, SOC_First_HS[t in T,  s in S; t==1], M[:hs_st][t,s] == M[:C_HS][s] - M[:hs_dh][t,s]/Heat_Storage_par["Discharging_eff"] + M[:hs_ch][t,s]*Heat_Storage_par["Charging_eff"])

    # Limit on the maximum charge state of charge of the battery
    @constraint(M, SOC_lim_up_HS[t in T,  s in S], M[:hs_st][t,s] <= Heat_Storage_par["Max_charge"]*M[:C_HS][s])
    # Limit on the minimum battery state of charge
    @constraint(M, SOC_lim_down_HS[t in T,  s in S], M[:hs_st][t,s] >= Heat_Storage_par["Min_charge"]*M[:C_HS][s])

    # Balance of the charging energy in hs
    @constraint(M, HS_ch_def[t in T,  s in S], M[:hs_ch][t,s] == M[:eb_hs][t,s] + M[:p2h_hs][t,s])
    # Balance of the discharging energy in hs
    @constraint(M, HS_dh_def[t in T,  s in S], M[:hs_dh][t,s] == M[:hs_dh_load][t,s])

    # Limit on the maximum hourly charging
    @constraint(M, HS_Charge_limit[t in T,  s in S], M[:hs_ch][t,s] <= Heat_Storage_par["Charging_lim"]*M[:C_HS][s])
    # Limit on the maximum hourly discharging
    @constraint(M, HS_Discharge_limit[t in T,  s in S], M[:hs_dh][t,s] <= Heat_Storage_par["Discharging_lim"]*M[:C_HS][s])




    return M
end



#***********************************************************************
#                Trying to get the vvariables structured
#***********************************************************************
function var_info(vars::VariableRef)
    split_conv = [":","]","[",","]
    x_str = name(vars)
    if occursin(":",x_str)
        x_str = replace(x_str, " " => "") #Deletes all spaces
        x_name,x_index = split(x_str,split_conv[1]) #splits raw variable name+ sets and index
        x_name = replace(x_name, split_conv[2] => "")
        x_name,s_set = split(x_name,split_conv[3])#splits raw variable name and sets
        x_set = split(s_set,split_conv[4])

        x_index = replace(x_index, split_conv[2] => "")
        x_index = replace(x_index, split_conv[3] => "")
        x_index = split(x_index,split_conv[4])
        return (x_name,x_set,x_index)
    else
        println("Var base_name not properly defined. Special Syntax required in form     var[s=set]:  ")
    end

end
function create_columns(x)
    col_ind=[String(var_info(x)[2][col]) for col in 1:size(var_info(x)[2])[1]]
    cols = append!(["Value"],col_ind)
    return cols
end

function create_index(x)
    col_ind=[String(var_info(x)[3][ind]) for ind in 1:size(var_info(x)[3])[1]]
    index = append!([string(value(x))],col_ind)
    return index
end

function create_sol_matrix(varss,model)
    nested_sol_array=[create_index(xx) for xx in all_variables(model) if varss[1]==var_info(xx)[1]]
    sol_array=hcat(nested_sol_array...)
    return sol_array
end

function create_var_dict(model)
    Variable_dict=Dict(vars[1]
            =>DataFrame(Dict(vars[2][1][cols]
                =>create_sol_matrix(vars,model)[cols,:] for cols in 1:size(vars[2][1])[1]))
                    for vars in unique([[String(var_info(x)[1]),[create_columns(x)]] for x in all_variables(model)]))
    return Variable_dict
end

#***********************************************************************
#                Exporting Results
#***********************************************************************
function results()
    XLSX.writetable("Resultsbaseheatlow.xlsx", overwrite=true, Cap_PV=(collect(DataFrames.eachcol(var_dict["C_PV"])),DataFrames.names(var_dict["C_PV"])),
                                                    Cap_BT=(collect(DataFrames.eachcol(var_dict["C_BT"])),DataFrames.names(var_dict["C_BT"])),
                                                    Cap_EV=(collect(DataFrames.eachcol(var_dict["C_EV"])),DataFrames.names(var_dict["C_EV"])),
                                                    Cap_P2H=(collect(DataFrames.eachcol(var_dict["C_P2H"])),DataFrames.names(var_dict["C_P2H"])),
                                                    Cap_EB=(collect(DataFrames.eachcol(var_dict["C_EB"])),DataFrames.names(var_dict["C_EB"])),
                                                    Cap_HS=(collect(DataFrames.eachcol(var_dict["C_HS"])),DataFrames.names(var_dict["C_HS"])),

                                                    Import_grid=(collect(DataFrames.eachcol(var_dict["g_im"])),DataFrames.names(var_dict["g_im"])),
                                                    Export_grid=(collect(DataFrames.eachcol(var_dict["g_ex"])),DataFrames.names(var_dict["g_ex"])),
                                                    Import_to_EV=(collect(DataFrames.eachcol(var_dict["g_im_ev"])),DataFrames.names(var_dict["g_im_ev"])),
                                                    #Import_to_BT=(collect(DataFrames.eachcol(var_dict["g_im_bat"])),DataFrames.names(var_dict["g_im_bat"])),
                                                #    Import_to_load=(collect(DataFrames.eachcol(var_dict["g_im_load"])),DataFrames.names(var_dict["g_im_load"])),
                                                    #Import_to_p2h=(collect(DataFrames.eachcol(var_dict["g_im_p2h"])),DataFrames.names(var_dict["g_im_p2h"])),
                                                    #Import_to_eb=(collect(DataFrames.eachcol(var_dict["g_im_eb"])),DataFrames.names(var_dict["g_im_eb"])),

                                                    PV_prod=(collect(DataFrames.eachcol(var_dict["p_PV"])),DataFrames.names(var_dict["p_PV"])),
                                                    PV_to_EV=(collect(DataFrames.eachcol(var_dict["p_PV_ev"])),DataFrames.names(var_dict["p_PV_ev"])),
                                                    #PV_to_BT=(collect(DataFrames.eachcol(var_dict["p_PV_bat"])),DataFrames.names(var_dict["p_PV_bat"])),
                                                    #PV_to_load=(collect(DataFrames.eachcol(var_dict["p_PV_load"])),DataFrames.names(var_dict["p_PV_load"])),
                                                    #PV_to_p2h=(collect(DataFrames.eachcol(var_dict["p_PV_p2h"])),DataFrames.names(var_dict["p_PV_p2h"])),
                                                    #PV_to_export=(collect(DataFrames.eachcol(var_dict["p_PV_ex"])),DataFrames.names(var_dict["p_PV_ex"])),
                                                    #PV_to_Eb=(collect(DataFrames.eachcol(var_dict["p_PV_eb"])),DataFrames.names(var_dict["p_PV_eb"])),
                                                    PV_to_Ec=(collect(DataFrames.eachcol(var_dict["p_PV_ec"])),DataFrames.names(var_dict["p_PV_ec"])),

                                                    BT_charge=(collect(DataFrames.eachcol(var_dict["b_ch"])),DataFrames.names(var_dict["b_ch"])),
                                                    BT_discharge=(collect(DataFrames.eachcol(var_dict["b_dh"])),DataFrames.names(var_dict["b_dh"])),
                                                    #BT_state=(collect(DataFrames.eachcol(var_dict["b_st"])),DataFrames.names(var_dict["b_st"])),
                                                    #BT_to_export=(collect(DataFrames.eachcol(var_dict["b_dh_ex"])),DataFrames.names(var_dict["b_dh_ex"])),
                                                #    BT_to_EV=(collect(DataFrames.eachcol(var_dict["b_dh_ev"])),DataFrames.names(var_dict["b_dh_ev"])),
                                                    BT_to_EC=(collect(DataFrames.eachcol(var_dict["b_dh_ec"])),DataFrames.names(var_dict["b_dh_ec"])),
                                                    #BT_to_load=(collect(DataFrames.eachcol(var_dict["b_dh_load"])),DataFrames.names(var_dict["b_dh_load"])),
                                                    #BT_to_p2h=(collect(DataFrames.eachcol(var_dict["b_dh_p2h"])),DataFrames.names(var_dict["b_dh_p2h"])),
                                                    #BT_to_eb=(collect(DataFrames.eachcol(var_dict["b_dh_eb"])),DataFrames.names(var_dict["b_dh_eb"])),

                                                    EV_charge=(collect(DataFrames.eachcol(var_dict["ev_ch"])),DataFrames.names(var_dict["ev_ch"])),
                                                    EV_discharge=(collect(DataFrames.eachcol(var_dict["ev_dh"])),DataFrames.names(var_dict["ev_dh"])),
                                                    EV_state=(collect(DataFrames.eachcol(var_dict["ev_st"])),DataFrames.names(var_dict["ev_st"])),
                                                    EV_to_export=(collect(DataFrames.eachcol(var_dict["ev_dh_ex"])),DataFrames.names(var_dict["ev_dh_ex"])),
                                                    EV_to_load=(collect(DataFrames.eachcol(var_dict["ev_dh_load"])),DataFrames.names(var_dict["ev_dh_load"])),
                                                    EV_to_ec=(collect(DataFrames.eachcol(var_dict["ev_dh_ec"])),DataFrames.names(var_dict["ev_dh_ec"])),

                                                    El_for_p2h=(collect(DataFrames.eachcol(var_dict["p2h_im"])),DataFrames.names(var_dict["p2h_im"])),
                                                    P2h_produced=(collect(DataFrames.eachcol(var_dict["p2h_p"])),DataFrames.names(var_dict["p2h_p"])),
                                                    #P2h_load=(collect(DataFrames.eachcol(var_dict["p2h_load"])),DataFrames.names(var_dict["p2h_load"])),
                                                    #P2h_storage=(collect(DataFrames.eachcol(var_dict["p2h_hs"])),DataFrames.names(var_dict["p2h_hs"])),

                                                    El_for_eb=(collect(DataFrames.eachcol(var_dict["eb_im"])),DataFrames.names(var_dict["eb_im"])),
                                                    eb_produced=(collect(DataFrames.eachcol(var_dict["eb_p"])),DataFrames.names(var_dict["eb_p"])),
                                                    #eb_load=(collect(DataFrames.eachcol(var_dict["eb_load"])),DataFrames.names(var_dict["eb_load"])),
                                                    #eb_storage=(collect(DataFrames.eachcol(var_dict["eb_hs"])),DataFrames.names(var_dict["eb_hs"])),

                                                    hs_charge=(collect(DataFrames.eachcol(var_dict["hs_ch"])),DataFrames.names(var_dict["hs_ch"])),
                                                    hs_discharge=(collect(DataFrames.eachcol(var_dict["hs_dh"])),DataFrames.names(var_dict["hs_dh"])),
                                                    #hs_state=(collect(DataFrames.eachcol(var_dict["hs_st"])),DataFrames.names(var_dict["hs_st"])),
                                                    #hs_to_load=(collect(DataFrames.eachcol(var_dict["hs_dh_load"])),DataFrames.names(var_dict["hs_dh_load"])),

                                                    EC_export=(collect(DataFrames.eachcol(var_dict["ec_ex"])),DataFrames.names(var_dict["ec_ex"])),
                                                    EC_import=(collect(DataFrames.eachcol(var_dict["ec_im"])),DataFrames.names(var_dict["ec_im"])),
                                                    EC_load=(collect(DataFrames.eachcol(var_dict["ec_im_load"])),DataFrames.names(var_dict["ec_im_load"])),
                                                    EC_bat=(collect(DataFrames.eachcol(var_dict["ec_im_bat"])),DataFrames.names(var_dict["ec_im_bat"])),
                                                    EC_ev=(collect(DataFrames.eachcol(var_dict["ec_im_ev"])),DataFrames.names(var_dict["ec_im_ev"])),
                                                    EC_p2h=(collect(DataFrames.eachcol(var_dict["ec_im_p2h"])),DataFrames.names(var_dict["ec_im_p2h"])),
                                                    EC_eb=(collect(DataFrames.eachcol(var_dict["ec_im_eb"])),DataFrames.names(var_dict["ec_im_eb"]))
                                                    )
end

#**************************************************************************
#
#************************************************************************
ModelDataImport()
# Initializng the model and variables
M = InitializeModel("new")
# Calculating the parameters
Demand, H_Demand, EV_avail, EV_demand, EV_SOC_goal =CalculatingParameters(30)
# Fixing the capacities of PV, Battery and EV
FixingCap(M, 30)
# Fixing the capacities
M = DefineConstraints(M, "new")
# Optimizing!
optimize!(M)
println("Solution:", value.(sum(M[:g_im][t,s] * EL_CO2[t,"EL average"] for s in 1 for t in T)))
println("Solution:", value.(sum(M[:g_im][t,s] * EL_CO2[t,"EL average"] for s in 2 for t in T)))
var_dict = create_var_dict(M)
results()

CSV.write("EV_demand.csv",EV_demand)


XLSX.writetable("EV_demand.xlsx", overwrite=true,(collect(DataFrames.eachcol(EV_demand)),DataFrames.names(EV_demand)))


println("Solution:", value.(sum(M[:g_im][t,s] * EL_CO2[t,"EL average"] for s in 1 for t in T)))

println("Solution:", value.((M[:C_HS][s] for s in S)))

println("Solution:", value.(sum(EV_demand[(t,h,s)] for s in 2 for h in H for t in T)))

println("Termination status: $(termination_status(M))")

#************************************************************************
#                   Plotting
#************************************************************************


global Cap_PV = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "Cap_PV")...)
global Cap_BT = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "Cap_BT")...)
global Cap_EV = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "Cap_EV")...)
global Cap_P2H = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "Cap_EV")...)
global Import_grid = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "Import_grid")...)
global Export_grid = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "Export_grid")...)
global Import_to_EV = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "Import_to_ev")...)
global Import_to_BT = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "Import_to_bt")...)
global Import_to_load = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "Import_to_load")...)
global Import_to_heat = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "Import_to_heat")...)
global PV_prod = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "PV_prod")...)
global PV_to_EV = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "PV_to_ev")...)
global PV_to_BT = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "PV_to_bt")...)
global PV_to_load = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "PV_to_load")...)
global PV_to_p2h = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "PV_to_p2h")...)
global PV_to_export = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "PV_to_export")...)
global BT_charge = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "BT_charge")...)
global BT_discharge = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "BT_discharge")...)
global BT_state = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "BT_state")...)
global BT_to_export = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "BT_to_export")...)
global BT_to_EV = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "BT_to_ev")...)
global BT_to_load = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "BT_to_load")...)
global BT_to_p2h = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "BT_to_p2h")...)
global EV_charge = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "EV_charge")...)
global EV_discharge = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "EV_discharge")...)
global EV_state = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "EV_state")...)
global EV_to_export = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "EV_to_export")...)
global EV_to_load = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "EV_to_load")...)
global El_for_p2h = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "El_for_p2h")...)
global Heat_produced = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "Heat_produced")...)
global Capacities = DataFrame(XLSX.readtable("Results/Struct_results1.xlsx", "Capacities")...)


Investments = [Capacities[!,"Cap_PV"]*Scalars["CRF"]*PV_par["Capital_cost"] Capacities[!,"Cap_BT"]*Scalars["CRF"]*Battery_par["Capital_cost"] Capacities[!,"Cap_EV"] Capacities[!,"Cap_P2H"]*Scalars["CRF"]*P2H_par["Capital_cost"] ]
Caps = [Capacities[!,"Cap_PV"] Capacities[!,"Cap_BT"] Capacities[!,"Cap_EV"] Capacities[!,"Cap_P2H"]]
Caps1 = [7.34742 0 270 82.7717]
Investments1 = [864.792, 0, 270, 11604.6]
ticklabel = ["Cap_PV", "Cap_BT", "Cap_EV", "Cap_P2H"]
Investment = Int.(Investments)

capplot = bar(["Cap_PV" "Cap_BT" "Cap_EV" "Cap_P2H"], [Caps1],bar_position = :dodge, ylabel= "KW /KWh installed", label = (["Cap_PV" "Cap_BT" "Cap_EV" "Cap_P2H"]) )
savefig(capplot,"capplot.png")

invplot = bar(["Inv_PV" "Inv_BT" "Inv_EV" "Inv_P2H"], [Investments1],bar_position = :dodge, ylabel= "Eur", label = (["Inv_PV" "Inv_BT" "Inv_EV" "Inv_P2H"]) )
savefig(invplot,"invplot.png")
EmmisionCost =

EL_CO2[t,"EL average"]*0.020
