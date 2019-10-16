# workspace();

include("OControllerProx.jl")
using .OControllerProx;
#include("Process.jl")
#using .Process;

using MathProgBase;
using LinearAlgebra;
using MAT;
using DataFrames;
using CSV;
using JuMP;
using ProxSDP;
using Dates;

# get the calculated switches status

path = "/Applications/projects/cLEAN/pnnl_control";

file = matopen("$(path)/input/one-time-input/switch.mat");
switches_status = read(file,"switches_status");

case=CSV.read("$(path)/input/one-time-input/case.csv");

#case = parse(Float64,CSV.read("$(path)/input/one-time-input/case.csv"));



#LoadInformation = CSV.File("$(path)/input/one-time-input/LoadInformation.xls"; normalizenames=true, types=[String, String, String, String, String, String, String, String]);
#Loads=CSV.read("$(path)/input/one-time-input/LoadInformation.xls"; normalizenames=true)

#using DataFrames
#using StringEncodings
#f=open("$(path)/input/one-time-input/LoadInformation.xls","r");
#s=StringDecoder(f,"LATIN1", "UTF-8");
#LoadInformation=CSV.read(s; normalizenames=true, comment="#");
#close(s)
#close(f)

#printf(Loads[2,2])


#(Bbus_processed,Gbus_processed) = Process.GBbus(Ybus, NodeIDs, NumberOfElements);
#Bbus_energizing=switches_status.*Bbus_processed;
#Gbus_energizing=switches_status.*Gbus_processed;

# get 5-minute actions

actions = CSV.read("$(path)/input/5-min-input/actions.csv"; comment="#");
optimal_injections=zeros(ComplexF64,8,1);
for i=1:8
optimal_injections[i]=parse(ComplexF64, actions[1,4*i+3]); # convert from string to complex number
end
number_of_energizing_segment=length(actions[:,1]);

# pnnl_actions csv file
global df = DataFrame(Dict("datatime" =>"", "generator.proposed_diesel_1.1"=> "", "generator.proposed_diesel_1.2" =>  "", "generator.proposed_diesel_1.3" =>  "", "generator.ouc_solar.1" =>  "","generator.ouc_solar.2" =>  "","generator.ouc_solar.3" =>  "","generator.gen_none_-60027406210.1" =>  "","generator.gen_none_-60027406210.2" => "", "generator.gen_none_-60027406210.3" =>  "","generator.proposed_solarpv_1.1" =>  "","generator.proposed_solarpv_1.2" =>  "","generator.proposed_solarpv_1.3" =>  "","generator.gen_none_6002693162.1" =>  "","generator.gen_none_6002693162.2" =>  "","generator.gen_none_6002693162.3" =>  "","storage.proposed_storage_1.1" =>  "","storage.proposed_storage_1.2" =>  "","storage.proposed_storage_1.3" =>  "","storage.proposed_fuelcell_1.1" =>  "","storage.proposed_fuelcell_1.2" =>  "","storage.proposed_fuelcell_1.3" =>  "","storage.proposed_fuelcell_2.1" =>  "","storage.proposed_fuelcell_2.2" =>  "","storage.proposed_fuelcell_2.3" =>  ""));
# save pnnl dynamic actions to the output folder
CSV.write("$(path)/output/pnnl_actions.csv", df);


nx=24; # 3 x number of (gen+storage)
mx=16;  # number of control actions
px=4;  # number of usable measurements
A=-Matrix{Float64}(I, 24,24);
B=[Matrix{Float64}(I, 16,16);zeros(8,16)];
C=zeros(4,24);C[1,3]=1;C[2,6]=1;C[3,9]=1;C[4,12]=1;
# finding controller and observer gains

Ax=A;
Bx=B;
Cx=C;

Delta1=1;
Delta2=1;
gamma = (Delta1+Delta2)/5;
                OptCtr = OControllerProx.solveControllerObserver(Ax,Bx,Cx,gamma, Delta1, Delta2,nx,mx,px);
                P_out=OptCtr[1];                       # Lyapunov matrix
                K_out=OptCtr[2];
                Q_out=OptCtr[3];                       # Lyapunov matrix
                L_out=OptCtr[4];                     # output control gain
##############################
# getting control signals to apply to control nodes in every 1 second for 5 minutes

current_observer=zeros(24,1);

global sensor_Mtr_382883
if case==1

for times = 1:number_of_energizing_segment
for i=1:300

    # remember to update actions when running for > 5 minute

# Step1: Get measurement data from measurement input folder
#sensor_Mtr_382883 =readtable("/Applications/projects/cLEAN/pnnl control/measurement input/sensor_Mtr_-382883.csv")
#df1 = CSV.File("/Applications/projects/cLEAN/pnnl control/measurement input/sensor_Mtr_-382883.csv", datarow=10);
#sensor_Mtr_382883 =CSV.read(df1)

current_measurement=zeros(4,1);
if isfile("$(path)/input/measurement-input/sensor_Mtr_-382883.csv")==false
  current_measurement[1]=0;
elseif isfile("$(path)/input/measurement-input/sensor_Mtr_1695832.csv")==false
  current_measurement[2]=0;
elseif isfile("$(path)/input/measurement-input/sensor_Mtr_118122225.csv")==false
  current_measurement[3]=0;
elseif isfile("$(path)/input/measurement-input/sensor_Mtr_119123112.csv")==false
  current_measurement[4]=0;
else

sensor_Mtr_382883 =CSV.read("$(path)/input/measurement-input/sensor_Mtr_-382883.csv"; comment="#");
sensor_Mtr_1695832 =CSV.read("$(path)/input/measurement-input/sensor_Mtr_1695832.csv"; comment="#");
sensor_Mtr_118122225 =CSV.read("$(path)/input/measurement-input/sensor_Mtr_118122225.csv"; comment="#");
sensor_Mtr_119123112 =CSV.read("$(path)/input/measurement-input/sensor_Mtr_119123112.csv"; comment="#");

# check if data at i-second instant is available; if not then wait until it is available
while (length(sensor_Mtr_382883[:,1])<(i+300*(times-1))) & (length(sensor_Mtr_1695832[:,1])<(i+300*(times-1))) & (length(sensor_Mtr_118122225[:,1])<(i+300*(times-1))) & (length(sensor_Mtr_119123112[:,1])<(i+300*(times-1)))
    sleep(0.001);
end

# calculate control input based on measurement at i-second instant

if typeof(sensor_Mtr_382883[2,2])==string
# convert from string to complex number measurement data
imcurrent_measurement=zeros(ComplexF64,4,1);
imcurrent_measurement[1]= (parse(ComplexF64, sensor_Mtr_382883[(i+300*(times-1)),2])+parse(ComplexF64, sensor_Mtr_382883[(i+300*(times-1)),3])+parse(ComplexF64, sensor_Mtr_382883[(i+300*(times-1)),4]))/3;
imcurrent_measurement[2]= (parse(ComplexF64, sensor_Mtr_1695832[(i+300*(times-1)),2])+parse(ComplexF64, sensor_Mtr_1695832[(i+300*(times-1)),3])+parse(ComplexF64, sensor_Mtr_1695832[(i+300*(times-1)),4]))/3;
imcurrent_measurement[3]= (parse(ComplexF64, sensor_Mtr_118122225[(i+300*(times-1)),2])+parse(ComplexF64, sensor_Mtr_118122225[(i+300*(times-1)),3])+parse(ComplexF64, sensor_Mtr_118122225[(i+300*(times-1)),4]))/3;
imcurrent_measurement[4]= (parse(ComplexF64, sensor_Mtr_119123112[(i+300*(times-1)),2])+parse(ComplexF64, sensor_Mtr_119123112[(i+300*(times-1)),3])+parse(ComplexF64, sensor_Mtr_119123112[(i+300*(times-1)),4]))/3;
  for j=1:4
    current_measurement[j]= abs(imcurrent_measurement[j]);
  end

elseif typeof(sensor_Mtr_382883[2,2])==Int64
current_measurement[1]= sqrt(sensor_Mtr_382883[(i+300*(times-1)),8]^2+sensor_Mtr_382883[(i+300*(times-1)),9]^2);
current_measurement[2]= sqrt(sensor_Mtr_1695832[(i+300*(times-1)),8]^2+sensor_Mtr_1695832[i(i+300*(times-1)),9]^2);
current_measurement[3]= sqrt(sensor_Mtr_118122225[(i+300*(times-1)),8]^2+sensor_Mtr_118122225[(i+300*(times-1)),9]^2);
current_measurement[4]= sqrt(sensor_Mtr_119123112[(i+300*(times-1)),8]^2+sensor_Mtr_119123112[(i+300*(times-1)),9]^2);
end
end

###############################
# find the next observer state
global current_observer
next_observer = current_observer + (A+B*K_out)*current_observer + L_out*(current_measurement-C*current_observer);

current_observer =next_observer;

################################
# Step 3: return PQ control signals bn
current_control = K_out*current_observer;
PQcurrent_control=zeros(ComplexF64,8,1);
for k=1:8
PQcurrent_control[k] = current_control[2*k-1] + current_control[2*k]*im;
end
################################
# Step 4: calculate the changes to apply to control nodes

dynamic_injections = optimal_injections + PQcurrent_control;

# write the file with the stringdata variable information
global df = DataFrame(Dict("datatime" => DateTime(2000,1,1,0,floor((i+300*(times-1))/60),(i+300*(times-1))-60*floor((i+300*(times-1))/60)), "generator.proposed_diesel_1.1" =>  [dynamic_injections[1]], "generator.proposed_diesel_1.2" =>  [dynamic_injections[1]], "generator.proposed_diesel_1.3" =>  [dynamic_injections[1]], "generator.ouc_solar.1" =>  [dynamic_injections[2]],"generator.ouc_solar.2" =>  [dynamic_injections[2]],"generator.ouc_solar.3" =>  [dynamic_injections[2]],"generator.gen_none_-60027406210.1" =>  [dynamic_injections[3]],"generator.gen_none_-60027406210.2" => [dynamic_injections[3]], "generator.gen_none_-60027406210.3" =>  [dynamic_injections[3]],"generator.proposed_solarpv_1.1" =>  [dynamic_injections[4]],"generator.proposed_solarpv_1.2" =>  [dynamic_injections[4]],"generator.proposed_solarpv_1.3" =>  [dynamic_injections[4]],"generator.gen_none_6002693162.1" =>  [dynamic_injections[5]],"generator.gen_none_6002693162.2" =>  [dynamic_injections[5]],"generator.gen_none_6002693162.3" =>  [dynamic_injections[5]],"storage.proposed_storage_1.1" =>  [dynamic_injections[6]],"storage.proposed_storage_1.2" =>  [dynamic_injections[6]],"storage.proposed_storage_1.3" =>  [dynamic_injections[6]],"storage.proposed_fuelcell_1.1" =>  [dynamic_injections[7]],"storage.proposed_fuelcell_1.2" =>  [dynamic_injections[7]],"storage.proposed_fuelcell_1.3" =>  [dynamic_injections[7]],"storage.proposed_fuelcell_2.1" =>  [dynamic_injections[8]],"storage.proposed_fuelcell_2.2" =>  [dynamic_injections[8]],"storage.proposed_fuelcell_2.3" =>  [dynamic_injections[8]]));
# save pnnl dynamic actions to the output folder
CSV.write("$(path)/output/pnnl_actions.csv", df,append=true);

end
end

elseif case==2 #incorrect measurement data
  for times = 1:number_of_energizing_segment
  for i=1:300

      # remember to update actions when running for > 5 minute

  # Step1: Get measurement data from measurement input folder
  #sensor_Mtr_382883 =readtable("/Applications/projects/cLEAN/pnnl control/measurement input/sensor_Mtr_-382883.csv")
  #df1 = CSV.File("/Applications/projects/cLEAN/pnnl control/measurement input/sensor_Mtr_-382883.csv", datarow=10);
  #sensor_Mtr_382883 =CSV.read(df1)
  current_measurement=zeros(4,1);
  if isfile("$(path)/input/measurement-input/sensor_Mtr_-382883.csv")==false
    current_measurement[1]=0;
  elseif isfile("$(path)/input/measurement-input/sensor_Mtr_1695832.csv")==false
    current_measurement[2]=0;
  elseif isfile("$(path)/input/measurement-input/sensor_Mtr_118122225.csv")==false
    current_measurement[3]=0;
  elseif isfile("$(path)/input/measurement-input/sensor_Mtr_119123112.csv")==false
    current_measurement[4]=0;
  else

  sensor_Mtr_382883 =CSV.read("$(path)/input/measurement-input/sensor_Mtr_-382883.csv"; comment="#");
  sensor_Mtr_1695832 =CSV.read("$(path)/input/measurement-input/sensor_Mtr_1695832.csv"; comment="#");
  sensor_Mtr_118122225 =CSV.read("$(path)/input/measurement-input/sensor_Mtr_118122225.csv"; comment="#");
  sensor_Mtr_119123112 =CSV.read("$(path)/input/measurement-input/sensor_Mtr_119123112.csv"; comment="#");

  # check if data at i-second instant is available; if not then wait until it is available
  while (length(sensor_Mtr_382883[:,1])<(i+300*(times-1))) & (length(sensor_Mtr_1695832[:,1])<(i+300*(times-1))) & (length(sensor_Mtr_118122225[:,1])<(i+300*(times-1))) & (length(sensor_Mtr_119123112[:,1])<(i+300*(times-1)))
      sleep(0.001);
  end

  # calculate control input based on measurement at i-second instant

  if typeof(sensor_Mtr_382883[2,2])==string
  # convert from string to complex number measurement data
  imcurrent_measurement=zeros(ComplexF64,4,1);
  imcurrent_measurement[1]= (parse(ComplexF64, sensor_Mtr_382883[(i+300*(times-1)),2])+parse(ComplexF64, sensor_Mtr_382883[(i+300*(times-1)),3])+parse(ComplexF64, sensor_Mtr_382883[(i+300*(times-1)),4]))/3;
  imcurrent_measurement[2]= (parse(ComplexF64, sensor_Mtr_1695832[(i+300*(times-1)),2])+parse(ComplexF64, sensor_Mtr_1695832[(i+300*(times-1)),3])+parse(ComplexF64, sensor_Mtr_1695832[(i+300*(times-1)),4]))/3;
  imcurrent_measurement[3]= (parse(ComplexF64, sensor_Mtr_118122225[(i+300*(times-1)),2])+parse(ComplexF64, sensor_Mtr_118122225[(i+300*(times-1)),3])+parse(ComplexF64, sensor_Mtr_118122225[(i+300*(times-1)),4]))/3;
  imcurrent_measurement[4]= (parse(ComplexF64, sensor_Mtr_119123112[(i+300*(times-1)),2])+parse(ComplexF64, sensor_Mtr_119123112[(i+300*(times-1)),3])+parse(ComplexF64, sensor_Mtr_119123112[(i+300*(times-1)),4]))/3;
     current_measurement[1]=(1.05-0.1*rand(Float64))*abs(imcurrent_measurement[1]);
    for j=2:4
      current_measurement[j]= abs(imcurrent_measurement[j]);
    end

  elseif typeof(sensor_Mtr_382883[2,2])==Int64
  current_measurement[1]= (1.05-0.1*rand(Float64))*sqrt(sensor_Mtr_382883[(i+300*(times-1)),8]^2+sensor_Mtr_382883[(i+300*(times-1)),9]^2);
  current_measurement[2]= sqrt(sensor_Mtr_1695832[(i+300*(times-1)),8]^2+sensor_Mtr_1695832[i(i+300*(times-1)),9]^2);
  current_measurement[3]= sqrt(sensor_Mtr_118122225[(i+300*(times-1)),8]^2+sensor_Mtr_118122225[(i+300*(times-1)),9]^2);
  current_measurement[4]= sqrt(sensor_Mtr_119123112[(i+300*(times-1)),8]^2+sensor_Mtr_119123112[(i+300*(times-1)),9]^2);
  end
  end

  ###############################
  # find the next observer state
  global current_observer
  next_observer = current_observer + (A+B*K_out)*current_observer + L_out*(current_measurement-C*current_observer);

  current_observer =next_observer;

  ################################
  # Step 3: return PQ control signals bn
  current_control = K_out*current_observer;
  PQcurrent_control=zeros(ComplexF64,8,1);
  for k=1:8
  PQcurrent_control[k] = current_control[2*k-1] + current_control[2*k]*im;
  end
  ################################
  # Step 4: calculate the changes to apply to control nodes

  dynamic_injections = optimal_injections + PQcurrent_control;

  # write the file with the stringdata variable information
  global df = DataFrame(Dict("datatime" => DateTime(2000,1,1,0,floor((i+300*(times-1))/60),(i+300*(times-1))-60*floor((i+300*(times-1))/60)), "generator.proposed_diesel_1.1" =>  [dynamic_injections[1]], "generator.proposed_diesel_1.2" =>  [dynamic_injections[1]], "generator.proposed_diesel_1.3" =>  [dynamic_injections[1]], "generator.ouc_solar.1" =>  [dynamic_injections[2]],"generator.ouc_solar.2" =>  [dynamic_injections[2]],"generator.ouc_solar.3" =>  [dynamic_injections[2]],"generator.gen_none_-60027406210.1" =>  [dynamic_injections[3]],"generator.gen_none_-60027406210.2" => [dynamic_injections[3]], "generator.gen_none_-60027406210.3" =>  [dynamic_injections[3]],"generator.proposed_solarpv_1.1" =>  [dynamic_injections[4]],"generator.proposed_solarpv_1.2" =>  [dynamic_injections[4]],"generator.proposed_solarpv_1.3" =>  [dynamic_injections[4]],"generator.gen_none_6002693162.1" =>  [dynamic_injections[5]],"generator.gen_none_6002693162.2" =>  [dynamic_injections[5]],"generator.gen_none_6002693162.3" =>  [dynamic_injections[5]],"storage.proposed_storage_1.1" =>  [dynamic_injections[6]],"storage.proposed_storage_1.2" =>  [dynamic_injections[6]],"storage.proposed_storage_1.3" =>  [dynamic_injections[6]],"storage.proposed_fuelcell_1.1" =>  [dynamic_injections[7]],"storage.proposed_fuelcell_1.2" =>  [dynamic_injections[7]],"storage.proposed_fuelcell_1.3" =>  [dynamic_injections[7]],"storage.proposed_fuelcell_2.1" =>  [dynamic_injections[8]],"storage.proposed_fuelcell_2.2" =>  [dynamic_injections[8]],"storage.proposed_fuelcell_2.3" =>  [dynamic_injections[8]]));
  # save pnnl dynamic actions to the output folder
  CSV.write("$(path)/output/pnnl_actions.csv", df,append=true);

  end
  end

elseif case==3 #failure actuator

  for times = 1:number_of_energizing_segment
  for i=1:300

      # remember to update actions when running for > 5 minute

  # Step1: Get measurement data from measurement input folder
  #sensor_Mtr_382883 =readtable("/Applications/projects/cLEAN/pnnl control/measurement input/sensor_Mtr_-382883.csv")
  #df1 = CSV.File("/Applications/projects/cLEAN/pnnl control/measurement input/sensor_Mtr_-382883.csv", datarow=10);
  #sensor_Mtr_382883 =CSV.read(df1)
  current_measurement=zeros(4,1);
  if isfile("$(path)/input/measurement-input/sensor_Mtr_-382883.csv")==false
    current_measurement[1]=0;
  elseif isfile("$(path)/input/measurement-input/sensor_Mtr_1695832.csv")==false
    current_measurement[2]=0;
  elseif isfile("$(path)/input/measurement-input/sensor_Mtr_118122225.csv")==false
    current_measurement[3]=0;
  elseif isfile("$(path)/input/measurement-input/sensor_Mtr_119123112.csv")==false
    current_measurement[4]=0;
  else

  sensor_Mtr_382883 =CSV.read("$(path)/input/measurement-input/sensor_Mtr_-382883.csv"; comment="#");
  sensor_Mtr_1695832 =CSV.read("$(path)/input/measurement-input/sensor_Mtr_1695832.csv"; comment="#");
  sensor_Mtr_118122225 =CSV.read("$(path)/input/measurement-input/sensor_Mtr_118122225.csv"; comment="#");
  sensor_Mtr_119123112 =CSV.read("$(path)/input/measurement-input/sensor_Mtr_119123112.csv"; comment="#");

  # check if data at i-second instant is available; if not then wait until it is available
  while (length(sensor_Mtr_382883[:,1])<(i+300*(times-1))) & (length(sensor_Mtr_1695832[:,1])<(i+300*(times-1))) & (length(sensor_Mtr_118122225[:,1])<(i+300*(times-1))) & (length(sensor_Mtr_119123112[:,1])<(i+300*(times-1)))
      sleep(0.001);
  end

  # calculate control input based on measurement at i-second instant

  if typeof(sensor_Mtr_382883[2,2])==string
  # convert from string to complex number measurement data
  imcurrent_measurement=zeros(ComplexF64,4,1);
  imcurrent_measurement[1]= (parse(ComplexF64, sensor_Mtr_382883[(i+300*(times-1)),2])+parse(ComplexF64, sensor_Mtr_382883[(i+300*(times-1)),3])+parse(ComplexF64, sensor_Mtr_382883[(i+300*(times-1)),4]))/3;
  imcurrent_measurement[2]= (parse(ComplexF64, sensor_Mtr_1695832[(i+300*(times-1)),2])+parse(ComplexF64, sensor_Mtr_1695832[(i+300*(times-1)),3])+parse(ComplexF64, sensor_Mtr_1695832[(i+300*(times-1)),4]))/3;
  imcurrent_measurement[3]= (parse(ComplexF64, sensor_Mtr_118122225[(i+300*(times-1)),2])+parse(ComplexF64, sensor_Mtr_118122225[(i+300*(times-1)),3])+parse(ComplexF64, sensor_Mtr_118122225[(i+300*(times-1)),4]))/3;
  imcurrent_measurement[4]= (parse(ComplexF64, sensor_Mtr_119123112[(i+300*(times-1)),2])+parse(ComplexF64, sensor_Mtr_119123112[(i+300*(times-1)),3])+parse(ComplexF64, sensor_Mtr_119123112[(i+300*(times-1)),4]))/3;
    for j=1:4
      current_measurement[j]= abs(imcurrent_measurement[j]);
    end

  elseif typeof(sensor_Mtr_382883[2,2])==Int64
  current_measurement[1]= sqrt(sensor_Mtr_382883[(i+300*(times-1)),8]^2+sensor_Mtr_382883[(i+300*(times-1)),9]^2);
  current_measurement[2]= sqrt(sensor_Mtr_1695832[(i+300*(times-1)),8]^2+sensor_Mtr_1695832[i(i+300*(times-1)),9]^2);
  current_measurement[3]= sqrt(sensor_Mtr_118122225[(i+300*(times-1)),8]^2+sensor_Mtr_118122225[(i+300*(times-1)),9]^2);
  current_measurement[4]= sqrt(sensor_Mtr_119123112[(i+300*(times-1)),8]^2+sensor_Mtr_119123112[(i+300*(times-1)),9]^2);
  end
  end

  ###############################
  # find the next observer state
  global current_observer;
  next_observer = current_observer + (A+B*K_out)*current_observer + L_out*(current_measurement-C*current_observer);

  current_observer =next_observer;

  ################################
  # Step 3: return PQ control signals bn
  current_control = K_out*current_observer;
  PQcurrent_control=zeros(ComplexF64,8,1);
  for k=1:8
  PQcurrent_control[k] = current_control[2*k-1] + current_control[2*k]*im;
  end

  PQcurrent_control[1]=0; #actuators at generator 1 fail

  ################################
  # Step 4: calculate the changes to apply to control nodes

  dynamic_injections = optimal_injections + PQcurrent_control;

  # write the file with the stringdata variable information
  global df = DataFrame(Dict("datatime" => DateTime(2000,1,1,0,floor((i+300*(times-1))/60),(i+300*(times-1))-60*floor((i+300*(times-1))/60)), "generator.proposed_diesel_1.1" =>  [dynamic_injections[1]], "generator.proposed_diesel_1.2" =>  [dynamic_injections[1]], "generator.proposed_diesel_1.3" =>  [dynamic_injections[1]], "generator.ouc_solar.1" =>  [dynamic_injections[2]],"generator.ouc_solar.2" =>  [dynamic_injections[2]],"generator.ouc_solar.3" =>  [dynamic_injections[2]],"generator.gen_none_-60027406210.1" =>  [dynamic_injections[3]],"generator.gen_none_-60027406210.2" => [dynamic_injections[3]], "generator.gen_none_-60027406210.3" =>  [dynamic_injections[3]],"generator.proposed_solarpv_1.1" =>  [dynamic_injections[4]],"generator.proposed_solarpv_1.2" =>  [dynamic_injections[4]],"generator.proposed_solarpv_1.3" =>  [dynamic_injections[4]],"generator.gen_none_6002693162.1" =>  [dynamic_injections[5]],"generator.gen_none_6002693162.2" =>  [dynamic_injections[5]],"generator.gen_none_6002693162.3" =>  [dynamic_injections[5]],"storage.proposed_storage_1.1" =>  [dynamic_injections[6]],"storage.proposed_storage_1.2" =>  [dynamic_injections[6]],"storage.proposed_storage_1.3" =>  [dynamic_injections[6]],"storage.proposed_fuelcell_1.1" =>  [dynamic_injections[7]],"storage.proposed_fuelcell_1.2" =>  [dynamic_injections[7]],"storage.proposed_fuelcell_1.3" =>  [dynamic_injections[7]],"storage.proposed_fuelcell_2.1" =>  [dynamic_injections[8]],"storage.proposed_fuelcell_2.2" =>  [dynamic_injections[8]],"storage.proposed_fuelcell_2.3" =>  [dynamic_injections[8]]));
  # save pnnl dynamic actions to the output folder
  CSV.write("$(path)/output/pnnl_actions.csv", df,append=true);

  end
  end

end
