clc
clear all
clear Rosenbrock GR_Rosenbrock H_Rosenbrock Powel GR_Powel H_Powel 
%%
x0 = [1;2];
[x_min, f_min, iter] = Newton_GSS(@Rosenbrock, @GR_Rosenbrock, @H_Rosenbrock, x0, 1e-3, 1e-5)
[~,NEWTONRosenEvals]=Rosenbrock(x0);
NEWTONRosenEvals = NEWTONRosenEvals-1;
[~,NEWTON_GR_RosenEvals]=GR_Rosenbrock(x0);
NEWTON_GR_RosenEvals = NEWTON_GR_RosenEvals-1;
[~,NEWTON_H_RosenEvals]=H_Rosenbrock(x0);
NEWTON_H_RosenEvals = NEWTON_H_RosenEvals-1;
clear Rosenbrock GR_Rosenbrock H_Rosenbrock Powel GR_Powel H_Powel 
%%
x0 = [1;2;2;2];
[x_min, f_min, iter] = Newton_GSS(@Powel, @GR_Powel, @H_Powel, x0, 1e-3, 1e-5)
[~,NEWTONPowelEvals]=Powel(x0);
NEWTONPowelEvals = NEWTONPowelEvals-1;
[~,NEWTON_GR_PowelEvals]=GR_Powel(x0);
NEWTON_GR_PowelEvals = NEWTON_GR_PowelEvals-1;
[~,NEWTON_H_PowelEvals]=H_Powel(x0);
NEWTON_H_PowelEvals = NEWTON_H_PowelEvals-1;
%%
clear Rosenbrock GR_Rosenbrock H_Rosenbrock Powel GR_Powel H_Powel 

x0 = [1;2];
[x_min, f_min, iter] = SD_GSS(@Rosenbrock, @GR_Rosenbrock, [1;2], 1e-3, 1e-5)
[~,SDRosenEvals]=Rosenbrock(x0);
SDRosenEvals = SDRosenEvals-1;
[~,SD_GR_RosenEvals]=GR_Rosenbrock(x0);
SD_GR_RosenEvals = SD_GR_RosenEvals-1;
[~,SD_H_RosenEvals]=H_Rosenbrock(x0);
SD_H_RosenEvals = SD_H_RosenEvals-1;
clear Rosenbrock GR_Rosenbrock H_Rosenbrock Powel GR_Powel H_Powel 

x0 = [1;2;2;2];
[x_min, f_min, iter] = SD_GSS(@Powel, @GR_Powel, x0, 1e-3, 1e-5)
[~,SDPowelEvals]=Powel(x0);
SDPowelEvals = SDPowelEvals-1;
[~,SD_GR_PowelEvals]=GR_Powel(x0);
SD_GR_PowelEvals = SD_GR_PowelEvals-1;
[~,SD_H_PowelEvals]=H_Powel(x0);
SD_H_PowelEvals = SD_H_PowelEvals-1;
