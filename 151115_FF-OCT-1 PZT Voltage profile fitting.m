clear all
cd('D:\MATLAB\141211_PZT waveform file generation');

%%
Data=dlmread('Waveform_f0.004_V10_O1.txt');        
plot(Data(:,1),Data(:,2));

Original_Velocity=1;

Time=Data(:,1);
Time(500000)=[];
Voltage=Data(:,2);
Voltage(end)=[];

plot(Time,Voltage);

%% chnage to ms

Time=Time*1000;

%% Change speed

V_Increasing_Ratio=0.812;

Time_New=Time/V_Increasing_Ratio;

%% Fitting

 g = fittype( @(a, b, c, x) a*x+b*(x.^0.5)+c);
 
 Fit_Result = fit(Time_New(1:floor(length(Time_New)/2)),Voltage(1:floor(length(Voltage)/2)),g);
 
 Voltage_fit=Fit_Result.a*Time_New+Fit_Result.b*(Time_New.^0.5)+Fit_Result.c;
 
 plot(Time_New,Voltage,Time_New,Voltage_fit);
Fit_Result.a
Fit_Result.b
Fit_Result.c

%%  filtering the signal