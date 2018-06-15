clear all
cd('D:\Users\TuanShu\141211_PZT waveform file generation');

%%
Data=dlmread('Waveform_f0.004_V10_O1.txt');        
plot(Data(:,1),Data(:,2));

Original_Velocity=1;

Time=Data(:,1);
Time(500000)=[];
Voltage=Data(:,2);
Voltage(end)=[];

plot(Time,Voltage);

[value index]=max(Voltage);


Voltage_Foward=Voltage(1:index);
Time_Foward=Time(1:index);
Voltage_Backward=Voltage(index+1:end);
Time_Backward=Time(index+1:end);
COEF_Foward=fit(Time_Foward,Voltage_Foward,'poly2');
Voltage_Foward_Fit=COEF_Foward.p1*Time_Foward.^2 + COEF_Foward.p2*Time_Foward + COEF_Foward.p3;

plot(Time_Foward,Voltage_Foward,Time_Foward,Voltage_Foward_Fit);
%%  filtering the signal
Ratio=2;      %i.e. 如果為0.5, 那就是減速一半

D_Time=Data(2,1)-Data(1,1);
Time_Max_Old=max(Time);
Time_Max_New=round(Time_Max_Old/Ratio);

Time_New_Waveform=D_Time:D_Time:Time_Max_New;

Voltage_New_Waveform=interp1(Time,Voltage,Time_New_Waveform*Ratio,'linear','extrap');

Voltage_New_Waveform(Voltage_New_Waveform<-1.999)=-1.999;

plot(Time_New_Waveform,Voltage_New_Waveform,'linewidth',2);
xlabel('Time (second)','fontsize',12);
ylabel('Voltage (V)','fontsize',12);

Output=[Time_New_Waveform' Voltage_New_Waveform'];

dlmwrite(sprintf('Waveform_SR%dHz_V%0.1dmicronsec_INTERPOLED_withbackwalking.txt',1/D_Time,Original_Velocity*Ratio),Output,'delimiter','\t','newline','pc','precision', '%.6f');
