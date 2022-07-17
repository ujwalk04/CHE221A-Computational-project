% Writing known Parameters
p=0.751;    % Total Pressure
% Antoine's Equation Parameters
A1=8.07131; 
B1=1730.630;
C1=233.426;
A2=7.16030;
B2=1447.700;
C2=210.000;
% Threshold values
eps=1e-4;
dT=1e-3;
% Wilson Parameters
A12=1675.0732;
A21=-1508.1875;
% Volume of pure liquid components
v1=18.07;  %water
v2=87.51;  %morpholine
% Initialising x1, y1 & Temp
x1=ones(50,1);
y1=ones(50,1);
Temp=ones(50,1);

x1(1)=0.02;
for i=2:50
    x1(i)=x1(i-1)+0.02;
end

T_sat1=(B1/(A1-log10(p*750.062)))-C1+273;
T_sat2=(B2/(A2-log10(p*750.062)))-C2+273;
% Calculating for each x1
for i=1:50
    x11=x1(i); % x1
    x12=1-x1(i); % x2=1-x1
    Temp(i)=T_sat1*x11+T_sat2*x12; 
    % Wilson parameters calculation
    lamda12=(v2/v1)*exp(-A12/(1.9872*Temp(i)));
    lamda21=(v1/v2)*exp(-A21/(1.9872*Temp(i)));
    gamma1=exp(-log(x11+(lamda12*x12))+x12*((lamda12/(x11+(lamda12*x12)))-(lamda21/(x12+(lamda21*x11)))));
    gamma2=exp(-log(x12+(lamda21*x11))+x11*((lamda21/(x12+(lamda21*x11)))-(lamda12/(x11+(lamda12*x12)))));
    P_sat2=0.00133322*(10^(A2-(B2/(C2+(Temp(i)-273)))));
    P_sat1=0.00133322*(10^(A1-(B1/(C1+(Temp(i)-273)))));
    P_t=P_sat1*gamma1*x11+P_sat2*gamma2*x12;
    y11=(gamma1*x11*P_sat1)/P_t;
    y12=(gamma2*x12*P_sat2)/P_t;

    while(abs((y11+y12)-1)>eps)

        if(y11+y12<1)    %Less vapour, hence increase Temperature
            Temp(i)=Temp(i)+dT;
        else             %More vapour, hence decrease Temperature
            Temp(i)=Temp(i)-dT;
        end
    
    lamda12=(v2/v1)*exp(-A12/(1.9872*Temp(i)));
    lamda21=(v1/v2)*exp(-A21/(1.9872*Temp(i)));
    gamma1=exp(-log(x1+(lamda12*x12))+x12*((lamda12/(x11+(lamda12*x12)))-(lamda21/(x12+(lamda21*x11)))));
    gamma2=exp(-log(x2+(lamda21*x11))+x11*((lamda21/(x12+(lamda21*x11)))-(lamda12/(x11+(lamda12*x12)))));
    P_t=P_sat1*gamma1*x11+P_sat2*gamma2*x12;
    y11=(gamma1*x11*P_sat1)/P_t;
    y12=(gamma2*x12*P_sat2)/P_t;

    end
y1(i)=y11;
end

% Reading experimental data
f="experimental.txt";
M=readmatrix(f);
Texp=M(:,1);
x_exp=M(:,2);
y_exp=M(:,3);
Texp=Texp+273;

% Plotting T vs x,y
figure
scatter(x_exp,Texp);
hold on;
scatter(y_exp,Texp);
hold on;
plot(x1,Temp,'-r');
xlim([0,1]);
hold on;
plot(y1,Temp,'-m')
xlabel('Mole fractions');
ylabel('Temp.(K)');
legend('x-exp','y-exp','x1','y1');
hold off;

% Plotting x vs y
figure
plot(x1,y1);
xlim([0,1]);
hold on;
scatter(x_exp,y_exp);
xlabel('x');
ylabel('y');
legend('y-calc','y-exp');
hold off;
