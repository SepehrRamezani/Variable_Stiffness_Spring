
% folder='C:\MyCloud\OneDriveUcf\Real\VariableSpring\Setup_Test';
% close all
% clear all
clc
lpvar=5;
ForceCo=[0.1491 0];
PressureCo=[6.6467 52.761];
rng default;
% glob Ap
% Ap=
% [file,path]=uigetfile('C:\MyCloud\OneDriveUcf\Real\VariableSpring\Setup_Test\*.csv');
o=0;
tt=0;
folderpath = 'C:\MyCloud\OneDriveUcf\Real\VariableSpring\Setup_Test\';
% folderpath = fullfile(folderpath, '**');    % What is the meaning of "/**/" ???
filelist   = dir(folderpath);
name       = {filelist.name};
name=name(strncmp(name, 'd7', 2)| strncmp(name, 'd16', 3)| strncmp(name, 'd30', 3));   % No files starting with '.'
name=erase(name, '.csv');
% name={'d7cc1','d7cc2','d7cc3','d16cc1','d16cc2','d16cc3','d30cc1','d30cc2','d30cc3'};
Name1={'0','5','10','15'};
Name2={'7','16','30'};
name={'s25p5v7','s25p5v16','s25p5v30','s25p10v7','s25p10v16','s25p10v30','s25p15v7','s25p15v16','s25p15v30'};
legendname=[];
figure
for k=2:length(Name1)
    for u=2:length(Name2)
     name=append("s25p",Name1(k),"v",Name2(u));   
    tt=tt+1;
%     if o==3
%          legend(num2str(Vel'),'Location','northwest')
%         
%         title(name(k-1));
%         figure;
%         o=0;
%        
%     end
  o=o+1;
    %% Imorting
    data=importdata(append(folderpath,name,".csv"));

    data.data(:,3)=data.data(:,3)-11.1;
    DStime=data.data(1,6);
    samfreq=1/DStime;
    %% filtering and crupping data
    [bb,aa] = butter(4, lpvar/(samfreq/2),'low');
    FilteredData=filtfilt(bb,aa,data.data(:,3));
    indx=(FilteredData<-0.1 & FilteredData >-6.5 );
    timeindx=find(indx);
    timeindx2=timeindx(timeindx>100);
    indx2=[timeindx2(1):1:timeindx2(end)]';
%        figure
%         hax=axes;
%         plot([FilteredData,data.data(:,[2,3,4])]);
%         title(name(k))
%         hold on
%         line([indx2(1) indx2(1)],get(hax,'YLim'),'Color',[.5 0.5 0])
%         line([indx2(end) indx2(end)],get(hax,'YLim'),'Color',[.5 0.5 0])
%         hold off
    
    newData=data.data(indx2,:);
    
    %% Position
    M_initial= mean(data.data(1:indx2-10,3));
    t=(newData(:,1)-newData(1,1)).*DStime;
    M=-1.*newData(:,3)./(11.1-4.5).*20.3;
    fun = @(w)LinearReg(w,t,M);
    x0 = 1;
    bestx = fminsearch(fun,x0);
%     Mc = polyfit(t,M,1);
%     M_l = polyval(Mc,t);
    M_l2=bestx(1)*t;

    %%velocity calculation
    Vel(tt)=bestx;
%     figure
    %     plot(t,M,t,M_l,t,M_l2); 
    %
    %     ylabel('Position (mm)')
    %     title(['Piston Position at Veocity ',num2str(Vel),' mm/s'])
    
    %% Force
    x=M_l2;
    F_initial= mean(data.data(1:indx2-10,4));
    F=-1*(newData(:,2)-F_initial);
    fun = @(w)sseval(w,x,F);
    x0 = [40;0.5;1;1;1];
    options = optimset('MaxFunEvals',100000,'MaxIter',100000);
    bestx = fminsearch(fun,x0,options)
    F_l2_cal=CurveFun(bestx,x).*149.08;
    Fc = polyfit(x,F,3);
%     F_l = polyval(Fc,x);
    F_cal= (F)*149.08;
%     F_l_cal=(F_l)*149.08;
%      figure
     yyaxis left
     g=plot(x, F_cal);
     g.Color(4)=0.2;
     ylabel('Force (N)');
     hold on
     plot(x,F_l2_cal);
     plot(s25p5v16.Displacment(:,2)*-1000,s25p5v16.Force,'LineWidth',2);
     
     yyaxis right
    
    
    %% pressure
%     P=1036.5*newData(:,4)*1000-54603;

    P=(newData(:,4)-0.055);
    fun = @(w)sseval(w,x,P);
    x0 = [40;0.5;1;1;1];
    options = optimset('MaxFunEvals',100000,'MaxIter',100000);
    bestxp = fminsearch(fun,x0,options);
    P_l2_cal=CurveFun(bestxp,x).*150.33;
    
%     Pc = polyfit(x,P,3);
%     P_l = polyval(Pc,x).*150.33;
%     figure
    g2=plot(x,P.*150.33);
    g2.Color(4)=0.2;
    plot(x,P_l2_cal);
    ylabel('Pressure(psi)');
    plot(s25p5v16.Displacment(:,2)*-1000,s25p5v16.P_T.signals(1).values.*0.000145038,'LineWidth',2);
%     plot(x,P_l2_cal,'LineWidth',2);
%     hold on
    %     ylabel('Pressure (Pa)')
    %     title(['Chamber Pressure Veocity ',num2str(Vel),' mm/s'])
%  legendname=[legendname,name(k)];
%% Stifness
% figure

K=diff(F_l2_cal)./diff(x)*1000;
% Kc = polyfit(x(2:end),K,1);
% P_l = polyval(Kc,x(2:end));
mdl = fitlm(x(2:end),K);
R(k,u)=mdl.Rsquared.Ordinary;
% plot(x(2:end),K,'LineWidth',2);

hold on
    end
    title(append('Stiffness at ',Name1(k),'psi'))
    legend(Name2,'Location','northwest');
    ylabel('K (N/m)');
    xlabel('Stroke (mm)');
    figure
end
ylabel('Force (N)');
xlabel('Stroke (mm)');
legend(legendname,'Location','northwest');
% title(['Dynamic test at ',num2str(Vel(1)),' mm/s'])
% title(name(end));
% legend(num2str(Vel'),'Location','northwest')
function sse = sseval(x,xdata,ydata)
A = x(1);
B = x(2);
C = x(3);
D = x(4);
E = x(5);
sse = sum((ydata -C.*(A./(A-D*xdata)).^B+E).^2);
end
function sse = LinearReg(x,xdata,ydata)
A = x(1);
sse = sum((ydata -A.*xdata).^2);
end
function ydata = CurveFun(x,xdata)
A = x(1);
B = x(2);
C = x(3);
D = x(4);
E=x(5);
ydata =C.*(A./(A-xdata*D)).^B-E;
end
