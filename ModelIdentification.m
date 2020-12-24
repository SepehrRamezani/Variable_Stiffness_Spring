
% folder='C:\MyCloud\OneDriveUcf\Real\VariableSpring\Setup_Test';
 close all
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
name={'d7cc1','d7cc2','d7cc3','d16cc1','d16cc2','d16cc3','d30cc1','d30cc2','d30cc3'};
for k=1:1:3
    tt=tt+1
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
    data=importdata(append(folderpath,name(k),".csv"));

    data.data(:,3)=data.data(:,3)-11.1;
    DStime=data.data(1,6);
    samfreq=1/DStime;
    %% filtering and crupping data
    [bb,aa] = butter(4, lpvar/(samfreq/2),'low');
    FilteredData=filtfilt(bb,aa,data.data(:,3));
    indx=(FilteredData<-0.05 & data.data(:,3)>-6.1 );
    timeindx=find(indx);
    indx2=timeindx(timeindx>175);
%        figure
%         hax=axes;
%         plot(data.data(:,[2,3,4]));
%         title(name(k))
%         hold on
%         line([indx2(1) indx2(1)],get(hax,'YLim'),'Color',[.5 0.5 0])
%         line([indx2(end) indx2(end)],get(hax,'YLim'),'Color',[.5 0.5 0])
    %     hold off
    %
    newData=data.data(indx2,:);
    
    %% Position
    t=(newData(:,1)-newData(1,1)).*DStime;
    M=-1.*newData(:,3)./(11.1-4.75).*20.3;
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
    F=-1*newData(:,2);
    fun = @(w)sseval(w,x,F);
    x0 = [40;0.5;1;1];
    
    options = optimset('MaxFunEvals',100000,'MaxIter',100000);
    bestx = fminsearch(fun,x0,options)
    F_l2_cal=CurveFun(bestx,x).*149.08;
    Fc = polyfit(x,F,3);
    F_l = polyval(Fc,x);
    F_cal= (F-F_l(1))*149.08;
    F_l_cal=(F_l-F_l(1))*149.08;
%         figure
      plot(x,[F_l2_cal],'LineWidth',2);
%         plot(x,F_cal,x,F_l_cal,x,F_l2_cal,'LineWidth',2);
%     plot(M_l_cal,F_l2_cal,'LineWidth',2);
     hold on
    
    
    
    %% pressure
%     P=1036.5*newData(:,4)*1000-54603;
    P=newData(:,4)-0.053;
    fun = @(w)sseval(w,x,P);
    x0 = [40;0.5;1;1];
     options = optimset('MaxFunEvals',100000,'MaxIter',100000);
    bestxp = fminsearch(fun,x0,options);
    P_l2_cal=CurveFun(bestxp,x).*1036.5;
    
    Pc = polyfit(x,P,3);
    P_l = polyval(Pc,x).*1036.5;
%     figure
%     plot(x,P,x,P_l2_cal,x,P_l,'LineWidth',2);
%     plot(x,P_l2_cal,'LineWidth',2);
%     hold on
    %     ylabel('Pressure (Pa)')
    %     title(['Chamber Pressure Veocity ',num2str(Vel),' mm/s'])
 
  
%     plot((name(1)).Displacment(:,2)*-1000,(name(1)).Force(:,2));
end
ylabel('Force (N)')
xlabel('Stroke (mm)')
title(['Dynamic test at ',num2str(Vel(1)),' mm/s'])
% title(name(end));
% legend(num2str(Vel'),'Location','northwest')
function sse = sseval(x,xdata,ydata)
A = x(1);
B = x(2);
C = x(3);
D = x(4);
sse = sum((ydata -C.*(A./(A-D*xdata)).^B+C).^2);
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
ydata =C.*(A./(A-xdata*D)).^B-C;
end
