%%%% Solves photometry decay equation

%clear t y

syms t y ;

% xa=[];
% ya=[];
% xb=[];
% yb=[];
% xc=[];
% yc=[];


xa=offsetNum;
ya=median(exampleTrace(offsetNum-medianRange:offsetNum+medianRange));

midTrace=floor(size(exampleTrace,1)/2);

xb=midTrace;
yb=median(exampleTrace(midTrace-medianRange:midTrace+medianRange));

endTrace=size(exampleTrace,1)-offsetNum;

xc=endTrace;
yc=median(exampleTrace(endTrace-medianRange:endTrace+medianRange));


eqn = t == -xc/(log(yc - (yb-ya*exp((xa-xb)/t)) / (1-exp((xa-xb)/t)))- log(ya - (yb-ya*exp((xa-xb)/t)) / (1-exp((xa-xb)/t)))-xa/t);

S = vpasolve(eqn,t);

%% Solving y0;

t1=double(S);

eqny0 = y == yb - (ya - y)*exp((xa-xb)/t1);

Sy0= vpasolve(eqny0, y);

y0=double(Sy0);


%% Calculating A1;

A1 = (ya-y0)*exp(xa/t1);


%% Plotting the curves
%figure()

%x=0:1/10:129;

x=1:1:size(exampleTrace,1);

if isempty(A1)

yy=nan(1,size(exampleTrace,1));

t1=nan;

else
    
yy=A1*exp(-x/t1)+y0;

end

plot(x,yy,'-m')
hold on
plot(x,exampleTrace,'-k')


