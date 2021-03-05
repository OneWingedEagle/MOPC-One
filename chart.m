function charter

%filename= 'D:\Works and Studies\Photonic\Lumerical\models\faradary_isolator\results-7-180s-108.txt';
[filename,filepath]=uigetfile('*.txt', 'Select results file')
file=strcat(filepath,strcat('\',filename));
data=csvread(file,3,0);
     

Fn=data(:,1);

ndata=size(Fn,1);

wn1=Fn(1);
wn2=Fn(ndata);

colT='-k';
colR='-k';


 figure(1)

 plot(Fn(:,1),data(:,3),colR,'LineWidth',2);
  hold on
 axis([wn1,wn2,0,1]);
 
 figure(2)
 plot(Fn(:,1),data(:,2),colT,'LineWidth',2);


 tmax=1.1*max(abs(data(:,2)));
   axis([wn1,wn2,-tmax,tmax]);

  hold on
end