function plot4me
  
[filename,filepath]=uigetfile('*.txt', 'Select results file')
file=strcat(filepath,strcat('\',filename));
data=csvread(file,3,0);


ndata=size(data,1);

Fn(:,1)=data(:,1);

wn1=Fn(1);
wn2=Fn(ndata);

colT='-r';
colR='-r';



  
if(size(data,2)==4)


figure(1)

Tz(:,1)=data(:,3)-data(:,4);
 %%%  Tz Transmittance
 plot(Fn(:,1),Tz(:,1),'-ob','LineWidth',2);
 
   axis([wn1,wn2,0,1]);
   
  hold on

 %%%  Tx Transmittance
% plot(Fn(:,1),data(:,4),'-or','LineWidth',2);
%  hold on
  

  
else 
  figure(1)

 %%%  Transmittance
 plot(Fn(:,1),data(:,3),'-ob','LineWidth',2);
 
    axis([wn1,wn2,0,1]);
    
  hold on 
  
end
 
 figure(3)
  %%%  Faraday Rotation
 plot(Fn(:,1),data(:,2),'-ob','LineWidth',2);

  axis([wn1,wn2,-90,90]);

  hold on
end