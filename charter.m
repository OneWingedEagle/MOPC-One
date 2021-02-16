function charter
  
  ndata=6;
%[filename1,filepath1]=uigetfile('*.txt', 'Selectinput file')
% 
 filepath1='D:\Works and Studies\Photonic\Lumerical\models\mopc16';
 filename1='D:\Works and Studies\Photonic\Lumerical\models\mopc16\Trbb7.txt';
 filename2='D:\Works and Studies\Photonic\Lumerical\models\mopc16\frbb.txt';
 filename3='D:\Works and Studies\Photonic\Lumerical\models\mopc16\results7_180x.txt';
 filename4='D:\Works and Studies\Photonic\Lumerical\models\mopc16\Trbb4.txt';
  filename5='D:\Works and Studies\Photonic\Lumerical\models\mopc16\Ex.txt';
    filename6='D:\Works and Studies\Photonic\Lumerical\models\mopc16\Ez.txt';
    
    dataEx=csvread(filename5,3,0);
    dataEz=csvread(filename6,3,0);
     
data=csvread(filename1,3,0);

data4=csvread(filename4,3,0);
ndata=size(data,1);

data2=csvread(filename2,3,0);
data3=csvread(filename3,3,0);

lams=data(:,1);

vals=data(:,2);

fr=data2(:,2);

a=1.17e-6;

Fn=a./lams;

wn1=Fn(1)
wn2=Fn(ndata);

colT='-r';
colR='-k';

 angs=dataEx(:,1);
#angs=180/pi*atan(dataEx(:,2)./dataEz(:,2));
for k=1:length(angs)
  Ex=dataEx(k,2);
  Ez=dataEz(k,2);
  TM_frac = abs(Ex)^2/(abs(Ex)^2+abs(Ez)^2);   
  A = sqrt(TM_frac);

angs(k)= 90-acos(A)/pi*180;

%angs(k)= atan(Ex/Ez)/pi*180;
end

 figure(1)
 %plot(Fn,vals,colT,'LineWidth',2);

  plot(dataEx(:,1),angs,colR,'LineWidth',2);

 %axis([wn1,wn2,0,1]);
 
 %hold on
  %  hold onplot(Fn,data4(:,2),'-g','LineWidth',2);

  %plot(data3(:,1),data3(:,2),colR,'LineWidth',2);
 % plot(Fn,fr/90,colR,'LineWidth',2);
  hold on
end