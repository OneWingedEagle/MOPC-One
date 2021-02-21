function charter
  
  ndata=6;
%[filename1,filepath1]=uigetfile('*.txt', 'Selectinput file')
% 
 filepath1='D:\Works and Studies\Photonic\Lumerical\models\mopc16';
 filename1='D:\Works and Studies\Photonic\Lumerical\models\mopc16\Trbb7.txt';
 filename2='D:\Works and Studies\Photonic\Lumerical\models\mopc16\fr.txt';
 filename3='D:\Works and Studies\Photonic\Lumerical\models\mopc16\fr_results7_180.txt';
  filename4='D:\Works and Studies\Photonic\Lumerical\models\mopc16\tr_results7_180.txt';
 %filename2='D:\Works and Studies\Photonic\Lumerical\models\mopc16\Trbb7.txt';

 %filename4='D:\Works and Studies\Photonic\Lumerical\models\mopc16\Trbb4.txt';
 % filename5='D:\Works and Studies\Photonic\Lumerical\models\mopc16\Ex.txt';
   % filename6='D:\Works and Studies\Photonic\Lumerical\models\mopc16\Ez.txt';
    
   % dataEx=csvread(filename5,3,0);
   % dataEz=csvread(filename6,3,0);
   
   filename8= 'D:\Works and Studies\Photonic\Lumerical\models\mopc16\resultsHF7.180.txt';
  data8=csvread(filename8,3,0);
     
data=csvread(filename1,3,0);


ndata=size(data,1);

data2=csvread(filename2,3,0);
data3=csvread(filename3,3,0);
data4=csvread(filename4,3,0);

lams=data(:,1);

vals=data(:,2);

fr=data2(:,2);

a=1.17e-6;
%a=1e-6;
Fn=a./lams;

wn1=Fn(1);
wn2=Fn(ndata);

colT='-b';
colR='-r';

Fn(1);
fr(1);
 figure(1)
plot(Fn,vals,'-.k','LineWidth',2);
 %plot(Fn,fr,colT,'LineWidth',2);
 %plot(data3(:,1),data3(:,2),colT,'LineWidth',2);


 %axis([wn1,wn2,0,1]);
 
 hold on
  %  hold onplot(Fn,data4(:,2),'-g','LineWidth',2);

 plot(data8(:,1),data8(:,3),colR,'LineWidth',2);
 % plot(Fn,fr/90,colR,'LineWidth',2);
  hold on
 % plot(data4(:,1),data4(:,2),'-g','LineWidth',2);
 axis([wn1,wn2,0,1]);
 % plot(Fn,fr,'-k','LineWidth',1);
 
 figure(2)
 plot(Fn,fr,'-.k','LineWidth',2);

 hold on
  %  hold onplot(Fn,data4(:,2),'-g','LineWidth',2);

 plot(data8(:,1),data8(:,2),'-b','LineWidth',2);
 % plot(Fn,fr/90,colR,'LineWidth',2);
   axis([wn1,wn2,0,90]);

  hold on
end