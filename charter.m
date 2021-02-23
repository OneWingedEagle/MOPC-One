function charter
  
  ndata=6;
%[filename1,filepath1]=uigetfile('*.txt', 'Selectinput file')
% D:\Works and Studies\Photonic\Lumerical\models\faradary_isolator
 %filename1='D:\Works and Studies\Photonic\Lumerical\models\mopc8\trslab3.txt';
 %filename2='D:\Works and Studies\Photonic\Lumerical\models\mopc8\frslab3.txt';

   filename1= 'D:\Works and Studies\Photonic\Lumerical\models\faradary_isolator\tr_mansuri.txt';
      filename2= 'D:\Works and Studies\Photonic\Lumerical\models\faradary_isolator\fr_mansuri.txt';
     filename8= 'D:/Works and Studies/Photonic/Faraday\results.txt';
       filename8= 'D:/Works and Studies/Photonic/Faraday\results_mansuri.txt';

  data8=csvread(filename8,3,0);
     
data=csvread(filename1,3,0);


ndata=size(data,1);

data2=csvread(filename2,3,0);


lams=data(:,1);

vals=data(:,2);

fr=data2(:,2);

a=1.17e-6;
a=1e-6;
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

 plot(data8(:,1),data8(:,3),'-ob','LineWidth',2);
 % plot(Fn,fr/90,colR,'LineWidth',2);
  hold on
 % plot(data4(:,1),data4(:,2),'-g','LineWidth',2);
 axis([wn1,wn2,0,1]);
 % plot(Fn,fr,'-k','LineWidth',1);
 
 figure(2)
 plot(Fn,fr,'-.k','LineWidth',2);

 hold on
  %  hold onplot(Fn,data4(:,2),'-g','LineWidth',2);

 plot(data8(:,1),data8(:,2),'-og','LineWidth',2);
 % plot(Fn,fr/90,colR,'LineWidth',2);
   axis([wn1,wn2,0,90]);

  hold on
end