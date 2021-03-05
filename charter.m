function charter
  
  ndata=6;
%[filename1,filepath1]=uigetfile('*.txt', 'Selectinput file')
% D:\Works and Studies\Photonic\Lumerical\models\faradary_isolator
 %filename1='D:\Works and Studies\Photonic\Lumerical\models\mopc8\trslab3.txt';
 %filename2='D:\Works and Studies\Photonic\Lumerical\models\mopc8\frslab3.txt';

 
   filename1= 'D:\personal\abbas\aa-paper\final results\FDTD_tr_gama0.txt';
      filename2= 'D:\personal\abbas\aa-paper\final results\FDTD_fr_216.txt';
%filename1= 'D:\personal\abbas\aa-paper\final results\1D\FDTD_tr.txt';
   %   filename2= 'D:\personal\abbas\aa-paper\final results\1D\FDTD_fr.txt';
      filename3= 'D:\personal\abbas\aa-paper\final results\results-7-180s-216.txt';
     % filename3= 'D:\personal\abbas\aa-paper\final results\1D\results800.txt';

  data3=csvread(filename3,3,0);
     
data=csvread(filename1,3,0);


ndata=size(data,1);

data2=csvread(filename2,3,0);


lams=data(:,1);

vals=data(:,2);

fr=data2(:,2);

a=1.17e-6;
%a=1e-6;
Fn=a./lams;

wn1=Fn(1);
wn2=Fn(ndata);

colT='-r';
colR='-b';

Fn(1);
fr(1);
 figure(1)
plot(Fn,vals,'-.k','LineWidth',2);
 %plot(Fn,fr,colT,'LineWidth',2);
 %plot(data3(:,1),data3(:,2),colT,'LineWidth',2);


 %axis([wn1,wn2,0,1]);
 
 hold on
  %  hold onplot(Fn,data4(:,2),'-g','LineWidth',2);

 plot(data3(:,1),data3(:,3),colT,'LineWidth',2);
 % plot(Fn,fr/90,colR,'LineWidth',2);
  hold on
 % plot(data4(:,1),data4(:,2),'-g','LineWidth',2);
 axis([wn1,wn2,0,1]);
 % plot(Fn,fr,'-k','LineWidth',1);
 
 figure(2)
 plot(Fn,fr,'-.k','LineWidth',2);

 hold on
  %  hold onplot(Fn,data4(:,2),'-g','LineWidth',2);

 plot(data3(:,1),data3(:,2),colR,'LineWidth',2);
 % plot(Fn,fr/90,colR,'LineWidth',2);
   axis([wn1,wn2,0,90]);

  hold on
  
  
fid = fopen('FDTD_results.txt','wt');  % Note the 'wt' for writing in text mode

fprintf(fid,' FDTD results\n');  
  
fprintf(fid,'******\n');
fprintf(fid,'[Fn *  Rotation * Transmitance ]\n');  
  for p=1:ndata
  
  result(p,1)= Fn(p);
  result(p,2)= fr(p);
  result(p,3)= vals(p);

  fprintf(fid,'%f, %f, %f\n',result(p,1),result(p,2),result(p,3));


  end
    fclose(fid);
end