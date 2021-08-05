function charter

%filename= 'D:\Works and Studies\Photonic\Lumerical\models\faradary_isolator\results-7-180s-108.txt';
[filename,filepath]=uigetfile('*.txt', 'Select results file')
file=strcat(filepath,strcat('\',filename));
data0=csvread(file,3,0);
data0

[filename1,filepath1]=uigetfile('*.txt', 'Select results file')
file1=strcat(filepath1,strcat('\',filename1));
data1=csvread(file1,3,0);
     

data=data0;

ndata=size(data(:,1),1);

reverse=1;
if(reverse==1)
for j=1:ndata
data(j,1)=1./data0(ndata-j+1,1);
data(j,2)=data0(ndata-j+1,2);
data(j,3)=data0(ndata-j+1,3);
end

data

fid = fopen('results_lam.txt','wt');  % Note the 'wt' for writing in text mode

fprintf(fid,'results vs lambda\n');  
  
fprintf(fid,'******\n');
fprintf(fid,'[Fn *  Rotation * Transmitance ]\n');  
  for p=1:ndata
  

  fprintf(fid,'%f, %f, %f\n',data(p,1),data(p,2),data(p,3));


  end
    fclose(fid);
end

Fn=data(:,1);


wn1=Fn(1);
wn2=Fn(ndata);

colT='-k';
colR='-r';


 figure(1)

 plot(Fn(:,1),data(:,3),colR,'LineWidth',2);
  hold on
 plot(data1(:,1),data1(:,3),'-.k','LineWidth',2);
 axis([wn1,wn2,0,1]);
 
 figure(2)
 plot(Fn(:,1),data(:,2),colT,'LineWidth',2);
   hold on
 plot(data1(:,1),data1(:,2),'-.k','LineWidth',2);

 tmax=1.1*max(abs(data(:,2)));
   axis([wn1,wn2,0,90]);

  hold on
end