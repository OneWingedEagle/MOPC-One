function convert_to_lambda

%filename= 'D:\Works and Studies\Photonic\Lumerical\models\faradary_isolator\results-7-180s-108.txt';
%[filename,filepath]=uigetfile('*.txt', 'Select results file')
%file=strcat(filepath,strcat('\',filename));
file ='D:\personal\abbas\aa-paper\final results\1D\results-slab_4000.txt';
data=csvread(file,3,0);
     


ndata=size(data,1);



colT='-g';
colR='-g';

output=zeros(ndata,4);
a=1;
for i=1:ndata

 output(i,1)=a/data(ndata+1-i,1); 
 output(i,2)=data(ndata+1-i,2); 
  output(i,3)=sqrt(data(ndata+1-i,3)-data(ndata+1-i,4)); 
   output(i,4)=sqrt(data(ndata+1-i,4)); 

end
 output(:,2);
fid = fopen('results_lam.txt','wt');  % Note the 'wt' for writing in text mode

fprintf(fid,'******* \n');  
  
fprintf(fid,'******\n');
fprintf(fid,'[Fn *  Rotation * Transmitance ]\n');  


  for p=1:ndata

  fprintf(fid,'%f, %f, %f, %f\n',output(p,1),output(p,2),output(p,3),output(p,4));


  end
    fclose(fid);


lam=output(:,1);

wn1=lam(1);
wn2=lam(ndata);

 figure(1)

 plot(lam(:,1),output(:,2),colR,'LineWidth',2);
  hold on

 
 tmax=1.1*max(abs(output(:,2)));
   axis([wn1,wn2,-tmax,tmax]);
 
 figure(2)
 plot(lam(:,1),output(:,3),colT,'LineWidth',2);

 axis([wn1,wn2,.1,.85]);

  hold on
   plot(lam(:,1),output(:,4),colT,'LineWidth',2);
 axis([wn1,wn2,.1,.85]);
end