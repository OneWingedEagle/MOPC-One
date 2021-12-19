function faraday_faster1D

clear all

colors = {'-or', '-ob', '-oc', '-ok','-*r', '-*b',...
 '-*c', '-*k','-xr', '-xb', '-xc', '-xk'};
%C = {'k','b','r','g','y'}; % Cell array of colros.


[filename1,filepath1]=uigetfile('*.txt', 'Selectinput file')
 cd(filepath1)
 fid= fopen(filename1)

line=getNewDataLine(fid);
numbs = str2num(line);
ddy1=numbs(1);
ddy2=numbs(2);
a2=ddy1+ddy2;

d1=0;
if(length(numbs)>2)
d1=numbs(3);
end


d2=d1;

line=getNewDataLine(fid);
numbs = str2num(line);
eps1=numbs(1);
eps3=numbs(2);

line=getNewDataLine(fid);
numbs = str2num(line);
ndatab=length(numbs);
gamab=0;
epsbx=numbs(1);
if(ndatab>1)
gamab=numbs(2);
end

line=getNewDataLine(fid);
numbs = str2num(line);
ndataa=length(numbs);
epsax=numbs(1);
gamaa=0;
if(ndataa>1)
gamaa=numbs(2);
end

if(ndatab>1)

epsb=[epsbx -1i*gamab;1i*gamab epsbx ];

epsa=[epsax -1i*gamaa;1i*gamaa epsax ];
else
epsb=[epsbx];

epsa=[epsax];
end

if(ndataa>1)

disp('epsilon tensor of background');
epsb
disp('epsilon tensor of holes');
epsa
end

line=getNewDataLine(fid);
numbs = str2num(line);
Na=numbs(1);
global ndef;
global defstart;
ndef=0;
defstart=0;
if(length(numbs)>1)
ndef=numbs(2);
defstart=floor(Na/2);
end
if(length(numbs)>2)
defstart=numbs(3);
end
line=getNewDataLine(fid);

numbs = str2num(line);
transmit=numbs(1);
rotation=numbs(2);

if(transmit>0)
colT=colors{transmit};
else
colT=colors{1};
end

if(rotation>0)
colR=colors{rotation+1};
else
colR=colors{2};
end

line=getNewDataLine(fid);
numbs = str2num(line);

Fn1=numbs(1);
Fn2=numbs(2);
ndiv=numbs(3);

Fn0=Fn1;

line=getNewDataLine(fid);
numbs = str2num(line);
nGy=numbs(1)

line=getNewDataLine(fid);
numbs = str2num(line);
plotFT=numbs(1);
plotWave=numbs(2);

%=================

t1=cputime;
if(ndiv==0)
dFn=0;
else
dFn=(Fn2-Fn1)/ndiv;
end

cf=(plotWave==0);

Fr_hom=zeros(1*cf*ndiv+1,1);


for p=1:1*cf*ndiv+1
      
     p
     
    Fn(p)=Fn1+dFn*(p-1);
    
%   if(ndivth>0)
 %   wn=Fn0;
 %   theta=Fn(p);
 % else
     wn=Fn(p);
 %   end

    k1=sqrt(eps1)*2*pi*wn/a2;
     Fr_hom(p)=-180*gamab/sqrt(epsbx)*Na*wn;

     pp=p;

    [Ts Rs,Fr, Tsx]=calculteFaraday(epsa,epsb,eps1,eps3,a2,ddy1,d1,d2,Na,nGy,k1,pp,plotFT,plotWave,colR);
           
      if(real(Ts)>1) 
      %  Ts=1;
      end
      
    FR(p)=real(Fr);
 
    Tt(p)=real(Ts);
    Ttx(p)=real(Tsx);
     Fn_fr_tr_trx=[Fn(p)  FR(p)  Tt(p) Ttx(p)]
    
    uu=Ts+Rs;

end

t2=cputime;

comptation_time=t2-t1;
comptation_time
	
result=zeros(1*cf*ndiv+1,4);


fidx = fopen('analyt.txt','wt');  % Note the 'wt' for writing in text mode
fid = fopen('results.txt','wt');  % Note the 'wt' for writing in text mode

fprintf(fid,'[nGx *  nGy]\n');  
  
fprintf(fid,'%d\t%d\n',0,nGy);
fprintf(fid,'[Fn *  Rotation * Transmitance ]\n');  


  for p=1:1*cf*ndiv+1
  
  result(p,1)= Fn(p);
  result(p,2)= FR(p);
  result(p,3)= Tt(p);
  result(p,4)= Ttx(p);


  fprintf(fid,'%f, %f, %f, %f\n',result(p,1),result(p,2),result(p,3),result(p,4));
  fprintf(fidx,'%f, %f\n',result(p,1),Fr_hom(p));


  end
    fclose(fid);
     fclose(fidx);
     
  disp('Results:');
  disp('[Fn *  Rotation * Transmitance ]');	
  disp(result);
  
if(rotation &&length(FR)>1)
                figure(1)
             plot(Fn,FR,colR);
             
            %  axis([Fn1,Fn2,-90,0]);

             axis([Fn1,Fn2,-90,90]);
             hold on
             
         %  plot(Fn,Fr_hom,'+k');
         %    hold on

end
            
if(transmit &&length(Tt)>1)
              figure(2)
               plot(Fn,Tt,colT);
                 axis([Fn1,Fn2,min(Tt),max(Tt)]);
                %axis([Fn1,Fn2,0,1]);
                 hold on
            
end

end



function [Ts Rs Fr Tsx]=calculteFaraday(epsa,epsb,eps1,eps3,a2,...
ddy1,d1,d2,Na,nGy,k1,p,plotFT,plotWave,colorAng)

dmm=size(epsb,2);
nk=1;
if(dmm>1)
nk=2;
end

single_percision=1; % using single or double percison numbers

d=d1+d2;

triang=0;

global MM;

global Kapa;
global bE;

w2c2=(k1)^2;

L=Na*a2+d;

Ly1=2*nGy+1;


numG=nGy;



E0=[0 1]';
if(nk==1)
E0=[1]';
end;


k3=sqrt(eps3/eps1*w2c2);

by=pi/L;

pph=2/pi;

dimx=nk*(numG+2);

if(p==1)
    disp('Computing Fourier series ..');
    Kapa=zeros(4*nGy+1,2)+1i*zeros(4*nGy+1,2);
 
   KapaRectangle(nGy,epsa,epsb,L,ddy1,Na,a2,d1);  

 
   if(plotFT)    
        
 	      ndy=40*Na;
        
        y=linspace(-L,L,ndy);
        Kapar=zeros(ndy);       

            for iy=1:ndy

                y1=y(iy);
                tt=0;
                    for m=-nGy:nGy
                        Gm=by*m;

                        TT=Kapa(m+2*nGy+1,plotFT);
                        
                        tt=tt+TT*exp(1i*(Gm*y1));
                    end

                Kapar(iy)=real(tt);
            end
   

        figure(4);

        plot(y,Kapar);

    end

    disp('Computing matrix, step 1...');
    
    dimx=nk*(numG+2);
   
   if(single_percision)
    MM=  single(zeros(dimx,dimx));
  else
     MM=  zeros(dimx,dimx);
    end
        
    bE=zeros(dimx,1);
    
    countG=0;
           
    for Gx=0:0
        for Gy=1:nGy
            
            countG=countG+1;
            countG1=0;
            
            for Gx1=0:0
                                           
                CT=zeros(nk,nk);
                CR=zeros(nk,nk);
                CE=zeros(nk,nk);
                
                for Gy1=1:nGy
                    kym=Gy1*by;
                    dGyp=Ly1+(Gy-Gy1);
                    sGyp=Ly1+Gy+Gy1;
                    countG1=countG1+1;
                         
                      if(nk==1) 
                    
                     KapTenss=[Kapa(sGyp,1)];
              
                    KapTensd=[Kapa(dGyp,1)];                 
                    else
                    
                    Kapas(1:2)=Kapa(sGyp,1:2);                                       
                    Kapad(1:2)=Kapa(dGyp,1:2);

                    KapTenss=[Kapas(1) 1i*Kapas(2); 
                        -1i*Kapas(2)  Kapas(1)];

                    
                    KapTensd=[Kapad(1) 1i*Kapad(2);
                        -1i*Kapad(2) Kapad(1)];
                    end
                    


                    CA=(KapTensd-KapTenss)*(Gy1*by)^2;
                   

                    
                    r1=(countG-1)*nk;
                    c1=(countG1-1)*nk;
                    for j=1:nk
                        for k=1:nk
                            MM(r1+j,c1+k)= -CA(j,k);
                            
                        end
                    end
                    
                    
                end
                          

                    for k=1:nk                          
                            MM((numG)*nk+k,(Gy-1)*nk+k)=pi*Gy;                         
                    end
                    
                    if(mod(Gy,2)==0)
                        for k=1:nk
                                MM((1+numG)*nk+k,(Gy-1)*nk+k)=pi*Gy;
                        end
                    else
                        for k=1:nk
                               MM((1+numG)*nk+k,(Gy-1)*nk+k)=-pi*Gy;
                        end
                    end

                
            end
            
            
        end
        
        
    end
    
end

MM;

bN=zeros(nk*(numG+2),1);

disp('Computing matrix, step 2...');

if(single_percision)
    NN=single(zeros(nk*(numG+2)));
else 
     NN=zeros(nk*(numG+2));
end  

    for Gx=1:size(NN,1)
        for Gy=1:size(NN,2)
          if(single_percision)
            NN(Gx,Gy)=single(0+1i*0);
          else
            NN(Gx,Gy)=0+1i*0;
           end
        end
    end


countG=0;

   
    for Gy=1:nGy
        
        countG=countG+1;
        for k=1:nk
                NN((countG-1)*nk+k,(countG-1)*nk+k)= w2c2;
        end
        
        if(mod(Gy,2)==0)
            for k=1:nk
                     NN((countG-1)*nk+k,(numG)*nk+k)=-pph/Gy*w2c2;
            end
        else
            for k=1:nk
                    NN((countG-1)*nk+k,(numG)*nk+k)=pph/Gy*w2c2;

            end
        end
        
        for k=1:nk

                NN((countG-1)*nk+k,(numG+1)*nk+k)=pph/Gy*w2c2;   
                bN((countG-1)*nk+k)=-pph/Gy*w2c2*E0(k);
                
        end

end





    krny=-k1;
    
     ktny=k3;
    
    for k=1:nk
        
                NN((numG)*nk+k,(numG)*nk+k)=1;
                
                
                NN((numG)*nk+k,(numG+1)*nk+k)=-(1i*L*krny+1);
                
                
                NN((numG+1)*nk+k,(numG)*nk+k)=-(1i*L*ktny-1);
                
                NN((numG+1)*nk+k,(numG+1)*nk+k)=-1;
        
   
             bN((numG)*nk+k)= (1i*L*k1+1)*E0(k);
             if(k==1 && nk~=1)
              bN((numG+1)*nk+k)=-E0(k);
              else
              bN((numG+1)*nk+k)=E0(k);
            end

    end

NN;
bN=bN+bE;

disp('solving matrix...');
bN;

  NN=NN+MM;   
NN;
gpu=0;

if(gpu>0)
 gpux=gpuArray(NN)\gpuArray(bN);
 x=gather(gpux);
else

 x=linsolve(NN,bN);
 end


Anm=zeros(nGy,nk);
Tn=zeros(nk,1);
Tn2=zeros(1,1);
Tn2x=zeros(1,1);

Rn=zeros(nk,1);
Rn2=zeros(1,1);

ix=0;

    for Gy=1:nGy
        for k=1:nk
            ix=ix+1;
            Anm(Gy,k)=x(ix);
        end
    end



kp=numG*nk;
kpp=(numG+1)*nk;


    for j=1:nk
      
        Tn(j)=x(kp+j);
        
        Rn(j)=x(kpp+j);
        
        if(j==1 && nk>1)
         Tn2x(1)= Tn2x(1)+Tn(j)*conj(Tn(j));
        end
        Tn2(1)= Tn2(1)+Tn(j)*conj(Tn(j));

        Rn2(1)= Rn2(1)+Rn(j)*conj(Rn(j));
        
    end

theta=0;
Ts=0;
Tsx=0;

Rs=0;

Ts0=0;
Tsx0=0;

 krny=-k1;
    
    
ktny=k3;
        
    %     Ts=Ts+Tn2(k);%/cosd(theta);
    Ts=Ts+ktny/k3*sqrt(eps3/eps1)*Tn2(1)/cosd(theta);
        %Ts=Ts+ktny/k3*sqrt(eps3/eps1)*abs(Tn(k,3))^2/cos(thetad);
    Tsx=Tsx+ktny/k3*sqrt(eps3/eps1)*Tn2x(1)/cosd(theta);

        Ts0=ktny/k3*sqrt(eps3/eps1)*Tn2(1)/cosd(theta);
        Tsx0=ktny/k3*sqrt(eps3/eps1)*Tn2x(1)/cosd(theta);
    
    Rs=Rs+abs(krny)/k1*Rn2(1)/cosd(theta);





Fr=0;

if(nk>1)

   
% TM_fract=abs(Tsx0/Ts0);
%ct=sqrt(TM_fract);

%angle=asin(ct);
Tn;
angle=atand(abs(Tn(1))/abs(Tn(nk)));

if(real(Tn(1))/real(Tn(nk))<0) 
%angle=-angle;
end;

Fr=angle;


if(plotWave)
    
    figure(7)
    
    plot(y,angs,colorAng);
    hold on
    
    %angsHomog=-0.5*y*k1*imag(epsb(1,3))/sqrt(epsb(1,1))*180/pi;
    %plot(y,angsHomog,'-k');
    
end

end

end


function KapaRectangle(nGy,epsa,epsb,L,ddy1,Na,a2,d1)
global Kapa;

disp('strips');
nk=2;
if(length(epsa)==1)
nk=1;
end

invepsa1=inv(epsa);
invepsb1=inv(epsb);

if(nk==1)

invepsa=[invepsa1(1,1)];
invepsb=[invepsb1(1,1)];

else
invepsa=[invepsa1(1,1)  imag(invepsa1(1,2))];
invepsb=[invepsb1(1,1)  imag(invepsb1(1,2))];

end

cell_size=a2*Na*2;
rect_size=a2;

filly=ddy1/a2;

by=pi/L;


KapaUnit=zeros(4*nGy+1,2)+1i*zeros(4*nGy+1,2);


for dGy=-2*nGy:2*nGy
        
        tty=dGy*by*ddy1/2;
     
        dGyp=dGy+1+2*nGy;
        
  for k=1:nk
        
         if(dGy==0)
       four_coefy= (invepsb(k)+filly*(invepsa(k)-invepsb(k)))/(2*Na);
        else
          four_coefy=(invepsa(k)-invepsb(k))*ddy1/(a2)*sin(tty)/(tty)/(2*Na);
       end
       
        KapaUnit(dGyp,k)=four_coefy;


  end

end

global ndef;
global defstart;

by=pi/L;


KapaDefect=zeros(4*nGy+1,2)+1i*zeros(4*nGy+1,2);
 


if(ndef>0)

     for dGy=-2*nGy:2*nGy
        
        dGyp=dGy+1+2*nGy;

        tt=dGy*by*a2/2;
          
          for k=1:nk
                if(dGy==0)
                  KapaDefect(dGyp,k)=invepsb(k)*a2/(2*L);
                else
                  KapaDefect(dGyp,k)=invepsb(k)*a2/(2*L)*sin(tt)/(tt);
             end
         
         end
    
        
    end
 
end
 

Kapa=zeros(4*nGy+1,2)+1i*zeros(4*nGy+1,2);

ndef
defstart

defcount=0;
isDef=zeros(2*Na,1);
for n=0:Na-1 
  np=n+Na+1;
 if(n+1>=defstart && defcount<ndef)
 defcount=defcount+1;
 isDef(np,1)=1;
 isDef(Na-n,1)=1;
 end
end

by=pi/L;

for n=-Na:Na-1 
 np=n+Na+1;

 for dGy=-2*nGy:2*nGy
       
      dGyp=dGy+1+2*nGy;
        twindle=exp(-1i*(n+.5)*by*dGy*a2);
  if(isDef(np,1))

    for k=1:nk
        Kapa(dGyp,k)=  Kapa(dGyp,k)+KapaDefect(dGyp,k)*twindle;
     end
   else
     for k=1:nk
         Kapa(dGyp,k)=Kapa(dGyp,k)+KapaUnit(dGyp,k)*twindle;
     end
   
   end    
    end
  end
 

end

function line=getNewDataLine(fid)

TF=1;
k=0;
while (k<100 && TF==1)
line=fgets(fid);
TF = strncmpi(line,'/',1);
k=k+1;

end


end


function [h]=arrow3d(x,y,z,head_frac,radii,radii2,colr)

%
% Written by Moshe Lindner , Bar-Ilan University, Israel.
% July 2010 (C)

if nargin==5
    radii2=radii*2;
    colr='blue';
elseif nargin==6
    colr='blue';
end
if size(x,1)==2
    x=x';
    y=y';
    z=z';
end

x(3)=x(2);
x(2)=x(1)+head_frac*(x(3)-x(1));
y(3)=y(2);
y(2)=y(1)+head_frac*(y(3)-y(1));
z(3)=z(2);
z(2)=z(1)+head_frac*(z(3)-z(1));
r=[x(1:2)',y(1:2)',z(1:2)'];

N=10;
dr=diff(r);
dr(end+1,:)=dr(end,:);
origin_shift=(ones(size(r))*(1+max(abs(r(:))))+[dr(:,1) 2*dr(:,2) -dr(:,3)]);
r=r+origin_shift;

normdr=(sqrt((dr(:,1).^2)+(dr(:,2).^2)+(dr(:,3).^2)));
normdr=[normdr,normdr,normdr];
dr=dr./normdr;
Pc=r;
n1=cross(dr,Pc);
normn1=(sqrt((n1(:,1).^2)+(n1(:,2).^2)+(n1(:,3).^2)));
normn1=[normn1,normn1,normn1];
n1=n1./normn1;
P1=n1+Pc;

X1=[];Y1=[];Z1=[];
j=1;
for theta=([0:N])*2*pi./(N);
    R1=Pc+radii*cos(theta).*(P1-Pc) +...
    radii*sin(theta).*cross(dr,(P1-Pc)) -origin_shift;
    X1(2:3,j)=R1(:,1);
    Y1(2:3,j)=R1(:,2);
    Z1(2:3,j)=R1(:,3);
    j=j+1;
end

r=[x(2:3)',y(2:3)',z(2:3)'];

dr=diff(r);
dr(end+1,:)=dr(end,:);
origin_shift=(ones(size(r))*(1+max(abs(r(:))))+[dr(:,1) 2*dr(:,2) -dr(:,3)]);
r=r+origin_shift;

normdr=(sqrt((dr(:,1).^2)+(dr(:,2).^2)+(dr(:,3).^2)));
normdr=[normdr,normdr,normdr];
dr=dr./normdr;
Pc=r;
n1=cross(dr,Pc);
normn1=(sqrt((n1(:,1).^2)+(n1(:,2).^2)+(n1(:,3).^2)));
normn1=[normn1,normn1,normn1];
n1=n1./normn1;
P1=n1+Pc;

j=1;
for theta=([0:N])*2*pi./(N);
    R1=Pc+radii2*cos(theta).*(P1-Pc) +...
    radii2*sin(theta).*cross(dr,(P1-Pc)) -origin_shift;
    X1(4:5,j)=R1(:,1);
    Y1(4:5,j)=R1(:,2);
    Z1(4:5,j)=R1(:,3);
    j=j+1;
end

X1(1,:)=X1(1,:)*0 + x(1);
Y1(1,:)=Y1(1,:)*0 + y(1);
Z1(1,:)=Z1(1,:)*0 + z(1);
X1(5,:)=X1(5,:)*0 + x(3);
Y1(5,:)=Y1(5,:)*0 + y(3);
Z1(5,:)=Z1(5,:)*0 + z(3);

h=surf(X1,Y1,Z1,'facecolor',colr,'edgecolor','none');

light('Color','r')

end
