function faradayX

clear all

colors = {'-or', '-ob', '-oc', '-ok','-*r', '-*b', '-*c', '-*k','-xr', '-xb', '-xc', '-xk'};
%C = {'k','b','r','g','y'}; % Cell array of colros.



[filename1,filepath1]=uigetfile('*.txt', 'Selectinput file')
 cd(filepath1)
 fid= fopen(filename1)
%fid = fopen('input.txt','rt');

line=getNewDataLine(fid);
numbs = str2num(line);
geometry=numbs(1);
rec=0;
if(geometry==1 && length(numbs)>1)
rec=numbs(2);
end
rec
line=getNewDataLine(fid);
numbs = str2num(line);
Rx=numbs(1);
Ry=numbs(2);

line=getNewDataLine(fid);
numbs = str2num(line);
fi=numbs(1);

line=getNewDataLine(fid);
numbs = str2num(line);
theta=numbs(1);
global inc_mode;
inc_mode=0;
if(length(numbs)>1)
inc_mode=numbs(2);
end

line=getNewDataLine(fid);
numbs = str2num(line);
ax=numbs(1);
ay=numbs(2);

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
ndata=length(numbs)

epsbx=numbs(1);
if(ndata>1)
epsby=numbs(2);
epsbz=numbs(3);
gamab=numbs(4);
end

line=getNewDataLine(fid);
numbs = str2num(line);
epsax=numbs(1);
if(ndata>1)
epsay=numbs(2);
epsaz=numbs(3);
gamaa=numbs(4);
end

if(ndata>1)

epsb=[epsbx 0 -1i*gamab;0 epsby 0;1i*gamab 0 epsbz ];

epsa=[epsax 0 -1i*gamaa;0 epsay 0;1i*gamaa 0 epsaz ];
else
epsb=[epsbx];

epsa=[epsax];
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

wn1=numbs(1);
wn2=numbs(2);
ndiv=numbs(3);

line=getNewDataLine(fid);
numbs = str2num(line);
nGx=numbs(1)
nGy=numbs(2)

line=getNewDataLine(fid);
numbs = str2num(line);
plotFT=numbs(1);
plotWave=numbs(2);



%=================

t1=cputime;

dwn=(wn2-wn1)/ndiv;

cf=(plotWave==0);

Fr_hom=zeros(1*cf*ndiv+1,1);

for p=1:1*cf*ndiv+1
    
    p
    

    
    Fn(p)=wn1+dwn*(p-1);
    k1=sqrt(eps1)*2*pi*Fn(p)/ay;
     Fr_hom(p)=-180*gamab/sqrt(epsbx)*Na*Fn(p);
    
    [Ts Rs,Fr]=calculteFaraday(geometry,epsa,epsb,eps1,eps3,ax,ay,Rx,Ry,d1,d2,Na,nGx,nGy,k1,p,plotFT,plotWave,colR,theta,fi,rec);
           
      if(real(Ts)>1) 
        Ts=1;
      end
      
    Tr(p)=real(Fr);
    Tt(p)=real(Ts);
    
    uu=Ts+Rs;

end


t2=cputime;

comptation_time=t2-t1;
comptation_time;
        
        Tr';
        Tt';
	
result=zeros(1*cf*ndiv+1,3);


fidx = fopen('analyt.txt','wt');  % Note the 'wt' for writing in text mode
fid = fopen('results.txt','wt');  % Note the 'wt' for writing in text mode

fprintf(fid,'[nGx *  nGy]\n');  
  
fprintf(fid,'%d\t%d\n',nGx,nGy);
fprintf(fid,'[wn *  Rotation * Transmitance ]\n');  


  for p=1:1*cf*ndiv+1
  
  result(p,1)= Fn(p);
  result(p,2)= Tr(p);
  result(p,3)= Tt(p);

  fprintf(fid,'%f\t%f\t%f\n',result(p,1),result(p,2),result(p,3));
  fprintf(fidx,'%f\t%f\n',result(p,1),Fr_hom(p));


  end
    fclose(fid);
     fclose(fidx);
     
  disp('Results:');
  disp('[wn *  Rotation * Transmitance ]');	
  disp(result);
  
if(rotation &&length(Tr)>1)
                figure(1)
             plot(Fn,Tr,colR);
             
 
             axis([wn1,wn2,-90,90]);
             hold on
             
        %   plot(Fn,Fr_hom,'+k');
        %     hold on

end
            
if(transmit &&length(Tt)>1)
              figure(2)
               plot(Fn,Tt,colT);
                 axis([wn1,wn2,min(Tt),max(Tt)]);
                %axis([wn1,wn2,0,1]);
                 hold on
            
end

end



function [Ts Rs Fr]=calculteFaraday(geometry,epsa,epsb,eps1,eps3,a1,a2,...
Rx,Ry,d1,d2,Na,nGx,nGy,k1,p,plotFT,plotWave,colorAng,theta,fi,rec)

nk=size(epsb,2);


d=d1+d2;

triang=0;

global MM;

global Kapa;
global bE;

w2c2=(k1)^2;

L=Na*a2+d;

Lx1=2*nGx+1;
Lx2=3*nGx+2;

Ly1=2*nGy+1;


numG=(Lx1)*nGy;

dd=Lx1;


E0=[0 0 1]';


kx=k1*sind(theta);
k1y=k1*cosd(theta);
k3=sqrt(eps3/eps1*w2c2);
bx=2*pi/a1;

by=pi/L;

pph=2/pi;

    dimx=nk*(numG+2*(Lx1));



if(p==1)
    disp('Computing Fourier series ..');
    Kapa=zeros(4*nGx+1,4*nGy+1,4)+1i*zeros(4*nGx+1,4*nGy+1,4);
    if(triang)
        FillKapaTriang(nGx,nGy,epsa,epsb,L,R,Na,a1,a2);  %triangular
        %lattivce not implemented properly.
    elseif (geometry==0 && Rx==Ry)
        FillKapaCylinderAntiSymDef(nGx,nGy,epsa,epsb,L,Rx,Na,a1,a2,d1); 

    else
    if(geometry==1 && rec==1)
    FillKapaRectangleAntiSym(nGx,nGy,epsa,epsb,L,Rx,Ry,Na,a1,a2,d1,fi);

  else

       FillKapaAntiSymNum(geometry,nGx,nGy,epsa,epsb,L,Rx,Ry,Na,a1,a2,d1,fi);
 end
 end


    
    if(plotFT)
        
        
        ndx=30;
 	      ndy=40*Na;
        
        x=linspace(-a1/2,a1/2,ndx);
        y=linspace(-L,L,ndy);
        
        [xx yy]=ndgrid(x,y);
        
        zz=zeros(size(xx));
        
        for ix=1:ndx
            x1=x(ix);
            for iy=1:ndy

                y1=y(iy);
                tt=0;
                for n=-nGx:nGx
                    for m=-nGy:nGy
                        Gn=bx*n;
                        Gm=by*m;

                        TT=Kapa(n+2*nGx+1,m+2*nGy+1,1);
                        
                        tt=tt+TT*exp(1i*(Gn*x1+Gm*y1));
                    end
                end
                
                zz(ix,iy)=real(tt);
            end
        end
        
        
    %    X=abs(Kapa(:,:,1));

        
        figure(4);

        surf(xx,yy,zz);

        %  surf(X);
        axis equal;
        set(gca,'DataAspectRatio',[1 1 .05]);
%        return;
    end


    disp('Computing matrix, step 1...');
    
    dimx=nk*(numG+2*(Lx1));
   
    MM=  single(zeros(dimx,dimx));
    

        
    bE=zeros(nk*(numG+2*(Lx1)),1);
    
    countG=0;
    
    
    for Gx=-nGx:nGx
        for Gy=1:nGy
            
            countG=countG+1;
            countG1=0;
            
            for Gx1=-nGx:nGx
                
                dGxp=Lx1+Gx-Gx1;
                Gx1p=nGx+1+Gx1;
                Gx1pp=Lx2+Gx1;


                
                kxn1=kx+Gx1*bx;
               % kxn1=0;
                
                kxn2=kxn1^2;
             %    kxn2=1;
                
                tempT=zeros(nk,nk);
                tempR=zeros(nk,nk);
                tempE=zeros(nk,nk);
                
                for Gy1=1:nGy
                    kym=Gy1*by;
                    dGyp=Ly1+(Gy-Gy1);
                    sGyp=Ly1+Gy+Gy1;
                    countG1=countG1+1;
                         
       
                    
                    Kapas(1:4)=Kapa(dGxp,sGyp,1:4);
      
                                        
                    Kapad(1:4)=Kapa(dGxp,dGyp,1:4);
                    
                    KapTenss=[Kapas(1) 0 1i*Kapas(4); 0 Kapas(2) 0 ;
                        -1i*Kapas(4) 0  Kapas(3)];

                    
                    KapTensd=[Kapad(1) 0 1i*Kapad(4); 0 Kapad(2) 0 ;
                        -1i*Kapad(4) 0  Kapad(3)];

                    
                    Akg=zeros(3,3);
                    
                    
                    Akg(1,1)=kym^2;
                    Akg(1,2)=-kxn1*kym;
                    Akg(2,1)= Akg(1,2);
                    Akg(2,2)=kxn2;
                    
                    Akg(3,3)=Akg(1,1)+ Akg(2,2);

                    AkdKap=(KapTenss*Akg-KapTensd*conj(Akg));
                    
                    Akg1=zeros(3,3);
                    
##                    if(mod(Gy1,2)==0)
##                        Akg1(1,2)=kxn1*pph/Gy1;
##                    else
##                        Akg1(1,2)=-kxn1*pph/Gy1;
##                    end
##                    
##                    Akg1(2,1)=Akg1(1,2);
                    
                    if(mod(Gy1,2)==0)
                        Akg1(2,2)=kxn2*pph/Gy1;
                    else
                        Akg1(2,2)=-kxn2*pph/Gy1;
                    end
                    
                    
                    Akg1(3,3)=Akg1(2,2);
                    
                    
                    dKapTens1=-1*(KapTenss*Akg1-KapTensd*conj(Akg1));
                    
              
                    tempT=tempT+dKapTens1;

                    
                   Akg2=zeros(3,3);
##                    
##                    Akg2(1,2)=-kxn1*pph/Gy1;
##         
##                    Akg2(2,1)=Akg2(1,2);
                    Akg2(2,2)=kxn2*pph/Gy1;
                   Akg2(3,3)=Akg2(2,2);


                   dKapTens2=-1*(KapTenss*Akg2-KapTensd*conj(Akg2));
           
                    
                   tempR=tempR-dKapTens2;
                    
                    if(Gx1==0)
                        
                       tempE=tempE-dKapTens2;
                    end
                    
                    r1=(countG-1)*nk;
                    c1=(countG1-1)*nk;
                    for j=1:nk
                        for k=1:nk
                            MM(r1+j,c1+k)= AkdKap(j,k);
                            
                        end
                    end
                    
                    
                end
                
                
                
                if(Gx1==0)
                    v=tempE*E0;
                    for k=1:nk
                        bE((countG-1)*nk+k)=v(k);
                    end
                    
                end
                
                
                for j=1:nk
                    for k=1:nk
                        
                        MM((countG-1)*nk+j,(numG+Gx1p-1)*nk+k)=tempT(j,k);
                        MM((countG-1)*nk+j,(numG+dd+Gx1p-1)*nk+k)=tempR(j,k);
                        
                    end
                end
                
                
                
                if(Gx==Gx1)
                    for k=1:nk
                        if(k~=2)
                            
                            MM((Gx1p+numG-1)*nk+k,((Gx1p-1)*nGy+Gy-1)*nk+k)=pi*Gy;
                        end
                    end
                    
                    if(mod(Gy,2)==0)
                        for k=1:nk
                            if(k~=2)
                                MM((Gx1pp+numG-1)*nk+k,((Gx1p-1)*nGy+Gy-1)*nk+k)=pi*Gy;
                            end
                        end
                    else
                        for k=1:nk
                            if(k~=2)
                                MM((Gx1pp+numG-1)*nk+k,((Gx1p-1)*nGy+Gy-1)*nk+k)=-pi*Gy;
                            end
                        end
                    end
                    
                end
                
            end
            
            
        end
        
        
    end
    
end


bN=zeros(nk*(numG+2*(Lx1)),1);

disp('Computing matrix, step 2...');

    NN=single(zeros(nk*(numG+2*(Lx1))));

    for Gx=1:size(NN,1)
        for Gy=1:size(NN,2)
            NN(Gx,Gy)=single(1e-10+1i*1e-10);
        end
    end


countG=0;

for Gx=-nGx:nGx
    Gxp=nGx+1+Gx;
    Gxpp=Lx2+Gx;
    
   
    for Gy=1:nGy
        
        countG=countG+1;
        for k=1:nk
                NN((countG-1)*nk+k,(countG-1)*nk+k)= w2c2;
        end
        
        if(mod(Gy,2)==0)
            for k=1:nk
                     NN((countG-1)*nk+k,(numG+Gxp-1)*nk+k)=-pph/Gy*w2c2;
            end
        else
            for k=1:nk
                    NN((countG-1)*nk+k,(numG+Gxp-1)*nk+k)=pph/Gy*w2c2;

            end
        end
        
        for k=1:nk

                NN((countG-1)*nk+k,(numG+Gxpp-1)*nk+k)=pph/Gy*w2c2;
      
            if(Gx==0)
                
                bN((countG-1)*nk+k)=-pph/Gy*w2c2*E0(k);
                
            end
        end
        
    end
end


for Gx=-nGx:nGx
    Gxp=nGx+1+Gx;
    Gxpp=Lx2+Gx;
    
    kxn1=kx+Gx*bx;
    kxn2=kxn1^2;
    
    
    if(k1>=abs(kxn1))
        krny=-sqrt(k1^2-kxn2);
    else
        krny=-1i*sqrt(kxn2-k1^2);
    end
    
    
    if(k3>=abs(kxn1))
        ktny=sqrt(k3^2-kxn2);
    else
        ktny=1i*sqrt(kxn2-k3^2);
    end
    

    
    for k=1:nk
        
       
            if(k~=2)
                NN((numG+Gxp-1)*nk+k,(numG+Gxp-1)*nk+k)=1;
                
                
                NN((numG+Gxp-1)*nk+k,(numG+Gxpp-1)*nk+k)=-(1i*L*krny+1);
                
                
                NN((numG+Gxpp-1)*nk+k,(numG+Gxp-1)*nk+k)=-(1i*L*ktny-1);
                
                NN((numG+Gxpp-1)*nk+k,(numG+Gxpp-1)*nk+k)=-1;
            else
                
                NN((numG+Gxp-1)*nk+k,(numG+Gxp-1)*nk+k)=1;
            
                NN((numG+Gxpp-1)*nk+k,(numG+Gxpp-1)*nk+k)=1;
            end
        
        if(Gx==0)
                      
            if(k~=2)
            bN((numG+Gxp-1)*nk+k)= (1i*L*k1y+1)*E0(k);
            bN((numG+Gxpp-1)*nk+k)=E0(k);
            end
        end
    end
end


bN=bN+bE;

  
disp('solving matrix...');


  NN=NN+MM;   


  x=linsolve(NN,bN);


Anm=zeros(2*nGx+1,nGy,nk);
Tn=zeros(2*nGx+1,nk);
Tn2=zeros(2*nGx+1,1);
Rn=zeros(2*nGx+1,nk);
Rn2=zeros(2*nGx+1,1);

ix=0;
for Gx=-nGx:nGx
    Gxp=nGx+1+Gx;
    for Gy=1:nGy
        for k=1:nk
            ix=ix+1;
            Anm(Gxp,Gy,k)=x(ix);
        end
    end
end


amx2=0;
if(p==-1)
A2=zeros(2*nGx+1,nGy);
uu=zeros(2*nGx+1,nGy);
vv=zeros(2*nGx+1,nGy);
for Gx=-nGx:nGx
    Gxp=nGx+1+Gx;
 for Gy=1:nGy
        for k=1:nk
         A2(Gxp,Gy)= A2(Gxp,Gy)+Anm(Gxp,Gy,k)*conj(Anm(Gxp,Gy,k));

      end
      
      if(amx2<A2(Gxp,Gy))
      amx2=A2(Gxp,Gy);
      end
  end
end

figure(14)
plot(A2(2,:),'-g');
%surf(sqrt(A2));
hold on;

end

kp=numG*nk;
kpp=(numG+2*nGx+1)*nk;

for Gx=-nGx:nGx
    k=nGx+Gx;
   ords(k+1)=k;
  
    for j=1:nk
        Tn(k+1,j)=x(kp+k*nk+j);
        Tn2(k+1)= Tn2(k+1)+Tn(k+1,j)*conj(Tn(k+1,j));
        
        Rn(k+1,j)=x(kpp+k*nk+j);
        Rn2(k+1)= Rn2(k+1)+Rn(k+1,j)*conj(Rn(k+1,j));
        
    end
    

end

Ts=0;
Rs=0;

for Gx=-nGx:nGx
    k=nGx+1+Gx;
    kxn1=kx+Gx*bx;
    kxn2=kxn1^2;
    
    if(k1>=abs(kxn1))
        krny=-sqrt(k1^2-kxn2);
    else
        krny=-1i*sqrt(kxn2-k1^2);
    end
    
    
    if(k3>=abs(kxn1))
        ktny=sqrt(k3^2-kxn2);
    else
        ktny=1i*sqrt(kxn2-k3^2);
    end
        
    Ts=Ts+ktny/k3*sqrt(eps3/eps1)*Tn2(k)/cosd(theta);
        %Ts=Ts+ktny/k3*sqrt(eps3/eps1)*abs(Tn(k,3))^2/cos(thetad);

    
    Rs=Rs+abs(krny)/k1*Rn2(k)/cosd(theta);
end

Fr=0;


nL=max(int32(L*k1/2/pi)*17*2,35);



yy=linspace(0,L,nL);

Nx=1;

if(Nx==1)
    xx=zeros(1,1);
else
    xx=linspace(-a1/2,a1/2,Nx);
end

ff=zeros(Nx,nL,nk);
si=zeros(Nx,nL,nk);

for ix=1:Nx
    
    x=xx(ix);
    
    for k=1:nL
        y=yy(k);
        
        for Gx=-nGx:nGx
            
            kxn=Gx*bx;
            
            Gxp=nGx+1+Gx;
            if(Gx==0)
            del=1;
          else 
            del=0;
            end
            for Gy=1:nGy
                for j=1:nk

		 ky0=Gy*by;

                    si(ix,k,j)=si(ix,k,j)+Anm(Gxp,Gy,j)*sin(ky0*y)*exp(1i*kxn*x);
                end
            end
            for j=1:nk
                ff(ix,k,j)=ff(ix,k,j)+y/L*Tn(Gxp,j)*exp(1i*kxn*x)*exp(1i*sqrt(k3^2-kxn2)*(y-L))+(1-y/L)*(Rn(Gxp,j)+del*E0(j))*exp(1i*kxn*x)*exp(1i*sqrt(k1^2-kxn2)*y);
            end
        end
        
    end
end

E1=ff+si;

E2=real(E1);
E2(:,:,1);

%writeMeshAndField(Nx,nL,1,E2,2,Na);


if(plotWave)
    figure(5)
    
    set(gca,'DataAspectRatio',[1 1 1]);
    axis([-2 2 0 L -2 2]);
    az = 40;
    el = 30;
    %  az=90;
    %  el=0;
    view(az, el);
    
    
    
    hold on
    
    for ix=1:Nx
        x=xx(ix);
        Vx=E2(ix,1,1);
        Vy=E2(ix,1,2);
        Vz=E2(ix,1,3);
        
        
        for k=1:nL
            
            y=yy(k);
            Vxp=Vx;
            Vx=E2(ix,k,1);
            
            Vyp=Vy;
            Vy=E2(ix,k,2);
            
            Vzp=Vz;
            Vz=E2(ix,k,3);
            
            
            color = 'r';
            
            
            if(k>1)
                h1(k)=plot3([Vxp+x Vx+x],[Vyp+yy(k-1) Vy+y],[Vzp Vz],'LineWidth',2,'Color','b');
                h2(k)= plot3([-Vxp+x -Vx+x],[yy(k-1) y],[-Vzp -Vz],'LineWidth',2,'Color','b');
                h3= plot3([0 0],[0 L],[0 0],'LineWidth',2,'Color','k');
                
                
            end
            
            
            arrow1(k) = arrow3d([x x+Vx],[y y+Vy],[0 Vz],.92,.01,.02,color);
            arrow1(k) = arrow3d([x x-Vx],[y y-Vy],[0 -Vz],.92,.01,.02,color);
            
            
            
        end
        
    end;
end


if(nk==3)

tn0=E2(1,1,1)/E2(1,1,3);


ang0=atan(tn0);
nan=0;
tans1=zeros(nL,1);
for k=1:nL
    
    Vx=E2(1,k,1);
    Vz=E2(1,k,3);
    
    magn=sqrt(Vx^2+Vz^2);
    
    if( magn>-.1)
    tans1(k)=Vx/Vz;
        nan=nan+1;
    else
        tans1(k)=1e10;
    end
end

tans=zeros(nan,1);
y=zeros(nan,1);
j=0;
for k=1:nL
    if(tans1(k)~=1e10)
        j=j+1;
        
        tans(j)=tans1(k);
        
        y(j)=yy(k);
    end
end;

angs=zeros(nan,1);
angs(1)=atan(tans(1));

for k=2:nan
    
    tnp=tans(k-1);
    tn=tans(k);
    angs(k)=-atan(tn);

    if(angs(k)>0)
  
%    angs(k)=-180+angs(k);
    
    end

    %  angs(k)=E2(1,k,3);
    % angs(k)=angs(k-1)+atan((tn-tnp)/(1+tn*tnp))*180/pi;
    
    
end



Fr=(angs(nan))*180/pi;

if(Fr>90)
 Fr=Fr-180;
 elseif(Fr<-90)
 Fr=Fr+180;
end


end
%Eout=sqrt(E2(1,nL,1)^2+E2(1,nL,2)^2+E2(1,nL,3)^2);

if(plotWave)
    
    figure(7)
    
    plot(y,angs,colorAng);
    hold on
    
    %angsHomog=-0.5*y*k1*imag(epsb(1,3))/sqrt(epsb(1,1))*180/pi;
    %plot(y,angsHomog,'-k');
    
end

end

function FillKapaCylinderAntiSymDef(nGx,nGy,epsa,epsb,L,R,Na,a1,a2,d1)
global Kapa;

disp('anti-sym Defect');
nk=4;
if(length(epsa)==1)
nk=1;
end

invepsa1=inv(epsa);
invepsb1=inv(epsb);

if(nk==1)

invepsa=[invepsa1(1,1)];
invepsb=[invepsb1(1,1)];

else
invepsa=[invepsa1(1,1)  invepsa1(2,2) invepsa1(3,3) imag(invepsa1(1,3))];
invepsb=[invepsb1(1,1) invepsb1(2,2) invepsb1(3,3)  imag(invepsb1(1,3))];

end

ff=pi*R*R/(a1*a2);
by=pi/L;
bx=2*pi/a1;

KapaUnit=zeros(4*nGx+1,4*nGy+1,4)+1i*zeros(4*nGx+1,4*nGy+1,4);


for dGx=-2*nGx:2*nGx
    Gn=bx*dGx;
      dGxp=dGx+1+2*nGx;
    for dGy=-2*nGy:2*nGy
        Gm=by*dGy;
        GnmR=sqrt(Gn*Gn+Gm*Gm)*R;
        
        dGyp=dGy+1+2*nGy;
        
        if(dGx==0&& dGy==0)
            for k=1:nk
                  KapaUnit(dGxp,dGyp,k)=(invepsb(k)+ff*(invepsa(k)-invepsb(k)))/(2*Na);
                  
            end
            
         else
        tt=ff*besselj(1,GnmR)/GnmR;

            factb=0;
             if(dGx==0)
            ttb=dGy*by*a2/2;
            factb=a2/(2*L)*sin(ttb)/(ttb);
            end
            
            for k=1:nk
                    KapaUnit(dGxp,dGyp,k)= factb*invepsb(k)+tt*(invepsa(k)-invepsb(k))/(Na);;
            end
            
            
        end
        
        
    end
end

global ndef;
global defstart;




KapaDefect=zeros(1,4*nGy+1,4)+1i*zeros(1,4*nGy+1,4);
 


if(ndef>0)
 
     for dGy=-2*nGy:2*nGy
        
        dGyp=dGy+1+2*nGy;

        tt=dGy*by*a2/2;
          
          for k=1:nk
                if(dGy==0)
                  KapaDefect(1,dGyp,k)=invepsb(k)*a2/(2*L);
                else
                  KapaDefect(1,dGyp,k)=invepsb(k)*a2/(2*L)*sin(tt)/(tt);
             end
         
         end
    
        
    end
 
end

eps1=1;
eps3=1;

d1

if(d1>0)
 KappaEnds=zeros(1,4*nGy+1,4)+1i*zeros(1,4*nGy+1,4);

     for dGy=-2*nGy:2*nGy
        
        dGyp=dGy+1+2*nGy;

        tt=dGy*by*d1/2;
          
          for k=1:nk
                if(dGy==0)
                  KappaEnds(1,dGyp,k)=invepsb(k)*d1/(2*L);
                else
                  KappaEnds(1,dGyp,k)=invepsb(k)*d1/(2*L)*sin(tt)/(tt);
             end
         
         end
    
        
    end
 
end
 
 
Kapa=zeros(4*nGx+1,4*nGy+1,4)+1i*zeros(4*nGx+1,4*nGy+1,4);

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

 dGxp=1+2*nGx;

by=pi/L;
for n=-Na:Na-1 
 np=n+Na+1;

 for dGy=-2*nGy:2*nGy
       
      dGyp=dGy+1+2*nGy;
      
      if(n<0)
      sd=-d1;
      else
      sd=d1;
      end
        twindle=exp(-1i*((n+.5)*a2+sd)*by*dGy);

     if(isDef(np,1))

    for k=1:nk
       if(k~=4 || n>=0)
       Kapa(dGxp,dGyp,k)=  Kapa(dGxp,dGyp,k)+KapaDefect(1,dGyp,k)*twindle;
        else
       Kapa(dGxp,dGyp,k)=  Kapa(dGxp,dGyp,k)-KapaDefect(1,dGyp,k)*twindle;
        end
     end
   else
     for k=1:nk
       if(k~=4 || n>=0)
         Kapa(:,dGyp,k)=Kapa(:,dGyp,k)+KapaUnit(:,dGyp,k)*twindle;
        else
         Kapa(:,dGyp,k)=Kapa(:,dGyp,k)-KapaUnit(:,dGyp,k)*twindle;
        end
     end
   
   end
        
        
    end
  end

  
  if(d1>0)
  
   for dGy=-2*nGy:2*nGy
       
      dGyp=dGy+1+2*nGy;
        twindle1=exp(-1i*(d1/2)*by*dGy)+exp(-1i*(L-d1/2)*by*dGy);
        twindle2=exp(1i*(d1/2)*by*dGy)+exp(1i*(L-d1/2)*by*dGy);
twindle=twindle1+twindle2;

     for k=1:nk
         Kapa(dGxp,dGyp,k)=Kapa(dGxp,dGyp,k)+KappaEnds(1,dGyp,k)*twindle;
     end
           
    
  end
 
  end
 
%  Kapa(:,:,:)=KappaEnds(:,:,:);


  
end

function FillKapaCylinderAntiSym(nGx,nGy,epsa,epsb,L,R,Na,a1,a2,d1)
global Kapa;

disp('anti-sym');
nk=4;
if(length(epsa)==1)
nk=1;
end

invepsa1=inv(epsa);
invepsb1=inv(epsb);

if(nk==1)

invepsa=[invepsa1(1,1)];
invepsb=[invepsb1(1,1)];

else
invepsa=[invepsa1(1,1)  invepsa1(2,2) invepsa1(3,3) imag(invepsa1(1,3))];
invepsb=[invepsb1(1,1) invepsb1(2,2) invepsb1(3,3)  imag(invepsb1(1,3))];

end

ff=Na*pi*R*R/(a1*L);
d11=d1/L;

for dGx=-2*nGx:2*nGx
    Gn=2*pi*dGx/a1;
    for dGy=-2*nGy:2*nGy
        Gm=dGy*pi/L;
        GnmR=sqrt(Gn*Gn+Gm*Gm)*R;
        dGxp=dGx+1+2*nGx;
        dGyp=dGy+1+2*nGy;
        
        if(dGx==0&& dGy==0)
            for k=1:nk
                if(k~=4)
                    Kapa(dGxp,dGyp,k)=1*(ff*invepsa(k)+(1-ff)*invepsb(k))+2*d11*(invepsa(k)-invepsb(k));
                end
                
            end
            
        elseif(dGy==0)
            tt=2*ff*besselj(1,GnmR)/GnmR;
            
            for k=1:nk
                if(k~=4)
                    Kapa(dGxp,dGyp,k)=tt*(invepsa(k)-invepsb(k));
                end
            end
            
        else
            tt=a2*dGy*pi/(2*L);
            
            
            for k=1:nk
                if(k==4)
                    factm=-1i*sin(dGy*pi/2)*sin(Na*tt)/sin(tt)/Na;
                    vv=2*ff*besselj(1,GnmR)/GnmR;
                    
                    Kapa(dGxp,dGyp,k)=1*factm*vv*(invepsa(k)-invepsb(k));
                    if(dGx==0)
                        vvL=2*sin(dGy*pi*.5*(L-2*d1)/L)/(dGy*pi);
                        factm=-1i*sin(dGy*pi*.5);
                        
                        Kapa(dGxp,dGyp,k)=Kapa(dGxp,dGyp,k)+factm*vvL*invepsb(k);
                    end
                else
                    factm=cos(dGy*pi/2)*sin(Na*tt)/sin(tt)/Na;
                    vv=2*ff*factm*besselj(1,GnmR)/GnmR;
                    
                    Kapa(dGxp,dGyp,k)=1*vv*(invepsa(k)-invepsb(k));
                end
                
                
                if(dGx==0 && d11>0 && k~=4)
                    aa= pi*dGy*d11;
                    
                    zz=1*sin(aa)/aa*d11+2*cos(pi*dGy*(L-d1/2)/L)*sin(aa/2)/aa*d11;
                    Kapa(dGxp,dGyp,k)= Kapa(dGxp,dGyp,k)+zz*(invepsa(k)-invepsb(k));
                end
                
            end
            
        end
        
        
    end
end
end



function FillKapaRectangleAntiSym(nGx,nGy,epsa,epsb,L,Rx,Ry,Na,a1,a2,d1,fi)
global Kapa;

disp('rectangle anti-sym');
nk=4;
if(length(epsa)==1)
nk=1;
end

invepsa1=inv(epsa);
invepsb1=inv(epsb);

if(nk==1)

invepsa=[invepsa1(1,1)];
invepsb=[invepsb1(1,1)];

else
invepsa=[invepsa1(1,1)  invepsa1(2,2) invepsa1(3,3) imag(invepsa1(1,3))];
invepsb=[invepsb1(1,1) invepsb1(2,2) invepsb1(3,3)  imag(invepsb1(1,3))];

end

cell_size=a1*a2*Na*2;
rect_size=2*Rx*2*Ry;

fillx=2*Rx/a1;
filly=2*Ry/a2;

bx=2*pi/a1;
by=pi/L;


KapaUnit=zeros(1,4*nGy+1,4)+1i*zeros(1,4*nGy+1,4);



for dGx=-2*nGx:2*nGx
      dGxp=dGx+1+2*nGx;
   ttx=dGx*bx*Rx;
   for dGy=-2*nGy:2*nGy
        
        tty=dGy*by*Ry;
     
        dGyp=dGy+1+2*nGy;
        
  for k=1:nk
      if(dGx==0)
       four_coefx= (invepsb(k)+fillx*(invepsa(k)-invepsb(k)));
        else
          four_coefx=(invepsa(k)-invepsb(k))*2*Rx/(a1)*sin(ttx)/(ttx);
       end
        
         if(dGy==0)
       four_coefy= (invepsb(k)+filly*(invepsa(k)-invepsb(k)))/(2*Na);
        else
          four_coefy=2*(invepsa(k)-invepsb(k))*2*Ry/(a1)*sin(tty)/(tty)/(2*Na);
       end
      
       
        KapaUnit(dGxp,dGyp,k)=four_coefx*four_coefy;


  end

end

end

Kapa=zeros(4*nGx+1,4*nGy+1,4)+1i*zeros(4*nGx+1,4*nGy+1,4);


global ndef;
global defstart;

by=pi/L;


KapaDefect=zeros(1,4*nGy+1,4)+1i*zeros(1,4*nGy+1,4);
 


if(ndef>0)

     for dGy=-2*nGy:2*nGy
        
        dGyp=dGy+1+2*nGy;

        tt=dGy*by*a2/2;
          
          for k=1:nk
                if(dGy==0)
                  KapaDefect(1,dGyp,k)=invepsb(k)*a2/(2*L);
                else
                  KapaDefect(1,dGyp,k)=invepsb(k)*a2/(2*L)*sin(tt)/(tt);
             end
         
         end
    
        
    end
 
end
 

Kapa=zeros(4*nGx+1,4*nGy+1,4)+1i*zeros(4*nGx+1,4*nGy+1,4);

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

 dGxp=1+2*nGx;


by=pi/L;

for n=-Na:Na-1 
 np=n+Na+1;

 for dGy=-2*nGy:2*nGy
       
      dGyp=dGy+1+2*nGy;
        twindle=exp(-1i*(n+.5)*by*dGy*a2);
     if(isDef(np,1))

    for k=1:nk
       if(k~=4 || n>=0)
       Kapa(dGxp,dGyp,k)=  Kapa(dGxp,dGyp,k)+KapaDefect(1,dGyp,k)*twindle;
        else
       Kapa(dGxp,dGyp,k)=  Kapa(dGxp,dGyp,k)-KapaDefect(1,dGyp,k)*twindle;
        end
     end
   else
     for k=1:nk
       if(k~=4 || n>=0)
         Kapa(:,dGyp,k)=Kapa(:,dGyp,k)+KapaUnit(:,dGyp,k)*twindle;
        else
         Kapa(:,dGyp,k)=Kapa(:,dGyp,k)-KapaUnit(:,dGyp,k)*twindle;
        end
     end
   
   end
        
        
    end
  end
  
end



function FillKapaCylinderSym(nGx,nGy,epsa,epsb,L,R,Na,a1,a2,d1)
global Kapa;



disp('sym');
nk=4;
if(length(epsa)==1)
nk=1;
end

invepsa1=inv(epsa);
invepsb1=inv(epsb);

if(nk==1)

invepsa=[invepsa1(1,1)];
invepsb=[invepsb1(1,1)];

else
invepsa=[invepsa1(1,1)  invepsa1(2,2) invepsa1(3,3) imag(invepsa1(1,3))];
invepsb=[invepsb1(1,1) invepsb1(2,2) invepsb1(3,3)  imag(invepsb1(1,3))];

end

ff=Na*pi*R*R/(a1*L);
d11=d1/L;

for dGx=-2*nGx:2*nGx
    Gn=2*pi*dGx/a1;
    for dGy=-2*nGy:2*nGy
        Gm=dGy*pi/L;
        GnmR=sqrt(Gn*Gn+Gm*Gm)*R;
        dGxp=dGx+1+2*nGx;
        dGyp=dGy+1+2*nGy;
        
        if(dGx==0&& dGy==0)
            for k=1:nk
                    Kapa(dGxp,dGyp,k)=1*(ff*invepsa(k)+(1-ff)*invepsb(k))+2*d11*(invepsa(k)-invepsb(k));
                
            end
            
        elseif(dGy==0)
            tt=2*ff*besselj(1,GnmR)/GnmR;
            
            for k=1:nk
                    Kapa(dGxp,dGyp,k)=tt*(invepsa(k)-invepsb(k));
            end
            
        else
            tt=a2*dGy*pi/(2*L);
            
            
            for k=1:nk

                    factm=cos(dGy*pi/2)*sin(Na*tt)/sin(tt)/Na;
                    vv=2*ff*factm*besselj(1,GnmR)/GnmR;
                    
                    Kapa(dGxp,dGyp,k)=1*vv*(invepsa(k)-invepsb(k));

                
                if(dGx==0 && d11>0 )
                    aa= pi*dGy*d11;
                    
                    zz=1*sin(aa)/aa*d11+2*cos(pi*dGy*(L-d1/2)/L)*sin(aa/2)/aa*d11;
                    Kapa(dGxp,dGyp,k)= Kapa(dGxp,dGyp,k)+zz*(invepsa(k)-invepsb(k));
                end
                
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
% The function plotting 3-dimensional arrow
%
% h=arrow3d(x,y,z,head_frac,radii,radii2,colr)
%
% The inputs are:
%       x,y,z =  vectors of the starting point and the ending point of the
%           arrow, e.g.:  x=[x_start, x_end]; y=[y_start, y_end];z=[z_start,z_end];
%       head_frac = fraction of the arrow length where the head should  start
%       radii = radius of the arrow
%       radii2 = radius of the arrow head (defult = radii*2)
%       colr =   color of the arrow, can be string of the color name, or RGB vector  (default='blue')
%
% The output is the handle of the surfaceplot graphics object.
% The settings of the plot can changed using: set(h, 'PropertyName', PropertyValue)
%
% example #1:
%        arrow3d([0 0],[0 0],[0 6],.5,3,4,[1 0 .5]);
% example #2:
%        arrow3d([2 0],[5 0],[0 -6],.2,3,5,'r');
% example #3:
%        h = arrow3d([1 0],[0 1],[-2 3],.8,3);
%        set(h,'facecolor',[1 0 0])
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
    R1=Pc+radii*cos(theta).*(P1-Pc) + radii*sin(theta).*cross(dr,(P1-Pc)) -origin_shift;
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
    R1=Pc+radii2*cos(theta).*(P1-Pc) + radii2*sin(theta).*cross(dr,(P1-Pc)) -origin_shift;
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


function FillKapaAntiSymNum(geometry,nGx,nGy,epsa,epsb,L,Rx,Ry,Na,a1,a2,d1,fi)

disp('Numerical anti-sym');
nk=4;
if(length(epsa)==1)
nk=1;
end

invepsa1=inv(epsa);
invepsb1=inv(epsb);

if(nk==1)

invepsa=[invepsa1(1,1)];
invepsb=[invepsb1(1,1)];

else
invepsa=[invepsa1(1,1)  invepsa1(2,2) invepsa1(3,3) imag(invepsa1(1,3))];
invepsb=[invepsb1(1,1) invepsb1(2,2) invepsb1(3,3)  imag(invepsb1(1,3))];

end

global Kapa;

 Kapa=zeros(4*nGx+1,4*nGy+1,4)+1i*zeros(4*nGx+1,4*nGy+1,4);


ngridx=40;
if(nGx==0) 
ngridx=2;
end
ngridy=40;

invep=zeros(ngridx,ngridy,nk);
rotMat=[cosd(fi) sind(fi);-sind(fi) cosd(fi)];


xyr=[0  0]';
%The following loop carries out the definition of the unit cell.

    nx=1;
    for countX=-a1/2:a1/ngridx:a1/2
        ny=1;
        for countY=-a2/2:a2/ngridy:a2/2
            %Following condition allows to define the circle with of radius r
			xy=[countX countY]';
			xyr=rotMat*xy;
			countXr=xyr(1);
           	countYr=xyr(2);
				
            if(geometry==0)
                inside=(countXr/Rx)^2+(countYr/Ry)^2<1;
            else
                 inside =false;
                  if(Rx>0 && Ry>0) 
                     inside=abs(countXr/Rx)<=1 && abs(countYr/Ry)<=1;
                  end
            end
            
            if(inside)
                
                ip=invepsa;

            else
                
                ip=invepsb;
 
            end
            
            
            
            for k=1:nk
                inveps(nx,ny,k)=ip(k);
            end
            
            xSet(nx)=countX;
            ySet(ny)=countY;
            ny=ny+1;
        end
        nx=nx+1;
    end



MtNt=(length(xSet)-1)*(length(ySet)-1)*2*Na;
%The next loop computes the Fourier expansion coefficients
bx=2*pi/a1;

by=pi/L;

KapaUnit=zeros(4*nGx+1,4*nGy+1,4)+1i*zeros(4*nGx+1,4*nGy+1,4);
 
for dGx=-2*nGx:2*nGx
    
    for dGy=-2*nGy:2*nGy
        
        
        dGxp=dGx+1+2*nGx;
        dGyp=dGy+1+2*nGy;
        
        for nx=1:length(xSet)-1
            for ny=1:length(ySet)-1
                x=xSet(nx);
                y=ySet(ny);
                tt=dGx*bx*x+dGy*by*y;
                
                for k=1:nk
                    KapaUnit(dGxp,dGyp,k)=KapaUnit(dGxp,dGyp,k)+inveps(nx,ny,k)*exp(-1i*(tt));
                end
            end
        end
        
        
    end
end

KapaUnit=  KapaUnit/MtNt;


global ndef;
global defstart;

KapaDefect=zeros(1,4*nGy+1,4)+1i*zeros(1,4*nGy+1,4);
 

if(ndef>0)
 
     for dGy=-2*nGy:2*nGy
        
        dGyp=dGy+1+2*nGy;

        tt=dGy*by*a2/2;
          
          for k=1:nk
                if(dGy==0)
                  KapaDefect(1,dGyp,k)=invepsb(k)*a2/(2*L);
                else
                  KapaDefect(1,dGyp,k)=invepsb(k)*a2/(2*L)*sin(tt)/(tt);
             end
         
         end
    
        
    end
 
end


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

 dGxp=1+2*nGx;


by=pi/L;
for n=-Na:Na-1 
 np=n+Na+1;

 for dGy=-2*nGy:2*nGy
       
      dGyp=dGy+1+2*nGy;
        twindle=exp(-1i*(n+.5)*by*dGy*a2);
     if(isDef(np,1))

    for k=1:nk
       if(k~=4 || n>=0)
       Kapa(dGxp,dGyp,k)=  Kapa(dGxp,dGyp,k)+KapaDefect(1,dGyp,k)*twindle;
        else
       Kapa(dGxp,dGyp,k)=  Kapa(dGxp,dGyp,k)-KapaDefect(1,dGyp,k)*twindle;
        end
     end
   else
     for k=1:nk
       if(k~=4 || n>=0)
         Kapa(:,dGyp,k)=Kapa(:,dGyp,k)+KapaUnit(:,dGyp,k)*twindle;
        else
         Kapa(:,dGyp,k)=Kapa(:,dGyp,k)-KapaUnit(:,dGyp,k)*twindle;
        end
     end
   
   end
        
        
        
    end

end

d1

if(d1>0)
 KappaEnds=zeros(1,4*nGy+1,4)+1i*zeros(1,4*nGy+1,4);

     for dGy=-2*nGy:2*nGy
        
        dGyp=dGy+1+2*nGy;

        tt=dGy*by*d1/2;
          
          for k=1:nk
                if(dGy==0)
                  KappaEnds(1,dGyp,k)=invepsb(k)*d1/(2*L);
                else
                  KappaEnds(1,dGyp,k)=invepsb(k)*d1/(2*L)*sin(tt)/(tt);
               end
         
         end
    
        
    end
 
end

  if(d1>0)
  
   for dGy=-2*nGy:2*nGy
       
      dGyp=dGy+1+2*nGy;
        twindle1=exp(-1i*(d1/2)*by*dGy)+exp(-1i*(L-d1/2)*by*dGy);
        twindle2=exp(1i*(d1/2)*by*dGy)+exp(1i*(L-d1/2)*by*dGy);
twindle=twindle1+twindle2;

     for k=1:nk
         Kapa(dGxp,dGyp,k)=Kapa(dGxp,dGyp,k)+KappaEnds(1,dGyp,k)*twindle;
     end
           
    
  end
 
  end


 

end

