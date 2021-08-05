function contour

[filename,filepath]=uigetfile('*.txt', 'Select contour file')
file=strcat(filepath,strcat('\',filename));

fid= fopen(filename);
line=getNewDataLine(fid);
X = str2num(line);
M=length(X);
line=getNewDataLine(fid);
Y = str2num(line);
N=length(Y);

[x,y] = meshgrid(X,Y);

z=zeros(N,M);

for i=1:N
line=getNewDataLine(fid);
numbs=str2num(line);;
for j=1:M
z(i,j) =numbs(j);
end
end

colormap('jet');                    % create a new figure
contourf(x,y,z);                  % contour plot
colorbar;                             % add a colorbar
axis on;  
set(gca,'Xtick',X)
set(gca,'Ytick',Y)


function line=getNewDataLine(fid)

TF=1;
k=0;
while (k<100 && TF==1)
line=fgets(fid);
TF = strncmpi(line,'/',1);
k=k+1;

end