M = readmatrix('E:/Ph.D/实验数据/分子动力学模拟/4ZLJ-dtac/2nd/gibbs~.txt');
x = M(:,1);
y = M(:,2);
z = M(:,3);
[xq,yq] = meshgrid(0.49:.00003:0.61, 3.586:.00003:3.693);
zq = griddata(x,y,z,xq,yq);
figure(1);
mesh(xq,yq,zq)
hold on
surf(xq,yq,-4+0*zq,zq,'EdgeColor','none');

cmap = [68 4 90;65 62 133;48 104 141;31 146 139;53 183 119; 145 213 66;255 246 143]/255;
cn = 256;
custom_cmap = interp1(linspace(0,1,size(cmap,1)),cmap,linspace(0,1,cn));
colormap(custom_cmap)
colorbar
hold on

set(gca,'FontName','Arial','FontSize',12)

xlabel('x');
ylabel('y');
zlabel('Gibbs free energy (kJ mol^{-1})');
xl=xlabel('RMSD (nm)');
set(xl,'rotation',13)
yl=ylabel('Gyrate (nm)');

f1 = figure(1);

print(f1,'E:/Ph.D/实验数据/分子动力学模拟/4ZLJ-dtac/2nd/gibbs-dtac.pdf','-dpdf','-r1600')



M = readmatrix('E:/Ph.D/实验数据/分子动力学模拟/4ZLJ-sds/2nd/gibbs~.txt');
x = M(:,1);
y = M(:,2);
z = M(:,3);
[xq,yq] = meshgrid(0.53:.00003:0.69, 3.609:.00003:3.687);
zq = griddata(x,y,z,xq,yq);
figure(2);
mesh(xq,yq,zq)
hold on
surf(xq,yq,-4+0*zq,zq,'EdgeColor','none');

cmap = [68 4 90;65 62 133;48 104 141;31 146 139;53 183 119; 145 213 66;255 246 143]/255;
cn = 256;
custom_cmap = interp1(linspace(0,1,size(cmap,1)),cmap,linspace(0,1,cn));
colormap(custom_cmap)
colorbar
hold on

set(gca,'FontName','Arial','FontSize',12)

xlabel('x');
ylabel('y');
zlabel('Gibbs free energy (kJ mol^{-1})');
xl=xlabel('RMSD (nm)');
set(xl,'rotation',13)
yl=ylabel('Gyrate (nm)');

f2 = figure(2);

print(f2,'E:/Ph.D/实验数据/分子动力学模拟/4ZLJ-sds/2nd/gibbs-sds.pdf','-dpdf','-r1600')



M = readmatrix('E:/Ph.D/实验数据/分子动力学模拟/4ZLJ-trix/gibbs~.txt');
x = M(:,1);
y = M(:,2);
z = M(:,3);
[xq,yq] = meshgrid(0.749:.00003:0.884, 3.605:.00003:3.695);
zq = griddata(x,y,z,xq,yq);
figure(3);
mesh(xq,yq,zq)
hold on
surf(xq,yq,-4+0*zq,zq,'EdgeColor','none');

cmap = [68 4 90;65 62 133;48 104 141;31 146 139;53 183 119; 145 213 66;255 246 143]/255;
cn = 256;
custom_cmap = interp1(linspace(0,1,size(cmap,1)),cmap,linspace(0,1,cn));
colormap(custom_cmap)
colorbar
hold on

set(gca,'FontName','Arial','FontSize',12)

xlabel('x');
ylabel('y');
zlabel('Gibbs free energy (kJ mol^{-1})');
xl=xlabel('RMSD (nm)');
set(xl,'rotation',13)
yl=ylabel('Gyrate (nm)');

f3 = figure(3);

print(f3,'E:/Ph.D/实验数据/分子动力学模拟/4ZLJ-trix/gibbs-trix.pdf','-dpdf','-r1600')