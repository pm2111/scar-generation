addpath(genpath('C:/Users/petnov/Dropbox/'));
addpath(genpath('C:/Users/petnov/Dropbox/qrs/'));
addpath(genpath('C:/Users/petnov/Dropbox/sim 3d scripts/'));
addpath(genpath('C:/Users/petnov/Dropbox/shared - copy/mesh scripts/'));
close all
addpath('C:\Users\petnov\Dropbox\sim3d scripts\');
addpath(genpath('C:\Users\petnov\Dropbox\mesh scripts\'));
addpath(genpath('C:\Users\petnov\Dropbox\qrs\'));
addpath(genpath('C:\Users\petnov\Dropbox\tools\'));

addpath('C:\Users\petnov\Dropbox\shared - copy\IO\');
patient_sim_folder='D:\ARVC meshing automatic\patients\DTI003\';

H_mesh=struct;
H_mesh.xyz=csvread('D:\ARVC meshing automatic\patients\DTI003\DTI003_mesh\DTI003_coarse_xyz.csv');
H_mesh.tri=csvread('D:\ARVC meshing automatic\patients\DTI003\DTI003_mesh\DTI003_coarse_tri.csv');
H_mesh.rv=csvread('D:\ARVC meshing automatic\patients\DTI003\DTI003_mesh\DTI003_coarse_rvface.csv');
H_mesh.lv=csvread('D:\ARVC meshing automatic\patients\DTI003\DTI003_mesh\DTI003_coarse_lvface.csv');
H_mesh.epi=csvread('D:\ARVC meshing automatic\patients\DTI003\DTI003_mesh\DTI003_coarse_epiface.csv');
H_mesh.face=vertcat(H_mesh.epi,H_mesh.lv,H_mesh.rv);
        

names=dir('D:/ARVC meshing automatic/patients/DTI003/results/plots/');

for i=3:size(names,1)
    no_lge=load(strcat('D:/ARVC meshing automatic/patients/DTI003/results/plots/',names(i).name));
 figure()
 hold on
 patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',no_lge,'facecolor','interp')
 colorbar
 axis off
   view([77 67])
 axis equal
 axis off
 colorbar off   
  headlight
saveas(gcf,strcat('C:\Users\petnov\Dropbox\figures_repo\figures_paperII\scars_LV_anterior\',names(i).name,'.png'))
end