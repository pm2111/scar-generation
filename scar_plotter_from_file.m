

names=dir('D:/ARVC meshing automatic/patients/DTI003/results/plots/');

for i=3:size(names,1)
    no_lge=load(strcat('D:/ARVC meshing automatic/patients/DTI003/results/plots/',names(i).name));
 figure()
 headlight
 hold on
 patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',no_lge,'facecolor','interp')
 colorbar
 axis off
   view([-125 -55])
 axis equal
 axis off
 colorbar off             
saveas(gcf,strcat('C:/Users/petnov/Dropbox/figures_repo/scar_intersheet_scar_CV_modulation_20%_RV_mid_size_',num2str(width_intersheet),'.png'))
endz