
addpath(genpath('C:\Users\petnov\Dropbox\shared - copy\'));
addpath(genpath('C:\Users\petnov\Dropbox\mesh scripts'))
addpath('C:\Users\petnov\Dropbox\shared - copy\MESHES\');
addpath('C:\Users\petnov\Dropbox\shared - copy\IO\');
addpath('C:\Users\petnov\Dropbox\shared - copy\IO\'); 
close all
enableVTK;

    %mkdir 'D:\ARVC meshing automatic\patients\patient01\regions';
    patient_sim_folder =  strcat('D:\ARVC meshing automatic\patients\DTI003\');

    %H0 = read_VTK( fullfile( TORSO_MODEL_DIR , 'HEART_with_UPS_louie.vtk' ) ); %this is louie's heart
H_mesh=struct;
H_mesh.xyz=csvread('D:\ARVC meshing automatic\patients\DTI003\DTI003_mesh\DTI003_coarse_xyz.csv');
H_mesh.tri=csvread('D:\ARVC meshing automatic\patients\DTI003\DTI003_mesh\DTI003_coarse_tri.csv');
H_mesh.rv=csvread('D:\ARVC meshing automatic\patients\DTI003\DTI003_mesh\DTI003_coarse_rvface.csv');
H_mesh.lv=csvread('D:\ARVC meshing automatic\patients\DTI003\DTI003_mesh\DTI003_coarse_lvface.csv');
H_mesh.epi=csvread('D:\ARVC meshing automatic\patients\DTI003\DTI003_mesh\DTI003_coarse_epiface.csv');
H_mesh.face=vertcat(H_mesh.epi,H_mesh.lv,H_mesh.rv);
  H_mesh.xyz=H_mesh.xyz;
%    mesh_cell_property2 =  load(strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'AHA_regions.mat'));
        
  %  mesh_cell_property2= mesh_cell_property2.mesh_cell_property2;
    
%     %scar type 5 - 1 conducting channel with
%     non uniform density borders and varying endo-epi surfaces
%      specify the params of the algorithm in terms of 
   % position of the scar, select  using cursor selection tool (select data tips (2nd row on the top right plot menu in Matlab
   % 2019, select point and right click on export cursor Data to workspace)
   % in matlab)
   % loop variable defines conduction slowing in surface of border of scar region 
     

 
 distance= @(vec,vec1) ((vec(:,1)-vec1(:,1)).^2+(vec(:,2)-vec1(:,2)).^2+(vec(:,3)-vec1(:,3)).^2).^0.5;

 figure()
 headlight
 hold on
 patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','facecolor','b')
 colorbar
 axis off
 
 H_mesh.triORTHOfull= load('D:\ARVC meshing automatic\patients\DTI003\DTI003_mesh\DTI003_coarse_tetrahedron_centreFibers.csv');
 H_mesh.triORTHO =H_mesh.triORTHOfull(:,1:3); 
 %use the apex-base dir to define scar thickness
 fibers_normal = H_mesh.triORTHOfull(:,7:9);
 %scar type 6 - 1 conducting channel with
%     non uniform density borders oriented along fiber direction

       
   
       
RV_mid=0;
RV_posterior=0;
RV_anterior=0;
%width=1/2;
rp=0.5; %1.5 for RV ant
Ni=1;
Np=15;
        for ri = [0.001]
        for j=6
            for density=2
                %for t=1:8
                  %  for r=1:8
                          k=0;
                        %  width=12;
                       
%                          if RV_anterior ==1
% %                             thickness_scar = thickness_scar*10;
% %                             width_intersheet=width_intersheet/5;
%                                 width =0.7
%                         end
                       %  scar_properties ={};
            if RV_posterior==1
                x = -0.989;
                 y = -9.9524;
                 z = -9.4597;
            elseif RV_mid ==1 
                x= 1.6; %RV post [x,y,z]=[-0.989, -9.9524,-9.4597]
                     y=    -11.6;
                     z= -6.19;
                     
            elseif RV_anterior==1
                x=-0.83;
                y=-11.14;
                z=-6.5;
            else
                x=cursor_info.Position(1);
                y=cursor_info.Position(2);
                z=cursor_info.Position(3);
               
            end
            
             ts_scaling_no_fib=ones(23,1);

               epi1= knnsearch(H_mesh.xyz,[x,y,z]);
         endo1=knnsearch(H_mesh.xyz(H_mesh.rv,:),H_mesh.xyz(epi1,:));
         
         %mid1= (H_mesh.xyz(H_mesh.epi(epi1),:)+H_mesh.xyz(H_mesh.rv(endo1),:))/2;
         epi1_coords= H_mesh.xyz(epi1,:);
         endo1_coords = H_mesh.xyz(endo1,:);
        AB = epi1_coords-endo1_coords;
        mid1 = endo1+AB/2; %finding a mid point between two points, A(endo point) and B (epi)

         mid1_id= knnsearch(H_mesh.xyz,mid1);
                  points = [mid1_id]; %can change this to something more specific if nec

         points = [epi1];
         mid1 =[x y z];
             for i=1:Np
                 [rw,col]=find(H_mesh.tri==mid1_id);
                  fiber_dir = -H_mesh.triORTHOfull(col(1),1:3);%use the direction in the fiber dir
                  if RV_anterior ==1
                      fiber_dir=-fiber_dir;
                  end
                % fiber_dir = -fibers_normal(col(1),1:3);
                 mid1 = mid1+(fiber_dir*i)/Np;
                  id_next=knnsearch(H_mesh.xyz,mid1);
                  
                  points = [points,id_next];
                 dists= distance(H_mesh.xyz,mid1);
                 keep=find(dists<rp);%keep the points which are a radial dist of rp away
                points = horzcat(points,keep');
                points=unique(points);
                %now for the same point, jump in the intersheet direction
                %and collect the points
                intersheet_dir = -H_mesh.triORTHOfull(id_next,7:9);
                for j=1:Ni
                    dummy_id_next = id_next;
                next_point_intersheet = H_mesh.xyz(dummy_id_next,:)+intersheet_dir*j;
                id_next_intersheet=knnsearch(H_mesh.xyz,next_point_intersheet);
                dists_intersheet = distance(H_mesh.xyz,H_mesh.xyz( id_next_intersheet,:));
                points_intersheet = find(dists_intersheet <ri);
                points =  horzcat(points,points_intersheet');
                points =unique(points);
                
                dummy_id_next = id_next_intersheet;
                end
                  mid1_id=id_next;
                  mid1= H_mesh.xyz(id_next,:);
             end
                
                  
                %find if there are any points in the LV
                no_lge= ones(size(H_mesh.xyz,1),1);

                indices_on_lv_surface=find(no_lge(H_mesh.lv) >1);
                %transform mesh to align it with z dir vertically
                
                
                              no_lge= ones(size(H_mesh.xyz,1),1);
                     no_lge(points)=density;
                     
                if exist('delete_points')==1
                 no_lge=points_cleaner(H_mesh,delete_points,no_lge);
               end
 figure()
 headlight
 hold on
 patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',no_lge,'facecolor','interp')
 colorbar
 axis off
   view([-124 -62])
 axis equal
 axis off
 colorbar off    
 %scatter3(mid1(1),mid1(2),mid1(3),105,'filled')
 %scatter3(endo1_coords(1),endo1_coords(2),endo1_coords(3),105,'filled','r')
 scatter3(epi1_coords(1),epi1_coords(2),epi1_coords(3),105,'filled','g')
saveas(gcf,strcat('C:\Users\petnov\Dropbox\figures_chapter8\scar_RV_post_rp_',num2str(rp),'_Np_',num2str(Np),'_ri_',num2str(ri),'_Ni_',num2str(Ni),'.png'))

           
            no_lge(points)=density;
            % no_lge(points)=i*dist_to_point_norm;
            dummy=find(no_lge)<1;
            no_lge(dummy)=1;
            
      
            if RV_anterior==1
                saveas(gcf,strcat('C:\Users\petnov\Dropbox\figures_chapter8\scar_RV_ant_intersheet_size_',num2str(width_intersheet),'.png'))
               %dlmwrite(strcat(patient_sim_folder,'scars_RV_anterior\scar_CV_slowing_factor',num2str(density),'_initial_length_','1','_growth_factor_',sprintf('%.1f',double(width_intersheet)),'.csv'), no_lge,'precision',10);
            elseif RV_mid==1
                saveas(gcf,strcat('C:\Users\petnov\Dropbox\figures_chapter8\scar_RV_mid_intersheet_size_',num2str(width_intersheet),'.png'))

               % dlmwrite(strcat(patient_sim_folder,'scars_RV_mid\scar_CV_slowing_factor',num2str(density),'_initial_length_','1','_growth_factor_',sprintf('%.1f',double(width_intersheet)),'.csv'), no_lge,'precision',10);
               
            elseif RV_posterior==1
                saveas(gcf,strcat('C:\Users\petnov\Dropbox\figures_chapter8\scar_RV_inf_intersheet_size_',num2str(width_intersheet),'.png'))
              %  dlmwrite(strcat(patient_sim_folder,'scars_RV_inferior\scar_CV_slowing_factor',num2str(density),'_initial_length_','1','_growth_factor_',sprintf('%.1f',double(width_intersheet)),'.csv'), no_lge,'precision',10);
            end
           
           % dlmwrite(strcat(patient_sim_folder,'scalings\eikonal06_coarse_fib_',num2str(0),'.csv'), no_lge,'precision',10);
        %   dlmwrite(strcat(patient_sim_folder,'scars_RV_inferior\scar_CV_slowing_factor',num2str(density),'_initial_length_','1','_growth_factor_',sprintf('%.1f',double(width)),'.csv'), no_lge,'precision',10);
         %   dlmwrite(strcat(patient_sim_folder,'scalings\eikonal06_coarse_fib_',num2str(2),'_scar.csv'),  int_lge,'precision',10);
          %  dlmwrite(strcat(patient_sim_folder,'scalings\eikonal06_coarse_fib_',num2str(3),'_scar.csv'),  adv_lge,'precision',10);
           end
            end
        end
%         
%         
%  figure()
%  headlight
%  hold on
%  patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',no_lge,'facecolor','interp')
%  colorbar
%  axis off
%  
%  H_meshvtk=struct;
%  H_meshvtk.xyz = H_mesh.xyz;
%  H_meshvtk.tri = H_mesh.tri;
%  H_meshvtk.xyzFibers = H_mesh.triORTHO(:,1:3);
%  write_VTK(H_meshvtk,'D:/ARVC meshing automatic/patients/DTI003/DTI003.vtk');
% 
%  
% 
% 
%  figure()
%  plot
%  
%  xs=linspace(0,1,50);
%  figure()
%  plot(xs,y(xs,2),'LineWidth',3)
%  hold on
%   plot(xs,y(xs,4),'LineWidth',3)
%  plot(xs,y(xs,8),'LineWidth',3)
%  legend({'typical ACM scar','pronounced ACM scar', 'severe ACM scar'},'Fontsize',35)
%  xlabel('normalised distance from middle of scar to border','Fontsize',35)
%  ylabel('conduction scaling factor','Fontsize',35)
%  set(gca,'Fontsize',25)
%  
%  
%  
%  
%  
%  scar_old_long=load('D:/ARVC meshing automatic/patients/patient06/scalings/eikonal06_coarse_fib_0_scar_along_fiber_slowing_4_length_10_radius_9.csv');
%  
%  positions=load('D:/ARVC meshing automatic/patients/patient06/mpp/ECG_ELECTRODES.mat');
%  
%  pos2=load('D:/ARVC meshing automatic/patients/patient06/DTI003_electrodePositions.csv');
%  
%  
%  
%  figure()
%  headlight
%  hold on
%  patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',scar_old_long,'facecolor','interp')
%  colorbar

       