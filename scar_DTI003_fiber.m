
addpath(genpath('C:\Users\petnov\Dropbox\shared - copy\'));
addpath(genpath('C:\Users\petnov\Dropbox\mesh scripts'))
addpath('C:\Users\petnov\Dropbox\shared - copy\MESHES\');
addpath('C:\Users\petnov\Dropbox\shared - copy\IO\');


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
 
 H_mesh.triORTHO= load('D:\ARVC meshing automatic\patients\DTI003\DTI003_mesh\DTI003_coarse_tetrahedron_centreFibers.csv');
 H_mesh.triORTHO =H_mesh.triORTHO(:,1:3); 
 
 %scar type 6 - 1 conducting channel with
%     non uniform density borders oriented along fiber direction
    for length=[3,10]
        for width =3:9
        for j=6
            for density=6:10
                %for t=1:8
                  %  for r=1:8
                          k=0;
                          width=12;
                          length=4;
                        length=length/9;
                        width=width/10;
                        
                       %  scar_properties ={};

                x= cursor_info.Position(1); %RV post [x,y,z]=[-0.989, -9.9524,-9.94597]
                     y=    cursor_info.Position(2);
                     z= cursor_info.Position(3);

             ts_scaling_no_fib=ones(23,1);

               epi1= knnsearch(H_mesh.xyz(  H_mesh.epi,:),[x,y,z]);
         endo1=knnsearch(H_mesh.xyz(H_mesh.rv,:),H_mesh.xyz(epi1,:));
         
         mid1= (H_mesh.xyz(H_mesh.epi(epi1),:)+H_mesh.xyz(H_mesh.rv(endo1),:))/2;
         
         mid1_id= knnsearch(H_mesh.xyz,mid1);
         
         points = [mid1_id];
         for i=1:50
             [rw,col]=find(H_mesh.tri==mid1_id);
             fiber_dir = -H_mesh.triORTHO(col(1),1:3);
             mid1 = mid1+fiber_dir*(length/3);
              id_next=knnsearch(H_mesh.xyz,mid1);
              points = [points,id_next];
             dists= distance(H_mesh.xyz,mid1);
             keep=find(dists<width);
            points = horzcat(points,keep');
            points=unique(points);

              mid1_id=id_next;
              mid1= H_mesh.xyz(id_next,:);
         end
         
         
         scar_dists=distance(H_mesh.xyz,mid1);
             scar_ids=find(scar_dists<4.75*j); %semicircle of radius 1cm
             border_ids=find(scar_dists>3*j & scar_dists <5*j);

             
             scar_ids2= find(scar_dists <3*j);
             border_ids2= find(scar_dists <3*j);
%            ts_scaling_no_fib(12:21)= ts_scaling_no_fib(12:21)*1/RV_scal(j);
%            ts_scaling_no_fib(1:11)= ts_scaling_no_fib(1:11)*1/LV_scal(i);
% 
%            ts_scaling([12:21])= ts_scaling([12:21])*1/RV_scal(j);
%                       ts_scaling([1:11])= ts_scaling([1:11])*1/LV_scal(i);
% 
%            ts_scaling_int([12:21])= ts_scaling_int([12:21])*1/RV_scal(j);
%                       ts_scaling_int([1:11])= ts_scaling_int([1:11])*1/LV_scal(i);
% 
%            ts_scaling_strong([12:21])= ts_scaling_strong([12:21])*1/RV_scal(j);
%            ts_scaling_strong([1:11])= ts_scaling_strong([1:11])*1/LV_scal(i);
% 
%            
           no_lge= ones(size(H_mesh.xyz,1),1);
           
            dist_to_point = distance(mid1,H_mesh.xyz(border_ids,:));
            dist_to_point_norm= dist_to_point/max(dist_to_point);
            
            
            no_lge(points)=density/2;
            % no_lge(points)=i*dist_to_point_norm;
            dummy=find(no_lge)<1;
            no_lge(dummy)=1;
            
            dist_to_point2 = distance([x,y,z],H_mesh.xyz(border_ids2,:));
            dist_to_point_norm2= dist_to_point2/max(dist_to_point2);
            
          %  no_lge(border_ids2)=i*dist_to_point_norm2;
            dummy=find(no_lge)<1;
            no_lge(dummy)=1;       
            
            
%             for l=1:size(border_ids)
%                 dummy=find(H_mesh.epi == border_ids(l));
%                 if isempty(dummy)==0
%                     no_lge(border_ids(l))=1;
%                 end
%             end

             %no_lge(H_mesh.epi)=1;
            
%              for l=1:size(border_ids)
%                 dummy=find(H_mesh.epi == border_ids(l));
%                 if isempty(dummy)==0
%                     no_lge(border_ids(l))=t*dist_to_point_norm(l);
%                 end
%             end
%                for l=1:size(border_ids)
%                 dummy=find(H_mesh.rv == border_ids(l));
%                 if isempty(dummy)==0
%                     no_lge(border_ids(l))=r*dist_to_point_norm(l);
%                 end
%             end
%             dummy=find(no_lge)<1;
%             no_lge(dummy)=1;
%             
%             if k==1
%                 no_lge(H_mesh.epi)=1;
%             end
%             light_lge(scar_ids)=4*dist_to_point_norm;
%             dummy=find(light_lge)<1;
%             light_lge(dummy)=1;
%             
%             
%             
%            % light_lge(border_ids)=20*i;
% 
%              int_lge(scar_ids)=4*max(int_lge);
%             adv_lge(scar_ids)=4*max(adv_lge);

%            no_lge(H_mesh.rv)=1;
%            light_lge(H_mesh.rv)=1;
%            int_lge(H_mesh.rv)=1;
%            adv_lge(H_mesh.rv)=1;
%            
%            no_lge(H_mesh.lv)=1;
%            light_lge(H_mesh.lv)=1;
%            int_lge(H_mesh.lv)=1;
%            adv_lge(H_mesh.lv)=1;
        
           
           % dlmwrite(strcat(patient_sim_folder,'scalings\eikonal06_coarse_fib_',num2str(0),'.csv'), no_lge,'precision',10);
            dlmwrite(strcat(patient_sim_folder,'scalings\eikonal06_coarse_fib_',num2str(0),'_scar_along_fiber_slowing_',num2str(density/2),'_length_',num2str(length),'_radius_',num2str(width),'.csv'), no_lge,'precision',10);
         %   dlmwrite(strcat(patient_sim_folder,'scalings\eikonal06_coarse_fib_',num2str(2),'_scar.csv'),  int_lge,'precision',10);
          %  dlmwrite(strcat(patient_sim_folder,'scalings\eikonal06_coarse_fib_',num2str(3),'_scar.csv'),  adv_lge,'precision',10);
           end
            end
        end
    end
      %  end
       % end
        

 figure()
 headlight
 hold on
 patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',no_lge,'facecolor','interp')
 colorbar
 axis off
 
 
 
 xs=linspace(0,1,50);
 figure()
 plot(xs,y(xs,2),'LineWidth',3)
 hold on
  plot(xs,y(xs,4),'LineWidth',3)
 plot(xs,y(xs,8),'LineWidth',3)
 legend({'typical ACM scar','pronounced ACM scar', 'severe ACM scar'},'Fontsize',35)
 xlabel('normalised distance from middle of scar to border','Fontsize',35)
 ylabel('conduction scaling factor','Fontsize',35)
 set(gca,'Fontsize',25)
 
 
 scar_old_long=load('D:/ARVC meshing automatic/patients/patient06/scalings/eikonal06_coarse_fib_0_scar_along_fiber_slowing_4_length_10_radius_9.csv');
 
 positions=load('D:/ARVC meshing automatic/patients/patient06/mpp/ECG_ELECTRODES.mat');
 
 pos2=load('D:/ARVC meshing automatic/patients/patient06/DTI003_electrodePositions.csv');
 
 
 
 figure()
 headlight
 hold on
 patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',scar_old_long,'facecolor','interp')
 colorbar

       