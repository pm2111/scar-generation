close all
addpath(genpath('C:/Users/petnov/Dropbox/'));
addpath(genpath('C:/Users/petnov/Dropbox/qrs/'));
addpath(genpath('C:/Users/petnov/Dropbox/sim 3d scripts/'));
addpath(genpath('C:/Users/petnov/Dropbox/shared - copy/mesh scripts/'));

addpath('C:\Users\petnov\Dropbox\sim3d scripts\');
addpath(genpath('C:\Users\petnov\Dropbox\mesh scripts\'));
addpath(genpath('C:\Users\petnov\Dropbox\qrs\'));
addpath(genpath('C:\Users\petnov\Dropbox\tools\'));

addpath('C:\Users\petnov\Dropbox\shared - copy\IO\');


sim_path='D:\ARVC meshing automatic\patients\patient06\results\Eikonal\exp30\';
plt_acts=0;
keep_models={};
patient_nr='06';
sim_name='Eikonal';
study='exp32';
patient_sim_folder =  strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\');
%H_mesh = read_CHASTE( strcat(patient_sim_folder,    'chaste',patient_nr,'\HEART'));
for z=[1]
    if z==1
        pop_path='D:\ARVC meshing automatic\patients\patient06\results\Eikonal_coarse\scar_fiber_dense_RV_post_ECGs\';
    else
        pop_path='D:\ARVC meshing automatic\patients\patient06\results\Eikonal_coarse\eikonal06_scar_border_ecgs\';
    end
end
names = dir(pop_path);
duration=[];
QRSs=[];
ratio=[];
R_peak=[];
    nb=size(names,1);
ecgs={};
upstroke={};
% ecg1=load(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\results\',sim_name,'\',sim_type,'.mat'));
% ecg1=ecg1.ecg_full;

%ecg biomarkers from sims
model_names ={};
R_peak=[];
QRSs=[];
duration=[];
ratio=[];
S_start=[];
for i=3:nb
    idx=[];
    ecg=load([pop_path names(i).name]);
     model_names{i-2}=names(i).name;

    ecg=ecg.vars;
       Ie=ecg(:,1)/10; 
        IIe=ecg(:,2)/10;
        IIIe=ecg(:,3)/10;
        aVRe=ecg(:,4)/10;
        aVLe=ecg(:,5)/10;
        aVFe=ecg(:,6)/10;
        V1e=ecg(:,7)/10;
        V2e=ecg(:,8)/10;
        V3e=ecg(:,9)/10;
        V4e=ecg(:,10)/10;
        V5e=ecg(:,11)/10;
        V6e=ecg(:,12)/10;
      % [Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e]=ecgi_ECG_sim_plot_mat(ecg1,V3e_f,0.1);
    ecgs{i-2}=[Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e];
    full_ecg={};
        full_ecg={Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e};

     %   plot_ecgs_2x6(full_ecg',names(3).name,0,z,1:size(V1e,1));
            hold on
   amplitude_V1(i-2)=max(V1e);
    precordial=[V1e,V2e,V3e,V4e,V5e,V6e];
    for j=1:6
       % [QRSs(j,i-2)]=qrs_estimator(precordial(:,j));
        [~,ratio(j,i-2),~,~,~,~,~,R_peak(j,i-2)]=qrs_estimator(precordial(:,j));
        [ ~,~,duration(j,i-2),swave(j,i-2),QRSs(j,i-2),S_start(j,i-2)]=s_wave_area(precordial(:,j));
    end
   
end



    labels={'slow Purkinje CV','intermediate Purkinje CV','control Purkinje CV'};



k=1;
model_names ={};
%for o=2:nb/2
figure('Position',[0 0 500 1000]); %left bottom width height
z=1;
start=66;
fig_end=72;
for i =start:fig_end%[15:16:256]%[3:nb]%[27,25,47,45]%o+(1:2)%o+122*(0:1)%o+(1:4)%[-5+o+ 5*(1:5)]%[2+o+36*(0:3)]% [-3+o+ 6*(1:6)]%[-5+o+ 5*(1:5)] %[o+(0:3)*25]%[-2+ 5*(1:5)] %3:nb
    idx=[];
    ecg=load([pop_path names(i).name]);
    ecg_control=load('D:\ARVC meshing automatic\patients\patient06\results\Eikonal_coarse\figure_ecgs_fib_lat\17.0x_0f_1.0fiber_roots_RV_ant5.0_RV_post_3.0_RV_septal0.0_NEW_pred_ATMap.csvpseudo_ECG.mat');
     model_names{z}=names(i).name(1:end-4);
scal=30;
    ecg=ecg.vars;
       Ie=ecg(:,1)/scal; 
        IIe=ecg(:,2)/scal;
        IIIe=ecg(:,3)/scal;
        aVRe=ecg(:,4)/scal;
        aVLe=ecg(:,5)/scal;
        aVFe=ecg(:,6)/scal;
        V1e=ecg(:,7)/scal;
        V2e=ecg(:,8)/scal;
        V3e=ecg(:,9)/scal;
        V4e=ecg(:,10)/scal;
        V5e=ecg(:,11)/scal;
        V6e=ecg(:,12)/scal;
      % [Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e]=ecgi_ECG_sim_plot_mat(ecg1,V3e_f,0.1);
    ecgs{i-2}=[Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e];
    full_ecg={};

        full_ecg={Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e};
        for n=1:12
            full_ecg_control{n} = ecg_control.vars(:,n)/scal;
        end
              if i==start
                       plot_ecgs_2x6(full_ecg_control',names(3).name,0,8,1:size(full_ecg_control{1},1));
              end
                        hold on

           plot_ecgs_2x6(full_ecg',names(3).name,0,z,1:size(full_ecg{1},1));


 precordial=[V1e,V2e,V3e,V4e,V5e,V6e];
    for j=1:6
       % [QRSs(j,i-2)]=qrs_estimator(precordial(:,j));
        %[~,ratio(j,i-2),~,~,~,~,~,R_peak(j,i-2)]=qrs_estimator(precordial(:,j));
      %  [ ~,~,duration(j,i-2),swave(j,i-2),QRSs(j,i-2),S_start(j,i-2)]=s_wave_area(precordial(:,j));
    end
   
%    
%     
%     else
%         dummy_name=model_names{k-1};
%         dummy_name(7)='2';
%         A=strfind({names.name},dummy_name);
%         idx=find(not(cellfun('isempty',A)));
z=z+1;
end
%labs={'control and 50% global INa reduction','control', 'slow purkinje and 50% INa reduction','slow Purkinje'};
%labs2={'control','control + LGE', 'slow Purkinje','slow Purkinje + LGE'};
%labs3={'Slow Purkinje Speed','Intermediate Purkinje Speed', 'Control Purkinje Speed'};
%labs4={'no LGE scaling','max LGE scaling ', 'Control Purkinje Speed'};
legend(model_names,'Interpreter', 'none','Fontsize',25);

for n=7:12
    [ ~,~,duration_con(n-6),swave_con(n-6),QRSs_con(n-6),S_start_con(n-6)]=s_wave_area(full_ecg_control{n});
end

labs5={};
labs5{1}='control';

for i=1:fig_end-start+1
    labs5{i+1} = [num2str(2+i), 'mm thick scarred fiber'];
end
legend(labs5,'Interpreter', 'none','Fontsize',25);


labs6={};
labs6{1}='control';
for i=1:fig_end-start+1
    labs6{i+1} = [num2str(2+i), 'mm thick scarred fiber'];
end

%plot the QRS feature for a posterior scarred fiber 
xlabels={'V1','V2','V3','V4','V5','V6'};
cols_colorbar={[0 0 0],[0 0 1], [0.54 0.82 0.99], [0.91 0.41 0.17], [0.32 0.41 0.17],[1 0 0], [1 0 1], [0.55 0 0 ], [0 0 0], [0.55 0 0 ]};
figure()
b=bar(horzcat(QRSs_con', QRSs(1:6,start-2:fig_end-2)),1,'FaceColor','flat')
%xticks([1.5 ,2.5 ,3.5])
xticklabels(xlabels);
ylabel(' QRS duration [ms]')
set(gca,'FontSize',35)
legend(labs6,'Fontsize',35,'Location','northeastoutside')
for i=1:(fig_end-start)+2
b(i).CData=cols_colorbar{i};
end
grid on

%plot the TAD feature for a posterior scarred fiber 
xlabels={'V1','V2','V3','V4'};
cols_colorbar={[0 0 0],[0 0 1], [0.54 0.82 0.99], [0.91 0.41 0.17], [0.32 0.41 0.17],[1 0 0], [1 0 1], [0.55 0 0 ], [0 0 0], [0.55 0 0 ]};
figure()
b=bar(horzcat(duration_con(1:4)', duration(1:4,start-2:fig_end-2)),1,'FaceColor','flat')
%xticks([1.5 ,2.5 ,3.5])
xticklabels(xlabels);
ylabel(' TAD [ms]')
set(gca,'FontSize',35)
legend(labs6,'Fontsize',35,'Location','northeastoutside')
for i=1:(fig_end-start)+2
b(i).CData=cols_colorbar{i};
end

scar=load('D:\ARVC meshing automatic\patients\patient06\scalings\eikonal06_coarse_fib_0_scar_along_fiber_slowing_4.5_length_10_radius_9.csv');
H_mesh2 = read_CHASTE( strcat(patient_sim_folder,    'chaste06_coarse','\HEART'));
figure()
patch('vertices',H_mesh2.xyz,'faces',H_mesh2.face,'edgecolor','none','FaceVertexCData',scar,'facecolor','interp')
headlight
axis off
axis equal

