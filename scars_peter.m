%% load up the mesh (takes a while)

addpath(genpath('C:\Users\petnov\Dropbox\shared - copy\'));
addpath(genpath('C:\Users\petnov\Dropbox\mesh scripts'))
addpath('C:\Users\petnov\Dropbox\shared - copy\MESHES\');
addpath('C:\Users\petnov\Dropbox\shared - copy\IO\');
addpath(genpath('C:\Users\petnov\Dropbox\scar generation\'))
enableVTK;
pats={'06'};
patient_nr=pats{1};
    %mkdir 'D:\ARVC meshing automatic\patients\patient01\regions';

patient_sim_folder =  strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\');
H_mesh = read_CHASTE( strcat(patient_sim_folder,    'chaste06_coarse','\HEART'));

tet = double(H_mesh.tri);
pts = double(H_mesh.xyz);

%% compute connectivity matrix (also takes a while)

p = nchoosek(1:4,2);
n = double(max(max(tet(:,2:end))));

C = sparse(n,n);
tic
for i = 1 : size(p,1)
    j = double(tet(:,p(i,:)));
    d = sum((pts(j(:,1),:)-pts(j(:,2),:)).^2,2).^.5;
    C = C+sparse([j(:,1); j(:,2)],[j(:,2); j(:,1)],[d; d],n,n);
end
C_time = toc

%% visualize heart (yes, this also takes a while)

tic
tr = TriRep(tet,pts);
tri = tr.freeBoundary;
boundary_time = toc
trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'linestyle','none');
headlight;
%% Select a number of datatips to mark the general area for the scar on the cardiac surface (e.g. the LAD)

% No code here, this is all manual labour! 
% Hold down alt to select multiple datatips.
% Once you're done, export the datatips to the matlab workspace as a
% structure named cursor_info (matlab defaults to this name).

%% Generate candidate nodes for the scar

% the cutoff is z-dependent, generating more candidates towards the apex
% to account for the coronary tree spreading in base->apex direction

cutoff = (1.2-(pts(:,3)-min(pts(:,3)))./(max(pts(:,3))-min(pts(:,3))) )*(max(pts(:,3))-min(pts(:,3)))*1.2;

clear ind
for i = 1 : length(cursor_info)
    p = cursor_info(i).Position;
    [~,ind(i)] = min(sqrt(sum([pts(:,1)-p(1) pts(:,2)-p(2) pts(:,3)-p(3)].^2,2)));
end
d = dijkstra2(C,sort(ind),ind*0);
scar = find(d'<cutoff);

expr = pts(:,1)*0;
expr(scar) = 1;

% visualize to check that all looks ok
trisurf(tri,pts(:,1),pts(:,2),pts(:,3),expr,'facecolor','interp','linestyle','none','facelighting','phong','AmbientStrength',.8,'DiffuseStrength',.5)

%% Generate core and border zones

% scar parameters - these might need playing around with to get what you want
ncores = 1;
core_ext_noise = 0;
core_post_noise = 0.05;
cz_threshold = .2;
bz_threshold = .25;

core = scar(randi(length(scar),ncores,1));

d = dijkstra2(C,sort(core),abs(randn(size(core,1),1))*core_ext_noise);
d = d/max(d)+randn(1,size(d,2))*core_post_noise;

% Final scar definition
cz_nodes = d<cz_threshold;
bz_nodes = d<bz_threshold;
cz_tets = cz_nodes(tet(:,1)) | cz_nodes(tet(:,2)) | cz_nodes(tet(:,3)) | cz_nodes(tet(:,4));
bz_tets = bz_nodes(tet(:,1)) | bz_nodes(tet(:,2)) | bz_nodes(tet(:,3)) | bz_nodes(tet(:,4));

% gradient between 0 in the CZ and 1 in the healthy tissue
bz_grad = d;
bz_grad(d<=cz_threshold) = 0;
bz_grad(d>cz_threshold & d<=bz_threshold) = (bz_grad(d>cz_threshold & d<bz_threshold)-cz_threshold)/(bz_threshold-cz_threshold);
bz_grad(d>bz_threshold) = 1;

%% visualize final result

trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor','w','linestyle','none','facealpha',.1);

% scar core
tr = TriRep(tet(cz_tets,:),pts);
cz_tri = tr.freeBoundary();
hold on
trisurf(cz_tri,pts(:,1),pts(:,2),pts(:,3),'linestyle','none','facecolor',[1 0 0],'ambientstrength',.5,'diffusestrength',.5,'facelighting','phong','facealpha',1);
hold off
light

% scar bz
tr = TriRep(tet(bz_tets,:),pts);
bz_tri = tr.freeBoundary();
hold on
trisurf(bz_tri,pts(:,1),pts(:,2),pts(:,3),'linestyle','none','facecolor',[0 1 0],'ambientstrength',.5,'diffusestrength',.5,'facelighting','phong','facealpha',.1);
hold off

