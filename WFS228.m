function [ F ] = WFS228( ori_file,dt,Tmax )
%Input:
% ori_file - the 3D *.off file to be investigated
% dt and Tmax - the parameter used to control the strength and number of iteraton of Laplacian smoothing
%Output:
% F - 1x228 Feature vector
%
% Contact: zheenyuli@gmail.com
%
% References:
% Li Z, Bors A G. Steganalysis of Meshes Based on 3D Wavelet multiresolution Analysis. 
% Information Sciences, 2020,522:164-179.
%
% Zhenyu Li 2020-March

addpath('functions');

[vertex_o, face_o] = read_off(ori_file);

%%
%smoothing
options.symmetrize = 0;
options.normalize = 1;

% laplacian_type = 'distance';
laplacian_type = 'combinatorial';
% laplacian_type = 'conformal';

disp('--> Computing laplacian');
L = compute_mesh_laplacian(vertex_o,face_o,laplacian_type,options);

options.dt = dt;
options.Tmax = Tmax;
vertex_s = perform_mesh_heat_diffusion(vertex_o,face_o,L,options);
face_s=face_o;

vertex_o=preprocess(vertex_o);
vertex_s=preprocess(vertex_s);


global e2f; %e2f(i,j) and e2f(j,i) are the number of the two faces adjacent to edge (i,j).

e2f = compute_edge_face_ring(face_o);

%%
%initial resolution features

edge_o=meshEdges(face_o');
edge_s=edge_o;

crossvector_o=zeros(3,length(edge_o));
crossvector_s=zeros(3,length(edge_o));
for i=1:length(edge_o)
    %get the face index adjacent to the edge_i
    face_1=e2f(edge_o(i,1),edge_o(i,2));
    face_2=e2f(edge_o(i,2),edge_o(i,1));
    %get the vertex index in the two adjacent faces
    v_f1=my_setdiff(face_o(:,face_1),edge_o(i,:));
    v_f2=my_setdiff(face_o(:,face_2),edge_o(i,:));
    %vector form v_f1 to v_f2
    crossvector_o(:,i)=vertex_o(:,v_f2)-vertex_o(:,v_f1);
    crossvector_s(:,i)=vertex_s(:,v_f2)-vertex_s(:,v_f1);
end

edge_vector_o=vertex_o(:,(edge_o(:,1)))-vertex_o(:,(edge_o(:,2)));
edge_vector_s=vertex_s(:,(edge_s(:,1)))-vertex_s(:,(edge_s(:,2)));

edgediff=abs(edge_vector_o-edge_vector_s)+eps;
edgediffN=vectorNorm3d(edgediff');
edgeAngle=vectorAngle3d(edge_vector_o',edge_vector_s')+eps;
edgeNormDiff=abs(vectorNorm3d(edge_vector_o')-vectorNorm3d(edge_vector_s'))+eps;

crossvector_diff=abs(crossvector_o'-crossvector_s')+eps;
crossvector_diffN=vectorNorm3d(crossvector_diff');
crossvector_angle=vectorAngle3d(crossvector_o',crossvector_s')+eps;
crossvector_normdiff=abs(vectorNorm3d(crossvector_o')-vectorNorm3d(crossvector_s'))+eps;

F=zeros(1,228);

F(1,1:12) = [mean(log(edgediff')),var(log(edgediff')),...
    skewness(log(edgediff')), kurtosis(log(edgediff'))];
F(1,13:16) = [mean(log(edgediffN)),var(log(edgediffN)),...
    skewness(log(edgediffN)), kurtosis(log(edgediffN))];
F(1,17:20) = [mean(log(edgeAngle)),var(log(edgeAngle)),...
    skewness(log(edgeAngle)), kurtosis(log(edgeAngle))];
F(1,21:24) = [mean(log(edgeNormDiff)),var(log(edgeNormDiff)),...
    skewness(log(edgeNormDiff)), kurtosis(log(edgeNormDiff))];
F(1,25:36) = [mean(log(crossvector_diff)),var(log(crossvector_diff)),...
    skewness(log(crossvector_diff)), kurtosis(log(crossvector_diff))];
F(1,37:40) = [mean(log(crossvector_angle)),var(log(crossvector_angle)),...
    skewness(log(crossvector_angle)), kurtosis(log(crossvector_angle))];
F(1,41:44) = [mean(log(crossvector_normdiff)),var(log(crossvector_normdiff)),...
    skewness(log(crossvector_normdiff)), kurtosis(log(crossvector_normdiff))];
F(1,45:48) = [mean(log(crossvector_diffN)),var(log(crossvector_diffN)),...
    skewness(log(crossvector_diffN)), kurtosis(log(crossvector_diffN))];

%%
%lower resolution features

face_num=length(face_o);
fea_o=zeros(6,face_num*3);
fea_s=fea_o;
c_edge_o=fea_o;
c_edge_s=fea_o;
for i=1:face_num % for each face in the mesh
    
    %get the vertex in the face(i) 
    v1=face_o(1,i); 
    v2=face_o(2,i);
    v3=face_o(3,i);
    
    %get the index of faces adjacent to the face(i) 
    if e2f(v1,v2)==i
        face_1=e2f(v2,v1);
    else
        face_1=e2f(v1,v2);
    end
    
    if e2f(v1,v3)==i
        face_2=e2f(v3,v1);
    else
        face_2=e2f(v1,v3);
    end
    
    if e2f(v3,v2)==i
        face_3=e2f(v2,v3);
    else
        face_3=e2f(v3,v2);
    end
    
    %get the vertex in the coarsed face
    v_f1=my_setdiff(face_o(:,face_1),face_o(:,i));
    v_f2=my_setdiff(face_o(:,face_2),face_o(:,i));
    v_f3=my_setdiff(face_o(:,face_3),face_o(:,i));
    
    if  ~isempty(v_f1) && ~isempty(v_f2) && ~isempty(v_f3)
    if abs((v_f1-v_f2)*(v_f1-v_f3)*(v_f2-v_f3))>0     
    
    %the vector from the vertex to its corresponding coarsed edge's middle point
    vector_v1=0.5*(vertex_o(:,v_f1)+vertex_o(:,v_f2))-vertex_o(:,v1);    
    vector_v2=0.5*(vertex_o(:,v_f1)+vertex_o(:,v_f3))-vertex_o(:,v2);
    vector_v3=0.5*(vertex_o(:,v_f3)+vertex_o(:,v_f2))-vertex_o(:,v3);
    
    vector_v1_s=0.5*(vertex_s(:,v_f1)+vertex_s(:,v_f2))-vertex_s(:,v1);    
    vector_v2_s=0.5*(vertex_s(:,v_f1)+vertex_s(:,v_f3))-vertex_s(:,v2);
    vector_v3_s=0.5*(vertex_s(:,v_f3)+vertex_s(:,v_f2))-vertex_s(:,v3);
    
    %get the coarsed edge
    coarsed_edge_1_o=vertex_o(:,v_f1)-vertex_o(:,v_f2);
    coarsed_edge_2_o=vertex_o(:,v_f1)-vertex_o(:,v_f3);
    coarsed_edge_3_o=vertex_o(:,v_f3)-vertex_o(:,v_f2);
    
    coarsed_edge_1_s=vertex_s(:,v_f1)-vertex_s(:,v_f2);
    coarsed_edge_2_s=vertex_s(:,v_f1)-vertex_s(:,v_f3);
    coarsed_edge_3_s=vertex_s(:,v_f3)-vertex_s(:,v_f2);
    
    %list the vector
    if v_f1<v_f2
        fea_o(:,i*3-2)=[v1;v_f1;v_f2;vector_v1];
        fea_s(:,i*3-2)=[v1;v_f1;v_f2;vector_v1_s];
        c_edge_o(:,i*3-2)=[v1;v_f1;v_f2;coarsed_edge_1_o];
        c_edge_s(:,i*3-2)=[v1;v_f1;v_f2;coarsed_edge_1_s];
    else
        fea_o(:,i*3-2)=[v1;v_f2;v_f1;vector_v1];
        fea_s(:,i*3-2)=[v1;v_f2;v_f1;vector_v1_s];
        c_edge_o(:,i*3-2)=[v1;v_f2;v_f1;coarsed_edge_1_o];
        c_edge_s(:,i*3-2)=[v1;v_f2;v_f1;coarsed_edge_1_s];
    end
    if v_f1<v_f3
        fea_o(:,i*3-1)=[ v2;v_f1;v_f3;vector_v2];
        fea_s(:,i*3-1)=[ v2;v_f1;v_f3;vector_v2_s];
        c_edge_o(:,i*3-1)=[v1;v_f1;v_f3;coarsed_edge_2_o];
        c_edge_s(:,i*3-1)=[v1;v_f1;v_f3;coarsed_edge_2_s];
    else
        fea_o(:,i*3-1)=[ v2;v_f3;v_f1;vector_v2];
        fea_s(:,i*3-1)=[ v2;v_f3;v_f1;vector_v2_s];
        c_edge_o(:,i*3-1)=[v1;v_f3;v_f1;coarsed_edge_2_o];
        c_edge_s(:,i*3-1)=[v1;v_f3;v_f1;coarsed_edge_2_s];
    end
    if v_f2<v_f3
        fea_o(:,i*3)=[ v3;v_f2;v_f3;vector_v3];
        fea_s(:,i*3)=[ v3;v_f2;v_f3;vector_v3_s];
        c_edge_o(:,i*3)=[v1;v_f2;v_f3;coarsed_edge_3_o];
        c_edge_s(:,i*3)=[v1;v_f2;v_f3;coarsed_edge_3_s];
    else
        fea_o(:,i*3)=[ v3;v_f3;v_f2;vector_v3];
        fea_s(:,i*3)=[ v3;v_f3;v_f2;vector_v3_s];
        c_edge_o(:,i*3)=[v1;v_f3;v_f2;coarsed_edge_3_o];
        c_edge_s(:,i*3)=[v1;v_f3;v_f2;coarsed_edge_3_s];
    end
    end
    end
   
end
[fea_o,ia,ic] = unique(fea_o','rows');
fea_s=fea_s(:,ia)';
c_edge_o=c_edge_o(:,ia)';
c_edge_s=c_edge_s(:,ia)';


fw_o=fea_o(:,4:6);
fw_s=fea_s(:,4:6);
c_edge_o=c_edge_o(:,4:6)+eps;
c_edge_s=c_edge_s(:,4:6)+eps;

Diff_Angle_w_edge=abs(vectorAngle3d(fw_o,c_edge_o)-vectorAngle3d(fw_s,c_edge_s))+eps;
Diff_NRatio_w_edge=abs(vectorNorm3d(fw_o)./vectorNorm3d(c_edge_o)-vectorNorm3d(fw_s)./vectorNorm3d(c_edge_s))+eps;
Diff_NRatio_w_edge(isnan(Diff_NRatio_w_edge))=eps;

Diff_w=abs(fw_o-fw_s)+eps; 
NDiff_w=abs(vectorNorm3d(fw_o)-vectorNorm3d(fw_s))+eps;
DiffN_w=vectorNorm3d(Diff_w);

Angle_w=vectorAngle3d(fw_o,fw_s)+eps;
Angle_w(isnan(Angle_w))=eps;

Diff_cedge=abs(c_edge_o-c_edge_s)+eps; 
NDiff_cedge=abs(vectorNorm3d(c_edge_o)-vectorNorm3d(c_edge_s))+eps;
DiffN_cedge=vectorNorm3d(Diff_cedge);

Angle_cedge=vectorAngle3d(c_edge_o,c_edge_s)+eps;
Angle_cedge(isnan(Angle_cedge))=eps;

 
F(1,49:60) = [mean(log(Diff_w)),var(log(Diff_w)),...
    skewness(log(Diff_w)), kurtosis(log(Diff_w))];

F(1,61:64) = [mean(log(NDiff_w)),var(log(NDiff_w)),...
    skewness(log(NDiff_w)), kurtosis(log(NDiff_w))];

F(1,65:68) = [mean(log(Angle_w)),var(log(Angle_w)),...
    skewness(log(Angle_w)), kurtosis(log(Angle_w))];

F(1,69:80) = [mean(log(Diff_cedge)),var(log(Diff_cedge)),...
    skewness(log(Diff_cedge)), kurtosis(log(Diff_cedge))];

 
F(1,81:84) = [mean(log(NDiff_cedge)),var(log(NDiff_cedge)),...
    skewness(log(NDiff_cedge)), kurtosis(log(NDiff_cedge))];

F(1,85:88) = [mean(log(Angle_cedge)),var(log(Angle_cedge)),...
    skewness(log(Angle_cedge)), kurtosis(log(Angle_cedge))];

F(1,89:92) = [mean(log(DiffN_w)),var(log(DiffN_w)),...
    skewness(log(DiffN_w)), kurtosis(log(DiffN_w))];

F(1,93:96) = [mean(log(DiffN_cedge)),var(log(DiffN_cedge)),...
    skewness(log(DiffN_cedge)), kurtosis(log(DiffN_cedge))];

F(1,97:100) = [mean(log(Diff_Angle_w_edge)),var(log(Diff_Angle_w_edge)),...
    skewness(log(Diff_Angle_w_edge)), kurtosis(log(Diff_Angle_w_edge))];

F(1,101:104) = [mean(log(Diff_NRatio_w_edge)),var(log(Diff_NRatio_w_edge)),...
    skewness(log(Diff_NRatio_w_edge)), kurtosis(log(Diff_NRatio_w_edge))];


%%
%higher resolution features

%generate a fine mesh for original mesh
vertex = {}; face = {};

vertex{1}=vertex_o;
face{1}=face_o;

options.name = ori_file;
% you can also try with 'sqrt3', 'butterfly', 'linear4'
options.sub_type = 'linear4';
% options.sub_type = 'butterfly';

options.spherical = 0;
options.verb = 0;
% clf;
for j=2:2
    if j>1
        [vertex{j},face{j}] = perform_mesh_subdivision(vertex{j-1}, face{j-1}, 1, options);
    end
    
end


f = vertex{end}';
% forward wavelet tranform
edge_vector_o=zeros(size(f));
[fw_o,edge_vector_o,vring_o] = perform_wavelet_mesh_transform(vertex,face, f, edge_vector_o,+1, options);

%vertex_bf_o is the vertex using the butterfly frame
vertex_bf_o=vertex{2}';
vertex_bf_o(1:length(vertex_o),:)=vertex_o';
vertex_bf_o(length(vertex_o)+1:end,:)=vertex_bf_o(length(vertex_o)+1:end,:)-fw_o(length(vertex_o)+1:end,:);
face_bf_o=face{2};
edge_bf_o=meshEdges(face_bf_o);
[Fout_o]=compute_wavelet_features(fw_o,vring_o,vertex);

%%
%generate a fine mesh for smoothed mesh
vertex = {}; face = {};
% [vertex{1},face{1}] = read_off(smo_file);
vertex{1}=vertex_s;
face{1}=face_s;
vertex{1}=preprocess(vertex{1});
options.name = ori_file;
% you can also try with 'sqrt3', 'butterfly', 'linear4'
options.sub_type = 'linear4';
% options.sub_type = 'butterfly';

options.spherical = 0;
options.verb = 0;

for j=2:2
    if j>1
        [vertex{j},face{j}] = perform_mesh_subdivision(vertex{j-1}, face{j-1}, 1, options);
    end
    
end

f = vertex{end}';
% forward wavelet tranform
edge_vector_s=zeros(size(f));
[fw_s,edge_vector_s,vring_s] = perform_wavelet_mesh_transform(vertex,face, f, edge_vector_s,+1, options);
vertex_bf_s=vertex{2}';
vertex_bf_s(1:length(vertex_s),:)=vertex_s';
vertex_bf_s(length(vertex_s)+1:end,:)=vertex_bf_s(length(vertex_s)+1:end,:)-fw_s(length(vertex_s)+1:end,:);
face_bf_s=face{2};
edge_bf_s=meshEdges(face_bf_s);

[Fout_s]=compute_wavelet_features(fw_s,vring_s,vertex);


%%
%forming the features
Num_vertex_coarse=length(vertex{1});

fw_o=fw_o(1+Num_vertex_coarse:end,:);
fw_s=fw_s(1+Num_vertex_coarse:end,:);
edge_vector_o=edge_vector_o(1+Num_vertex_coarse:end,:)+eps;
edge_vector_s=edge_vector_s(1+Num_vertex_coarse:end,:)+eps;

Diff_Angle_w_edge=abs(vectorAngle3d(fw_o,edge_vector_o)-vectorAngle3d(fw_s,edge_vector_s))+eps;
Diff_NRatio_w_edge=abs(vectorNorm3d(fw_o)./vectorNorm3d(edge_vector_o)-vectorNorm3d(fw_s)./vectorNorm3d(edge_vector_s))+eps;
Diff_NRatio_w_edge(isnan(Diff_NRatio_w_edge))=eps;

Diff_w=abs(fw_o-fw_s)+eps;  %difference between the wc vectors of the fine (original and smooth) mesh
Diff_waf=abs(Fout_o.wcv_ave_fine-Fout_s.wcv_ave_fine)+eps;  %difference between the averaged wc vectors of the fine (original and smooth) mesh
Diff_w_f=abs((fw_o-Fout_o.wcv_ave_fine)-(fw_s-Fout_s.wcv_ave_fine))+eps; %wc vector - average wc vector in the fine mesh 
 

DiffN_w=vectorNorm3d(Diff_w);
 
DiffN_waf=vectorNorm3d(Diff_waf);
DiffN_w_f=vectorNorm3d(Diff_w_f);
 

NDiff_w=abs(vectorNorm3d(fw_o)-vectorNorm3d(fw_s))+eps;
 
NDiff_waf=abs(vectorNorm3d(Fout_o.wcv_ave_fine)-vectorNorm3d(Fout_s.wcv_ave_fine))+eps;
NDiff_w_f=abs(vectorNorm3d(fw_o-Fout_o.wcv_ave_fine)-vectorNorm3d(fw_s-Fout_s.wcv_ave_fine))+eps;
 


Angle_w=vectorAngle3d(fw_o,fw_s)+eps;
Angle_w(isnan(Angle_w))=eps;
 
Angle_waf=vectorAngle3d(Fout_o.wcv_ave_fine,Fout_s.wcv_ave_fine)+eps;
Angle_waf(isnan(Angle_waf))=eps;
Angle_w_f=vectorAngle3d((fw_o-Fout_o.wcv_ave_fine),(fw_s-Fout_s.wcv_ave_fine))+eps;
Angle_w_f(isnan(Angle_w_f))=eps;
 

Diff_angle_ave_fine=abs(vectorAngle3d(fw_o,Fout_o.wcv_ave_fine)-vectorAngle3d(fw_s,Fout_s.wcv_ave_fine))+eps;

Diff_wcv_angle_mean=abs(Fout_o.wcv_angle_mean-Fout_s.wcv_angle_mean)+eps;
Diff_wcv_angle_var=abs(Fout_o.wcv_angle_var-Fout_s.wcv_angle_var)+eps;

Diff_wcv_norm_mean=abs(Fout_o.wcv_fine_norm_mean-Fout_s.wcv_fine_norm_mean)+eps;
Diff_wcv_norm_var=abs(Fout_o.wcv_fine_norm_var-Fout_s.wcv_fine_norm_var)+eps;



edge_vector_bf_o=vertex_bf_o(edge_bf_o(:,1),:)-vertex_bf_o(edge_bf_o(:,2),:);
edge_vector_bf_s=vertex_bf_s(edge_bf_s(:,1),:)-vertex_bf_s(edge_bf_s(:,2),:);

edgediff=abs(edge_vector_bf_o-edge_vector_bf_s)+eps;
edgediffN=vectorNorm3d(edgediff);
edgeNormDiff=abs(vectorNorm3d(edge_vector_bf_o)-vectorNorm3d(edge_vector_bf_s))+eps;
edgeAngle=vectorAngle3d(edge_vector_bf_o,edge_vector_bf_s)+eps;

%Calculate the stastics of the vertors

F(1,105:116) = [mean(log(Diff_w)),var(log(Diff_w)),...
    skewness(log(Diff_w)), kurtosis(log(Diff_w))];
F(1,117:120) = [mean(log(NDiff_w)),var(log(NDiff_w)),...
    skewness(log(NDiff_w)), kurtosis(log(NDiff_w))];
F(1,121:124) = [mean(log(Angle_w)),var(log(Angle_w)),...
    skewness(log(Angle_w)), kurtosis(log(Angle_w))];
F(1,125:128) = [mean(log(DiffN_w)),var(log(DiffN_w)),...
    skewness(log(DiffN_w)), kurtosis(log(DiffN_w))];

F(1,129:140) = [mean(log(Diff_waf)),var(log(Diff_waf)),...
    skewness(log(Diff_waf)), kurtosis(log(Diff_waf))];
F(1,141:144) = [mean(log(NDiff_waf)),var(log(NDiff_waf)),...
    skewness(log(NDiff_waf)), kurtosis(log(NDiff_waf))];
F(1,145:148) = [mean(log(Angle_waf)),var(log(Angle_waf)),...
    skewness(log(Angle_waf)), kurtosis(log(Angle_waf))];
F(1,149:152) = [mean(log(DiffN_waf)),var(log(DiffN_waf)),...
    skewness(log(DiffN_waf)), kurtosis(log(DiffN_waf))];


F(1,153:164) = [mean(log(Diff_w_f)),var(log(Diff_w_f)),...
    skewness(log(Diff_w_f)), kurtosis(log(Diff_w_f))];
F(1,165:168) = [mean(log(NDiff_w_f)),var(log(NDiff_w_f)),...
    skewness(log(NDiff_w_f)), kurtosis(log(NDiff_w_f))];
F(1,169:172) = [mean(log(Angle_w_f)),var(log(Angle_w_f)),...
    skewness(log(Angle_w_f)), kurtosis(log(Angle_w_f))];
F(1,173:176) = [mean(log(DiffN_w_f)),var(log(DiffN_w_f)),...
    skewness(log(DiffN_w_f)), kurtosis(log(DiffN_w_f))];


F(1,177:180) = [mean(log(Diff_angle_ave_fine)),var(log(Diff_angle_ave_fine)),...
    skewness(log(Diff_angle_ave_fine)), kurtosis(log(Diff_angle_ave_fine))];
F(1,181:184) = [mean(log(Diff_wcv_angle_mean)),var(log(Diff_wcv_angle_mean)),...
    skewness(log(Diff_wcv_angle_mean)), kurtosis(log(Diff_wcv_angle_mean))];
F(1,185:188) = [mean(log(Diff_wcv_angle_var)),var(log(Diff_wcv_angle_var)),...
    skewness(log(Diff_wcv_angle_var)), kurtosis(log(Diff_wcv_angle_var))];

F(1,189:192) = [mean(log(Diff_wcv_norm_mean)),var(log(Diff_wcv_norm_mean)),...
    skewness(log(Diff_wcv_norm_mean)), kurtosis(log(Diff_wcv_norm_mean))];
F(1,193:196) = [mean(log(Diff_wcv_norm_var)),var(log(Diff_wcv_norm_var)),...
    skewness(log(Diff_wcv_norm_var)), kurtosis(log(Diff_wcv_norm_var))];


F(1,197:208) = [mean(log(edgediff)),var(log(edgediff)),...
    skewness(log(edgediff)), kurtosis(log(edgediff))];
F(1,209:212) = [mean(log(edgediffN)),var(log(edgediffN)),...
    skewness(log(edgediffN)), kurtosis(log(edgediffN))];
F(1,213:216) = [mean(log(edgeNormDiff)),var(log(edgeNormDiff)),...
    skewness(log(edgeNormDiff)), kurtosis(log(edgeNormDiff))];
F(1,217:220) = [mean(log(edgeAngle)),var(log(edgeAngle)),...
    skewness(log(edgeAngle)), kurtosis(log(edgeAngle))];

F(1,221:224) = [mean(log(Diff_Angle_w_edge)),var(log(Diff_Angle_w_edge)),...
    skewness(log(Diff_Angle_w_edge)), kurtosis(log(Diff_Angle_w_edge))];

F(1,225:228) = [mean(log(Diff_NRatio_w_edge)),var(log(Diff_NRatio_w_edge)),...
    skewness(log(Diff_NRatio_w_edge)), kurtosis(log(Diff_NRatio_w_edge))];


end

