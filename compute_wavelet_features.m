function [Out]=compute_wavelet_features(fw,vring,vertex)
n=length(vertex{1}); %number of the vertex in the coarse mesh
Out.bf_coarse=fw(1:n,:);

Out.bf_ave_coarse=zeros(n,3);
Out.bf_var_coarse=zeros(n,3);
Out.bf_coarse_angle_mean=zeros(n,1); 
Out.bf_coarse_angle_var=zeros(n,1);
Out.bf_coarse_norm_mean=zeros(n,1);
Out.bf_coarse_norm_var=zeros(n,1);

Out.bf_med_coarse=zeros(n,3);
Out.bf_coarse_med_angle_mean=zeros(n,1);
Out.bf_coarse_med_angle_var=zeros(n,1);
Out.bf_coarse_med_norm_mean=zeros(n,1);
Out.bf_coarse_med_norm_var=zeros(n,1);
% pca_wcv=zeros(n,3);
% pca_wcv_evr=zeros(n,1);
for i=1:n
    Out.bf_ave_coarse(i,:)=mean(fw(vring{i},:)); %each vertex in the coarse mesh has an vector which is average of the wcv of the one ring neighbors
    Out.bf_med_coarse(i,:)=median(fw(vring{i},:)); 
    Out.bf_var_coarse(i,:)= var([fw(vring{i},:);fw(i,:)]);
    bf_coarse_angle=vectorAngle3d(Out.bf_ave_coarse(i,:),fw(vring{i},:));% the angles between the averg wcv and its 1 ring wcv
    bf_coarse_angle(isnan(bf_coarse_angle))=eps; 
    Out.bf_coarse_angle_mean(i,:)=mean(bf_coarse_angle);
    Out.bf_coarse_angle_var(i,:)=var(bf_coarse_angle);
    
%     bf_coarse_med_angle=vectorAngle3d(Out.bf_med_coarse(i,:),fw(vring{i},:));
%     bf_coarse_med_angle(isnan(bf_coarse_med_angle))=eps;
%     Out.bf_coarse_med_angle_mean(i,:)=mean(bf_coarse_med_angle);
%     Out.bf_coarse_med_angle_var(i,:)=var(bf_coarse_med_angle);
%     wcv_angle_coarse(isnan(wcv_angle_coarse))=eps;
    bf_coarse_norm=abs(vectorNorm3d(Out.bf_ave_coarse(i,:))-vectorNorm3d(fw(vring{i},:))); % the difference between the norm of the averaged wcv and the norm of its 1 ring wcv
    Out.bf_coarse_norm_mean(i,:)=mean(bf_coarse_norm);
    Out.bf_coarse_norm_var(i,:)=var(bf_coarse_norm);
    
%     bf_coarse_med_norm=abs(vectorNorm3d(Out.bf_med_coarse(i,:))-vectorNorm3d(fw(vring{i},:)));
%     Out.bf_coarse_med_norm_mean(i,:)=mean(bf_coarse_med_norm);
%     Out.bf_coarse_med_norm_var(i,:)=var(bf_coarse_med_norm);
%     [coeff,score,latent]=pca(fw(vring{i},:));
%     pca_wcv(i,:)=coeff(:,1)';
%     pca_wcv_evr(i,:)=latent(1)/latent(2);
end

m=length(vertex{2});
Out.wcv_ave_fine=zeros(m-n,3);
Out.wcv_angle_mean=zeros(m-n,1);
Out.wcv_angle_var=zeros(m-n,1);
Out.wcv_var=zeros(m-n,3);
Out.wcv_fine_norm_mean=zeros(m-n,1);
Out.wcv_fine_norm_var=zeros(m-n,1);
% wcv_pca_fine=zeros(m-n,3);
% wcv_pca_fine_evr=zeros(m-n,1);

Out.wcv_med_fine=zeros(m-n,3);
Out.wcv_med_angle_mean=zeros(m-n,1);
Out.wcv_med_angle_var=zeros(m-n,1);
Out.wcv_fine_med_norm_mean=zeros(m-n,1);
Out.wcv_fine_med_norm_var=zeros(m-n,1);

for j=n+1:m
    fw_ave=mean(fw(vring{j}(vring{j}>n),:));
    Out.wcv_ave_fine(j-n,:)=fw_ave; %average of the wc vectors of the 1ring neighbors
%     Out.wcv_med_fine(j-n,:)=median(fw(vring{j}(vring{j}>n),:));
%     [fw_pca,score,latent]=pca(fw(vring{j}(vring{j}>n),:));
%     wcv_pca_fine(j-n,:)=fw_pca(:,1)'; %first principal component of the wc vectors of the 1 ring neighbors
%     wcv_pca_fine_evr(j-n,:)=latent(1)/latent(2);
%     o=length(vring{j}(vring{j}>n));
%     fw_ring=vring{j}(vring{j}>n);
%     angle=zeros(o,1);
    Out.wcv_var(j-n,:)=var([fw(vring{j}(vring{j}>n),:);fw(j,:)]);
    wcv_angle=vectorAngle3d(fw(j,:),fw(vring{j}(vring{j}>n),:)); % the angle between the wcv in fine mesh and its 1 ring wcv        
    wcv_angle(isnan(wcv_angle))=eps;
%     wcv_med_angle=vectorAngle3d(Out.wcv_med_fine(j-n,:),fw(vring{j}(vring{j}>n),:));
%     wcv_med_angle(isnan(wcv_med_angle))=eps;
%     Out.wcv_med_angle_mean(j-n)=mean(wcv_med_angle);
%     Out.wcv_med_angle_var(j-n)=var(wcv_med_angle);
%     wcv_fine_med_norm=abs(vectorNorm3d(Out.wcv_med_fine(j-n,:))-vectorNorm3d(fw(vring{j}(vring{j}>n),:)));
%     Out.wcv_fine_med_norm_mean(j-n)=mean(wcv_fine_med_norm);
%     Out.wcv_fine_med_norm_var(j-n)=var(wcv_fine_med_norm);
    
    Out.wcv_angle_mean(j-n)=mean(wcv_angle);%average angle between the wcv and its 1 ring neighbors
    Out.wcv_angle_var(j-n)=var(wcv_angle);
    wcv_fine_norm=abs(vectorNorm3d(fw(j,:))-vectorNorm3d(fw(vring{j}(vring{j}>n),:))); %the difference between the norm of the fine wcv and the norm of its 1 ring wcv
    Out.wcv_fine_norm_mean(j-n)=mean(wcv_fine_norm);
    Out.wcv_fine_norm_var(j-n)=var(wcv_fine_norm);
%     index=vring{j}(vring{j}<n+1);
%     edges_angle(j-n)=vectorAngle3d((vertex{2}(:,index(1))-vertex{2}(:,j))',(vertex{2}(:,index(2))-vertex{2}(:,j))');
end
