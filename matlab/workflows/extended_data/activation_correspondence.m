function result = activation_correspondence(configFile)
% Archived projection/Pearson formula for retained ED7e/f; no spatial null.
% Requires original 997-participant HCP activation file and frozen motif t map.
c=jsondecode(fileread(configFile));assert(~isfile(c.output_file),'Use a new output file.');
[mask_L,~,~,~,V_mask_L]=load_cortical_atlas('Yeo2011_7networks','L');
[mask_R,~,~,~,V_mask_R]=load_cortical_atlas('Yeo2011_7networks','R');
mask=[mask_L(:,1:end-1) mask_R];V_mask=[V_mask_L V_mask_R];
L=gifti(template_path('HCP_S1200_Atlas_Z4_pkXDZ/S1200.L.flat.32k_fs_LR.surf.gii'));
R=gifti(template_path('HCP_S1200_Atlas_Z4_pkXDZ/S1200.R.flat.32k_fs_LR.surf.gii'));
L.vertices(:,1)=L.vertices(:,1)-250;R.vertices(:,1)=R.vertices(:,1)+250;
flatVertices=[L.vertices;R.vertices];
activation_map=ft_read_cifti(c.activation_file);
Vec=activation_map.tfmri_wm_2bk_0bk(activation_map.brainstructure<=2);idx=~isnan(Vec);
Map=interpolate_to_grid(Vec(idx),flatVertices(idx,1),flatVertices(idx,2),-500:2:500,-200:2:200,mask);
V_map=conv2(Map,ones(2)/4,'valid');
S=load(c.motif_t_file,'t_values');t_values=S.t_values;
if ndims(t_values)==4
    assert(isfield(c,'motif_contrast_index'),'A 4-D cache requires an explicit motif_contrast_index matching 2-back minus 0-back.');
    validateattributes(c.motif_contrast_index,{'numeric'},{'scalar','integer','>=',1,'<=',size(t_values,4)});
    t_values=t_values(:,:,:,c.motif_contrast_index);
end
assert(ndims(t_values)==3 && size(t_values,2)==4 && size(t_values,3)==3,'Expected preselected vertices x 4 motifs x 3 methods.');
assert(size(t_values,1)==nnz(V_mask),'Original motif-mask ordering must match.');
R_with_activation=nan(size(t_values,3),size(t_values,2));
for i_method=1:size(t_values,3)
    for i_pattern=1:size(t_values,2)
        R_with_activation(i_method,i_pattern)=corr(surrogate_grid_to_vector(V_map,V_mask),t_values(:,i_pattern,i_method));
    end
end
result=struct('activation_native',Map,'activation_vvf',V_map,'R_with_activation',R_with_activation, ...
    'scope','Recovered projection plus archived commented Pearson expression; original activation input and exact historical correspondence unresolved');
save(c.output_file,'result','-v7');
end
