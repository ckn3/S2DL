function super_label=seg_ERS(HSI,L_base, NumSP)

[no_rows,no_lines, no_bands] = size(HSI);
img=reshape(HSI,[no_rows*no_lines,no_bands]);
p=3;           
[~,score,~] = pca(img);
PCA_data = score(:,1:p);
base_image=reshape(PCA_data, no_rows,no_lines, p);
for i = 1: p
    base_image(:,:,i)=mat2gray(base_image(:,:,i));        
    base_image(:,:,i)=im2uint8(base_image(:,:,i));        
end
BW(:,:,1) = edge(base_image(:,:,1),'Sobel');
BW(:,:,2) = edge(base_image(:,:,2),'Sobel');
BW(:,:,3) = edge(base_image(:,:,3),'Sobel');

BW = sum(BW,3);
Rtexture = size(find(BW~=0),1)/(no_rows*no_lines);
if L_base == 0
    L_superpixel=NumSP;
else
    L_superpixel=round(Rtexture*L_base);
end
super_label = mex_ers(double(base_image),L_superpixel);