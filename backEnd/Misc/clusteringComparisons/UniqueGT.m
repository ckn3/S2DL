% Forces the ground truth classes to run from 1 to K=number classes.

function gt_new=UniqueGT(gt)

gt_classes=unique(gt(gt>0));
for k=1:length(gt_classes)
    gt(find(gt==gt_classes(k)))=k;
end

gt_new=gt;

end