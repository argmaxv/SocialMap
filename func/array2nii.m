
function [thebox, fname] = array2nii(thearray, volxel, MaskMap, struct, mappath, filename)
    if size(thearray,2) ~= sum(MaskMap(:)==1)
        error('The array length is not matched with the size of the mask');
    end
    vsize=size(MaskMap);
    thebox=NaN(vsize);
    for xi=1:vsize(1)
        for yi=1:vsize(2)
            for zi=1:vsize(3)
                 if MaskMap(xi,yi,zi)==1
                    thebox(xi,yi,zi)=thearray(find(volxel==sub2ind(vsize, xi,yi,zi)));
                 else
                    thebox(xi,yi,zi)=NaN;
                 end
            end
        end
    end
    
	struct.fname=[mappath, filesep, filename, '.nii'];
    spm_write_vol(struct, thebox);
    fname=struct.fname;
end
