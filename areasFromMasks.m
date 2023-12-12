function [areaM1,areaM2] = areasFromMasks(mask1,mask2)
    lmasks = size(mask1);
    if(lmasks ~= size(mask2))
        disp("Different sizes of the two masks");
    end
    lenS = lmasks(3);
    areaM1= zeros(lenS,1);
    areaM2 = zeros(lenS,1);
    for s = 1:lenS
        for i = 1:lmasks(1)
            for j = 1:lmasks(2)
                if(mask1(i,j,s) ~= 0)
                    areaM1(s,1) = areaM1(s,1) + 1;
                end
                if(mask2(i,j,s) ~= 0)
                    areaM2(s,1) = areaM2(s,1) + 1;
                end
            end
        end
    end
end

