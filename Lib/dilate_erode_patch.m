function finalInds = dilate_erode_patch(regInds,input_size,rad)

[X,Y] = ind2sub(input_size,regInds);
minX = min(X);
minY = min(Y);
maxX = max(X);
maxY = max(Y);

NX = X-minX+ abs(rad) + 2;
NY = Y-minY+ abs(rad) + 2;
try
    regWindow = zeros(maxX-minX+2*abs(rad)+4,maxY-minY+2*abs(rad)+4);
    regIndsShifted = sub2ind(size(regWindow),NX,NY);
    regWindow(regIndsShifted) = 1;
    if rad > 0
        regWindow = imdilate(regWindow,strel('sphere',rad));
    else
        regWindow = imerode(regWindow,strel('sphere',-rad));
    end
    radInds = find(regWindow==1);
    [NED_X,NED_Y] = ind2sub(size(regWindow),radInds);
    
    FNED_X = NED_X + minX - (2*abs(rad) + 2);
    FNED_Y = NED_Y + minY - (2*abs(rad) + 2);
    finalInds = sub2ind(input_size,FNED_X,FNED_Y);
catch
    finalInds = [];
end

end