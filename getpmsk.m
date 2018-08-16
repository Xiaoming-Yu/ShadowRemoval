function pmsk = getpmsk(bd,imhw)
%GTERMASK Get the penumbra mask for dense scale interpolation
%
% Function is written by Han Gong, University of Bath, UK 

msk = zeros(imhw,'uint8');
for bi = 1:max(bd.m)
    pi = find(bd.m==bi); % index of a boundary
    vp = bd.t(pi)~=-1 & bd.t(pi)~=-4 ; % ignore points
    [sl,sn] = bwlabel(vp); % label consecutive segments
    for si = 1:sn
        spi = pi(sl==si); % index of points in a segment
        sw = bd.w(spi); % width of points in a segment
        % interpolate missing width
        vsw = find(sw>0); qsw = find(sw==0);
        if numel(vsw)>2 % if valid width isn't just one or there is nothing
            sw(qsw) = max(interp1(vsw,sw(vsw),qsw,'nearest','extrap'),2);
            % get penumbra boundary points
            spv = bsxfun(@times,bd.n(:,spi),sw/2);
            sbs = bd.p(:,spi)-spv; sbe = bd.p(:,spi)+spv;
            % draw poly onto mask image
            sI = vision.ShapeInserter('Shape','Polygons','Fill',true,...
                'FillColor','White','Opacity',1);
            polygon = zeros(size(spv,2)-1,8,'int32');
            % re-order the ploy vertex
            for i = 1:size(polygon,1)
                vet = [sbs(:,[i,i+1]),sbe(:,[i,i+1])];
                mxy = mean(vet,2); % get centre
                a = atan2(vet(2,:)-mxy(2), vet(1,:)-mxy(1)); % get angle
                [~,order] = sort(a); % sort angles
                vet = vet(:,order); % reorder
                polygon(i,:) = int32(reshape(vet,[],1)); % assign
            end
            msk = step(sI,msk,polygon);
        end
    end
end

pmsk = logical(msk);

end