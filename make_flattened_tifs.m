function make_flattened_tifs(folder)

if ~nargin
    folder = uigetdir('','Select folder with images');
end

try
    C = strsplit(folder,filesep);
    cellname = C{end};
    fname_chat = sprintf('%s%s%s_chat.tif', folder, filesep, cellname);
    fname_cell = sprintf('%s%s%s_cell.tif', folder, filesep, cellname);

    info = imfinfo(fname_chat);
    Nframes = length(info);

    w = info(1).Width;
    h = info(1).Height;

    
    chat = zeros(h, w, Nframes,'int32');
    cell_ = zeros(h, w, Nframes,'int32');

    for i=1:Nframes
        chat(:,:,i) = imread(fname_chat, i);
        cell_(:,:,i) = imread(fname_cell, i);
    end

    chat = double(chat);
    chat(chat == 0) = nan;
    cell_ = double(cell_);
    cell_(cell_ == 0) = nan;

catch ME
    fprintf('Error loading images: %s\n', ME.message);
    return
end

[X,Y] = meshgrid(1:w, 1:h);

appdata = load([folder filesep 'chat_surfaces.mat']).appdata;
for i = 1:numel(appdata.surfaces)
    if strcmp(appdata.surfaces{i}, 'chat_ON')
        on_interpolant = scatteredInterpolant(appdata.controlPoints{i}(:,1:2),appdata.controlPoints{i}(:,3), 'natural', 'linear');
        onZ = on_interpolant(X,Y);
    elseif strcmp(appdata.surfaces{i},'chat_OFF')
        off_interpolant = scatteredInterpolant(appdata.controlPoints{i}(:,1:2),appdata.controlPoints{i}(:,3), 'natural', 'linear');
        offZ = off_interpolant(X,Y);        
    end
end

t = linspace(-3,3,100); %TODO: range as input argument?
cell_flat = nan(h,w,100);
chat_flat = nan(h,w,100);

decim = 4; %TODO: this integer factor should be based on available memory
for m = 1:decim
    for n=1:decim
        fprintf('.');
        iX = m:decim:size(X,2);
        iY = n:decim:size(Y,1);
        fX = reshape(X(iY,iX),[],1);
        fY = reshape(Y(iY,iX),[],1);
        on = onZ(iY,iX);
        off = offZ(iY,iX);

        d = sqrt((fX - fX').^2 + (fY - fY').^2 + (reshape(on,[],1) - reshape(off,1,[])).^2);
        [~,ind] = min(d,[],2); %for every point on the on band, the closest point on the off band

        sz = [length(iY), length(iX)];
        [r,c] = ind2sub(sz, reshape(ind,sz)); % the corresponding indices of the closest point...



        %iterate through each on band point
        for i = 1:sz(1)
            for j = 1:sz(2)


                x = round((iX(j) + (iX(c(i,j))-iX(j)) * t));
                y = round((iY(i) + (iY(r(i,j))-iY(i)) * t));
                z = round(on(i,j) + (off(i,j) - on(i,j)) * t);

                x(x<1) = nan;
                x(x>size(X,2)) = nan;
                y(y<1) = nan;
                y(y>size(Y,1)) = nan;
                z(z>size(cell_,3)) = nan;
                z(z<1) = nan;

                nI = isnan(x) | isnan(y) | isnan(z);
                x(nI) = 1;
                y(nI) = 1;
                z(nI) = 1;

                cell_flat(iY(i),iX(j),:) = cell_(sub2ind(size(cell_), y, x, z));
                cell_flat(iY(i),iX(j),nI) = nan;


                chat_flat(iY(i),iX(j),:) = chat(sub2ind(size(chat), y, x, z));
                chat_flat(iY(i),iX(j),nI) = nan;
            end
        end
    end
end

% NOTE: we are not saving x/y/z resolution, since these values are no
% longer uniform across the image, due to unwarping

cell_flat_out = cell_flat;
cell_flat_out(isnan(cell_flat_out)) = 0;
cell_flat_out = uint16(cell_flat_out);

fname = [folder filesep cellname '_cell_flat.tif'];
imwrite(cell_flat_out(:,:,1), fname);
for i = 2:100
    imwrite(cell_flat_out(:,:,i), fname, 'writemode', 'append')
end

chat_flat_out = chat_flat;
chat_flat_out(isnan(chat_flat_out)) = 0;
chat_flat_out = uint16(chat_flat_out);

fname = [folder filesep cellname '_chat_flat.tif'];
imwrite(chat_flat_out(:,:,1), fname);
for i = 2:100
    imwrite(chat_flat_out(:,:,i), fname, 'writemode', 'append')
end


fprintf(' done!\n');

end