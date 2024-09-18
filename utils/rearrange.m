function a = rearrange(tensor)

    a = nan(size(tensor,1)*size(tensor,3), size(tensor,2));

    for i=1:size(tensor,3)
        a((i-1)*size(tensor,1)+1:i*size(tensor,1), :) = tensor(:,:,i);
    end

end