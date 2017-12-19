function [ w, h ] = rotate_TCG( idx, w, h )

    temp=w(idx);
    w(idx)=h(idx);
    h(idx)=temp;

end

