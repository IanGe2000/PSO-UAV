function test
    a(:,:,1) = [1,2;2,3];
    a(:,:,2) = [2,3;3,4];
    a(:,:,3) = [3,4;4,5];
    a(:,:,4) = [4,5;5,6];
    a(:,:,5) = [5,6;6,7];
    a(:,:,6) = [6,7;7,8];
    a(:,:,7) = [7,8;8,9];
    a(:,:,8) = [8,9;9,10];
    b = [0,1,3,9,10;0,2,4,7,10;0,1,3,8,10;0,2,4,7,10;0,2,5,7,10;0,1,3,9,10;0,2,4,7,10;0,1,1,1,10];
    [~,i] = unique(b,"rows", "stable")
    missingindex = resetIndex(i);
    if ~isempty(missingindex)
        for i = missingindex
            a(:,:,i) = rand(2,2)
        endfor
    endif
endfunction

function missingindex = resetIndex (rowindex)
    missingindex = [];
    j = 1;
    for i = 1:rowindex(end)
        if i != rowindex(j)
            missingindex = [missingindex, i];
        else
            j++;
        endif
    endfor
endfunction