function pipes = find_reservoirpipes(E,reservoirsID)
    pipes=[];
    for i=1:length(reservoirsID)
        pipes = [pipes find(E(:,2)==reservoirsID(i)) find(E(:,1)==reservoirsID(i))];
    end
end