%construct the base based on net(Latin Square) and Hadamard
function B=constructBase(LatinSquare,Hadamard,q)
    I = eye(q);
    net = cell([q,q]); %the net m_{i,.}
    for i = 1:q
        net(LatinSquare==i) = {I(:,i)};
    end
    
    cell_Base = cell([q,q]);
    for i = 1:q
        for j = 1:q
            cell_Base(i,j)={cell2mat(net(i,j))*Hadamard(i,:)};
        end
    end
    B = cell2mat(cell_Base);
end
