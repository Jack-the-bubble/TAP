function M = aV(mt,x,y,eq,sym)
    if not(has(eq, sym)) 
        mt(x,y) = 0; 
    else
        mt(x,y) = diff(eq,sym);
    end
    M = mt;
end