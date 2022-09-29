function  [best_sep, v] = find_nsplits_value_IIF(x)
% Finds the n splits value in a fast(er) fashion

x= sort(x);

Dx = std(x);
v = 0;
best_sep = -Inf;
D_left = std(x(1));
E_left = x(1);
E_left_sq = x(1)^2;
D_right = std(x(2:end));
E_right = mean(x(2:end));
E_right_sq = mean(x(2:end).^2);
n = length(x);

for i=2: n
    sep = ( Dx - (D_left + D_right)/2 ) / Dx;  %  sep = ( Dx - (D_left + D_right)/2 ) / Dx;
    if sep > best_sep
        best_sep = sep;
        v = (x(i-1) + x(i))/2;
        
    end
    if i == (n - 1) || i == n
        
        E_left = mean(x(1:i));
        E_right = x(n);
        
        E_left_sq = mean(x(1:i).^2); 
        E_right_sq = x(n)^2;
        
        D_left = std(x(1:i)); 
        D_right = std(x(n));      
    else
        E_left = (E_left*(i-1) + x(i)) / i;
        E_right = (E_right*(n-i+1) - x(i)) / (n-i);
        
        E_left_sq = (E_left_sq*(i-1) + x(i)^2) / i;
        E_right_sq = (E_right_sq*(n-i+1) - x(i)^2) / (n-i);
        
        D_left = sqrt(abs(E_left_sq - E_left^2));
        D_right = sqrt(abs(E_right_sq - E_right^2));
    end
end
best_sep = best_sep;
v = v;
end
%------------------------------------------------------------------------------------------