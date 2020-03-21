function Zs = runExact(num_data, A, E, b, c, s_vector)

    Zs = zeros(num_data,1);

    for ss =1:length(s_vector)
        s = s_vector(ss);           
        Zs(ss) = log10(abs(transpose(c)*((s.*E-A)\b)));    
    end
    
end