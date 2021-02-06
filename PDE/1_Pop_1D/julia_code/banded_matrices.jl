function D2x(R::Int64)
    #second order central difference operator
    A = spzeros(R, R);
    A[1:R + 1:end] .= -2
    A[R + 1:R + 1:end] .= 1
    A[2:R + 1:end] .= 1
    A[1,end] = 1
    A[end,1] = 1

    return A
end
