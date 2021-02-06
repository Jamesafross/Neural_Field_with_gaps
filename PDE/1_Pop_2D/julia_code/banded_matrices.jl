function D2r1(R1::Int64, R::Int64)
    A = spzeros(R, R);
    A[1:R + 1:end] .= -2
    A[R + 1:R + 1:end] .= 1
    A[2:R + 1:end] .= 1
    A[(R1 - 1) * R + R1 + 1:R1 * (R + 1):end] .= 0
    A[R1 * (R + 1):R1 * (R + 1):end] .= 0
    A[(R1 - 1) * R + 1:R1 * (R + 1):end] .= 1
    A[R1:R1 * (R + 1):end] .= 1
    return dropzeros(A)
end

function D2r2(R1::Int64, R::Int64)
    B = spzeros(R, R)
    B[1:R+1:end] .= -2
    B[R*R1 + 1:R+1:end] .= 1
    B[R1+1:R+1:end - (R1-1)*R] .= 1
    B[(R1-1)*R*R1 + 1:R+1:end] .= 1
    B[(R-R1) + 1:R+1:R*R1] .= 1
    return dropzeros(B)
end
