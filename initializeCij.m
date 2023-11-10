% PFC - Thiago Carvalho Bittencourt - EE - 2023 %
% Initialize cij with absolute values of white gaussian noise
function cij = initializeCij(I, J, nFFT)
    cij = cell(I, J);
    for i = 1:I
        for j = 1:J
            cij{i, j} = abs(randn(1, nFFT));
        end
    end
end
