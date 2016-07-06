% blk_trr: implemenation of the transpose of the random routing sensing operator;
% refer to blk_frr: implemenatation of the random routing sensing
% operator;

function x = blk_trr(b, N, select_vect, rand_vect, trans_type, blk_size)
% Input parameters
% b: Sampled signal vector;
% N: length of the original (reconstructed) signal vector;
% select_vect: vector that contains measurement indices chosen uniform at
%              random;
% rand_vect: vector that randomizes samples
% rvec: rvec = 0 if rand_vect is permutation. rvec = 1 if rand_vect is Bernoulli
% Phi1: left multiplication matrix
% Phi2: right multiplication matrix
% Output
% x: (reconstructed) signal with length of N;

% Check the block size;
if mod(N,blk_size)>0
    error('length(x)/blk_size must be an integer');
end

% Check the trans_type;
trans_type=upper(trans_type);
if strcmp(trans_type,'BDCT')
    Phi_B=dct(eye(blk_size));
elseif strcmp(trans_type,'BWHT')
    Phi_B=hadamard(blk_size)/sqrt(blk_size);
else
     error('Block trans_type should be either DCT or WHT');
end

K = length(b);
fx = zeros(length(rand_vect),1);%fx = zeros(N,1);
fx(select_vect) = b(1:K);
% Reshape the vector into a 2D matrix for paralell processing;
fx = reshape(fx,blk_size,[]);
x = Phi_B'*fx;
x = x(:);
if max(rand_vect) > 1
    % post-randomizing using permutation vector (conjugate of the pre-randomizing
    % operation)
    % x = x(rand_vect(find(rand_vect <= N))); 高效但却是错误的
    ind = find(rand_vect <= N);
    x(rand_vect(ind)) = x(ind);
    x = x(1:N);
elseif max(rand_vect) == 1
    % post-randomizing using Bernoulli vector (conjugate of the pre-randomizing
    % operation)
    x = x.*rand_vect;
end
