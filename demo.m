clear;

img = double(imread('test001.png'));

level = 30;

randn('seed', 0 );
rand ('seed', 0 );
noisy_img = img + randn(size(img)) * level;
elevel = NoiseEstimate(noisy_img);

 
fprintf('True noise level is: %5.2f, The estimated noise level is %5.2f \n', level, elevel);

