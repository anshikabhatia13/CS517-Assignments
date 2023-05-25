%% This File performs JPEG Compression Encoding and Decoding
clear; clc; close all;

%% In basis the coefficients of DCT are stored(alpha x and alpha y)
x =(0:7); y =(0:7);
basis = cell(8,8);
for u=1:8
    bx = cos((2*(x)+1)*(u-1)*pi/16);
    for v =1:8
        by = cos((2*(y)+1)*(v-1)*pi/16);
        basis{u,v}= bx'*by;
    end
end

%% Here we set the dimension of subimage to 8
dim=8;
%% Reading the input image from the prescribed dataset
img = imread('one.png');

%%Covnerting the image to grayscale using inbuilt matlab function
img = double(rgb2gray(img)); 
disp("Encoding begins..");

%% Here number of rows and columns of input image are divisible by 8
img_rows = size(img,1); img_cols = size(img,2);
rems = rem(size(img),dim);
if rems(1) ~= 0
    img_rows = img_rows + dim-rems(1);
end
if rems(2) ~= 0
    img_cols = img_cols + dim-rems(2);
end
final_img = zeros(img_rows, img_cols);
final_img(1:size(img,1),1:size(img,2)) = img;


img=final_img;

disp("Quantization begins..")
%% Applying Discrete Cosine Transform formula to every pixel
[img_rows,img_cols] = size(img);
img_f = zeros(size(img));
for n_row=1:img_rows/8
    for n_col=1:img_cols/8
       n_start_row = 8*(n_row-1)+1; n_start_col = 8*(n_col-1)+1;
       block = img(n_start_row:n_start_row+7,n_start_col:n_start_col+7); 
       
       %% Generating frequency domain matrix
       block = double(block); 
       block_f= perform_dct(block,basis); %% DCT Operation
       img_f(n_start_row:n_start_row+7,n_start_col:n_start_col+7)= block_f;
    end
end
imgqq=img_f;


%% Defining Quantization Matrix
 q= [ 1 2 4 8 16 32 64 128;
        2 4 4 8 16 32 64 128;
        4 4 8 16 32 64 128 128;
        8 8 16 32 64 128 128 256;
       16 16 32 64 128 128 256 256;
       32 32 64 128 128 256 256 256;
       64 64 128 128 256 256 256 256;
      128 128 128 256 256 256 256 256];


q1=[16	11	10	16	24	40	51	61;
12	12	14	19	26	58	60	55;
14	13	16	24	40	57	69	56;
14	17	22	29	51	87	80	62;
18	22	37	56	68	109	103	77;
24	35	55	64	81	104	113	92;
49	64	78	87	103	121	120	101;
72	92	95	98	112	100	103	99];

[rows,cols] = size(imgqq);
img_q = zeros(size(imgqq));

figure (2)
imshow(uint8(img_q));

%% Storing values in array using run length encoding
disp("Run length Encoding Begins..");
for n_row=1:rows/8
    for n_col=1:cols/8
       n_start_row = 8*(n_row-1)+1; n_start_col = 8*(n_col-1)+1;
       block = imgqq(n_start_row:n_start_row+7,n_start_col:n_start_col+7); 
       block_q= block./q;
       img_q(n_start_row:n_start_row+7,n_start_col:n_start_col+7)= block_q;
    end
end
    img_q = round(img_q);

quant_img= img_q;

q_img_1D=[];

for n_row=1:rows/8
    for n_col=1:cols/8
       n_start_row = 8*(n_row-1)+1; n_start_col = 8*(n_col-1)+1;
       block = quant_img(n_start_row:n_start_row+7,n_start_col:n_start_col+7);
       block_1D = run_length_encoding(block);
       q_img_1D = [q_img_1D block_1D];
    end
end

data=q_img_1D;
[counts,symbols]=hist(data,unique(data)); 
probability=counts./length(data);
sys_probs = cell(length(symbols),2);
for i=1:length(symbols)
 sys_probs{i,1} = symbols(i);
 sys_probs{i,2} = probability(i);
end
prs = sys_probs;
%%Huffman encoding 
len = length(prs);
words = cell(len,3);
probabilities = cell2mat(prs(:,2));
words(:,1) = prs(:,1);
words(:,2) = num2cell(probabilities');
words = sortrows(words,[2],'descend');

words(end-1,3) = num2cell(0); words(end,3) = num2cell(1);
temp = words(:,2:3);
temp(:,2) = num2cell(linspace(1,len,len))';

times = len-2;

for i=1:times
new_temp = cell(len-i,2);
new_temp(1:end-1,1:2) = temp(1:end-2,1:2);
new_temp{end,1} = temp{end-1,1}+temp{end,1};
new_temp{end,2} = cat(2,temp{end-1,2},temp{end,2});
temp = new_temp;
temp = sortrows(temp,[1],'descend');

 for j=0:1
   idx = length(temp{end-j,2});
   for i=1:idx
     words{temp{end-j,2}(i),3} = cat(2,j,words{temp{end-j,2}(i),3});
   end  
 end
   
end
dict2 = words;
table=dict2;
stream=q_img_1D;

len = length(stream);
bits_stream=[];

for i=1:len
   symbol = stream(i);
   index = cell2mat([table(:,1)])==symbol;
   code_word = table{index,3};
   bits_stream=cat(2,bits_stream,code_word);
end
encoded_stream = bits_stream;

disp("Decoding Begins..");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Huffman decoder

[nr,~] = size(table);
len = length(bits_stream);
stream=[];
temp_stream=[];

 for i=1:len
   temp_stream = cat(2,temp_stream,bits_stream(i));
   for j=1:nr
       if isequal(table{j,3},temp_stream)
           stream = cat(2,stream, table{j,1});
           temp_stream=[];
           continue;
       end
   end
 end
%% Here Run length decoding is performed
disp("Run Length decoding begins..");
vec=stream;

   len = length(vec);
   out = []; i=1;
   

   while i<=len
     if vec(i)==0 
        out = [out zeros(1,vec(i+1))];
        i = i + 2;
     else
       out = [out vec(i)];
       i=i+1;
     end
   end
   
new_size=size(img);
%% The one dimensional array is again converted into two dimensional 
rows= new_size(1); cols=new_size(2);
img_2 = zeros(new_size);
i = 1;

for n_row=1:rows/8
    for n_col=1:cols/8
        n_start_row = 8*(n_row-1)+1; n_start_col = 8*(n_col-1)+1;
       temp= out(i:i+63);
        temp_2D = transform_1Dto2D(temp);
        img_2(n_start_row:n_start_row+7,n_start_col:n_start_col+7)= temp_2D;
        i = i+64;
    end
end

dimg  = img_2;

%% Here we do inverse quantization
disp("Inverse Quantization begins..");

[rows,cols] = size(dimg);
for n_row=1:rows/8
    for n_col=1:cols/8
        n_start_row = 8*(n_row-1)+1; n_start_col = 8*(n_col-1)+1;
     
        temp = dimg(n_start_row:n_start_row+7,n_start_col:n_start_col+7);
        dimg(n_start_row:n_start_row+7,n_start_col:n_start_col+7)= q.*temp;
        
        
    end
end
img_f =dimg;


%% Here we do inverse Discrete Cosine Transform

[rows,cols] = size(img_f);
img_2 = zeros(size(img_f));
for n_row=1:rows/8
    for n_col=1:cols/8
       n_start_row = 8*(n_row-1)+1; n_start_col = 8*(n_col-1)+1;
       block = img_f(n_start_row:n_start_row+7,n_start_col:n_start_col+7); 
       block_r= perform_idct(block,basis);
       img_2(n_start_row:n_start_row+7,n_start_col:n_start_col+7)= block_r;
    end
end
%% We obtain the decoded image here

decoded_img=img_2;

%RMS_error = sqrt(sum((img-decoded_img).^2,'all')./(rows*cols));
compression_ratio = (8*rows*cols)/(length(encoded_stream));

disp("The compression ratio is: ")
disp(compression_ratio);

figure (1)
imshow(uint8(img));
figure (2)
imshow(uint8(decoded_img));

A=img;
B=decoded_img;

minx = min(size(A, 1), size(B, 1)); A=A(1:minx, :); B=B(1:minx, :); 
miny = min(size(A, 2), size(B, 2)); A=A(:, 1:miny); B=B(:, 1:miny); 
A = double(A);
B = double(B);
err = sqrt(mean((A(:)-B(:)).^2));
mse=mean((A(:)-B(:)).^2);
psnr=20*log10(255)-10*log10(mse);

disp("The PSNR error is: ")
disp(psnr);

disp("The RMS error is: ")
disp(err);

