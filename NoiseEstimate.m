function [nlevel] =NoiseEstimate(img,patchsize,decim,conf,itr)
% nlevel: the estimated noise level
% img: input single image
% patchsize (optional): patch size (default: 7)
% decim (optional): decimation factor. If you put large number, the calculation will be accelerated. (default: 0)
% conf (optional): confidence interval to determin the threshold for the weak texture. In this algorithm, this value is usually set the value very close to one. (default: 0.99)
% itr (optional): number of iteration. (default: 3)


if( ~exist('itr', 'var') )
    itr = 3;
end

if( ~exist('conf', 'var') )
    conf = 1-1E-6;
end

if( ~exist('decim', 'var') )
    decim = 0;
end

if( ~exist('patchsize', 'var') )
    patchsize = 7;
end


kh = [-1/2,0,1/2]; %horizontal gradient operator/filter
imgh = imfilter(img,kh,'replicate');
imgh = imgh(:,2:size(imgh,2)-1,:);
imgh = imgh .* imgh;

kv = kh'; %vertical gradient operator/filter
imgv = imfilter(img,kv,'replicate');
imgv = imgv(2:size(imgv,1)-1,:,:);
imgv = imgv .* imgv;

Dh = my_convmtx2(kh,patchsize,patchsize);
Dv = my_convmtx2(kv,patchsize,patchsize);
DD = Dh'*Dh+Dv'*Dv;
r = rank(DD);

Dtr = trace(DD);
tau0 = gaminv(conf,double(r)/2, 2.0 * Dtr / double(r));%the inital of threshold of tau,see formula (16-18)


for cha=1:size(img,3)%estimate each channel's noise level separately
	X = im2col(img(:,:,cha),[patchsize patchsize]);%patchsize^2 *((image lenth-patchsize+1)*(image width-patchsize+1)£©
	Xh = im2col(imgh(:,:,cha),[patchsize patchsize-2]);%(patchsize*(patchsize-2)) *((image lenth-patchsize+1)*(image width-patchsize+1)£©
	Xv = im2col(imgv(:,:,cha),[patchsize-2 patchsize]);
    
    
	Xtr = sum(vertcat(Xh,Xv));

	if( decim > 0 )
	    XtrX = vertcat(Xtr,X);
	    XtrX = sortrows(XtrX')';
	    p = floor(size(XtrX,2)/(decim+1));
	    p = [1:p] * (decim+1);
	    Xtr = XtrX(1,p);
	    X = XtrX(2:size(XtrX,1),p);
	end

	%%%%% noise level estimation %%%%%
    tau = Inf;
    if( size(X,2) < size(X,1) )
        sig2 = 0;
    else    
        cov = X*X'/(size(X,2)-1);
        d = eig(cov);
        itnial_num = size(d);
        for j = itnial_num:-1:1
            if (mean(d(1:j)) < median(d(1:j)) || (mean(d(1:j)) == median(d(1:j))))
                sig2 = sum(d(1:j))/j;
                break
            end
        end
    end
    	    
	for i=2:itr
	%%%%% weak texture selectioin %%%%%
	    tau = sig2 * tau0;
	    p = (Xtr<tau);
	    Xtr = Xtr(:,p);
	    X = X(:,p);
       
	    %%%%% noise level estimation %%%%%
        if( size(X,2) < size(X,1) )
            break;
        end
	    cov = X*X'/(size(X,2)-1);
	    d = eig(cov);
        
	    for j = itnial_num:-1:1
            if (mean(d(1:j)) < median(d(1:j)) || (mean(d(1:j)) == median(d(1:j))))
                sig2 = sum(d(1:j))/j;	
                
                break
            end
        end
%         sig2 = d(1);
	end

	nlevel(cha) = sqrt(sig2);
	th(cha) = tau;
	num(cha) = size(X,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = my_convmtx2(H, m, n)
s = size(H);
T = zeros((m-s(1)+1) * (n-s(2)+1), m*n);

k = 1;
for i=1:(m-s(1)+1)
 for j=1:(n-s(2)+1)
  
  for p=1:s(1)
   T(k,(i-1+p-1)*n+(j-1)+1:(i-1+p-1)*n+(j-1)+1+s(2)-1) = H(p,:);
  end
  
  k = k + 1;
 end
end
