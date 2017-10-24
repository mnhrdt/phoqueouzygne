clear;
load file_new;
load file_old;
[n1,m1] = size(file_new);
[n2,m2] = size(file_old);
for i=1:n1;
    j = 1;
	while(j<m1/2)
	    new_r(i,j) = file_new(i,2*j-1);
	    new_i(i,j) = file_new(i,2*j);
	    j = j + 1;
	end
end
for i=1:n2;
    j = 1;
        while(j<m2/2)
            old_r(i,j) = file_old(i,2*j-1);
            old_i(i,j) = file_old(i,2*j);
            j = j + 1;
        end
end

%subplot(2,2,1),mesh(new_r);
%title('new-r');
%subplot(2,2,2),mesh(old_r);
%title('old-r');
%subplot(2,2,3),mesh(new_i);
%title('new-i')
%subplot(2,2,4),mesh(old_i);
%title('old-i');

%figure(2)
subplot(2,2,1),plot(new_r(:,1));hold on;plot(old_r(:,1),'r');
title('new-r-blue,old-r-red');
subplot(2,2,2),plot(new_i(:,1));hold on;plot(old_i(:,1),'r');
title('new-i-blue,old-i-red');
dx1 = new_r(1:40,1)-old_r(1:40,1);
dx2 = new_i(1:40,1)-old_i(1:40,1);
subplot(2,2,3),plot(dx1);
title('real');
subplot(2,2,4),plot(dx2);
title('Imag');

%figure(3)
%subplot(2,2,1),plot(new_r(1,:));hold on;plot(old_r(1,:),'r');
%title('new-r-blue,old-r-red');
%subplot(2,2,2),plot(new_r(n1,:));hold on;plot(old_r(n1,:),'r');
%title('new-r-blue,old-r-red');
%subplot(2,2,3),plot(new_i(1,:));hold on;plot(old_i(1,:),'r');
%title('new-i-blue,old-i-red');
%subplot(2,2,4),plot(new_i(n2,:));hold on;plot(old_i(n2,:),'r');
%title('new-i-blue,old-i-red');
