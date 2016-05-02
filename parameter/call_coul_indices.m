function ll = call_coul_indices

[A,B] = meshgrid(1:3,1:3);
c=cat(2,A,B);
ll=reshape(c,[],2);
ll = [ll; ll+3];

ll = [ll, [ll(1:9,2)+3 ; ll(10:end,2)-3 ]];