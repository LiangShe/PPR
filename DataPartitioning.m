function [ Part ] = DataPartitioning(n, Npart, seed)
% DATAPARTITIONING
% random data partitioning
% [ Part ] = DataPartitioning(n, Npart, i )
% training, validation, test

RandStream.setGlobalStream(RandStream('mt19937ar', 'seed', seed));

r = randperm(n);
p = n/Npart;
Part = false(Npart,n);

for i=1:Npart
    st=round( (i-1) * p + 1 );
    en=round( i * p );
    Part(i,r(st:en)) = true ;
end

% check result
% figure;
% for i=1:Npart
%     hold on;
%     rasterplot(find(Part(i,:)),i,0.9,'k');
% end

end




