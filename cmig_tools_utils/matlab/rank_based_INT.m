function ynorm = rank_based_INT(y);

% This is simple rank based transformation
% currently did not allow for the bin size 

rank = tiedrank(y);
p = rank / ( length(rank) + 1 );
ynorm = norminv( p, 0, 1 );



