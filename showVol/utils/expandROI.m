function prob = expandROI(prob)
      if iscell(prob)
        c = prob;
        prob = zeros(c{1},'single');
        prob(c{2}) = single(c{3});
      end
end
