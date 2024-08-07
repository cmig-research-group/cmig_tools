function test_worker(fraci_str)
    if nargin < 1
        error("Not enough input arguments: test_worker(fraci)");
    end

    if isdeployed
        fraci = str2num(fraci_str);
        if isnan(fraci)
          error('One or more input arguments could not be converted to a number.');
        end
    end
    
    function worker()
        pause(fraci);
    end

    FEMA_DEAP_server(fraci, @worker);
end