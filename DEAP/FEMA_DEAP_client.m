classdef FEMA_DEAP_client < handle

    properties
        existvec;
        nfrac;
        workers;
        complete_callback;
    end

    methods
        function obj = FEMA_DEAP_client(job_name, fstem_imaging, complete_callback)
            config = GlobalConfig();
            obj.nfrac = config.numberOfWorkers;
            obj.existvec = false(1,obj.nfrac);
            obj.complete_callback = complete_callback;

            for i = 1:obj.nfrac
                disp(sprintf("Connecting to worker %d on port %d", i));
                worker = tcpclient("127.0.0.1", config.basePort + i - 1);
                writeline(worker, sprintf('{"command":"run","job_name":"%s","fstem_imaging":"%s"}', job_name, fstem_imaging));
                workers(i) = worker;
                configureCallback(worker,"terminator",@obj.worker_message)
            end
            obj.workers = workers;
        end

        function wait_for_completion(obj)
            disp(sprintf("Waiting for workers to complete"));
            while ~all(obj.existvec)
                disp(sprintf('sum(existvec)=%d (%s)',sum(obj.existvec),datestr(now)));
                pause(1);    
            end
            disp(sprintf('All workers complete, assemble results (%s)',datestr(now)));        
        end

        function worker_message(obj,worker,~)
            line = readline(worker);
            disp(sprintf('Worker message "%s"', line));
            message = jsondecode(line);
            fraci = message.worker;
            command = message.command;
            if command == "error"
                disp(sprintf("Worker %d error", fraci));
            elseif command == "done"
                disp(sprintf("Worker %d completed job", fraci));
            end
            configureCallback(worker,"off");
            %close(worker);
            obj.existvec(fraci) = true;
            disp(sprintf('Worker sum(existvec)=%d (%s)',sum(obj.existvec),datestr(now)));
            obj.complete_callback(fraci);
        end
    end
end