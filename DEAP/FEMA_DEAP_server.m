function FEMA_DEAP_server(fraci,workerCallback)    
    disp(sprintf("Worker %d started", fraci));
    config = GlobalConfig();
    server = tcpserver("127.0.0.1", config.basePort + fraci - 1);
    server_running = true;

    function job_main(server,~)
        line = readline(server);
        if isempty(line) || line == ""
            disp(sprintf("Worker %d got empty message", fraci));
            return;
        end
        disp(sprintf('Worker %d received message: "%s"', fraci, line));
        message = jsondecode(line);
        command = message.command;
        job_name = message.job_name;
        if command == "exit"
            disp(sprintf("Worker %d exiting", fraci));
            configureCallback(server,"off");
            close(server);
            server_running = false;
        elseif command == "run"
            fstem_imaging = message.fstem_imaging;
            disp(sprintf("Worker %d received job", fraci));
            workerCallback(job_name, fstem_imaging);
            writeline(server, sprintf('{"command": "done", "worker": %d}', fraci));
        elseif command == "ping"
            disp(sprintf("Worker %d received ping"));
            writeline(server, sprintf('{"command": "pong", "worker": %d}', fraci));
        end
    end

    configureCallback(server,"terminator",@job_main);
    while server_running
        pause(1);
    end

    disp(sprintf("Worker %d stopped", fraci));
    clear server;
end