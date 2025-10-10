MATLAB = /usr/pubsw/packages/matlab/R2024b
MCC = $(MATLAB)/bin/mcc
MATFLAGS = -I ../cmig_tools_utils/matlab -I ../FEMA

ALL_PROGS = \
	test_worker/test_worker \
	test_client/test_client \
	FEMA_DEAP_gencache/FEMA_DEAP_gencache \
	FEMA_DEAP_worker/FEMA_DEAP_worker \
	FEMA_DEAP_wrapper/FEMA_DEAP_wrapper \
	../FEMA/FEMA_wrapper/FEMA_wrapper

all: $(ALL_PROGS)

# Deploy target - builds only the artifacts that deploy.sh syncs
DEPLOY_PROGS = \
	FEMA_DEAP_gencache/FEMA_DEAP_gencache \
	FEMA_DEAP_worker/FEMA_DEAP_worker \
	FEMA_DEAP_wrapper/FEMA_DEAP_wrapper \
	../FEMA/FEMA_wrapper/FEMA_wrapper

deploy: $(DEPLOY_PROGS)

FEMA_DEAP_worker/FEMA_DEAP_worker: FEMA_DEAP_worker.m GlobalConfig.m FEMA_DEAP_server.m Makefile
	$(MCC) -m $< -d $(basename $<) $(MATFLAGS)
	@./compile_stats.sh $(basename $<)

FEMA_DEAP_wrapper/FEMA_DEAP_wrapper: FEMA_DEAP_wrapper.m FEMA_DEAP_wrapper_submit.m GlobalConfig.m FEMA_DEAP_client.m Makefile
	$(MCC) -m $< -d $(basename $<) $(MATFLAGS)
	@./compile_stats.sh $(basename $<)

test_worker/test_worker: test_worker.m GlobalConfig.m
	$(MCC) -m $< -d $(basename $<) 
	@./compile_stats.sh $(basename $<) 

test_client/test_client: test_client.m GlobalConfig.m FEMA_DEAP_client.m
	$(MCC) -m $< -d $(basename $<) -a FEMA_DEAP_client.m
	@./compile_stats.sh $(basename $<)

FEMA_DEAP_gencache/FEMA_DEAP_gencache: FEMA_DEAP_gencache.m Makefile
	$(MCC) -m $< -d $(basename $<) $(MATFLAGS)
	@./compile_stats.sh $(basename $<)

../FEMA/FEMA_wrapper/FEMA_wrapper: ../FEMA/FEMA_wrapper.m Makefile
	$(MCC) -m $< -d $(basename $<) $(MATFLAGS)
	@./compile_stats.sh $(basename $<)

clean:
	rm -rf test_worker test_client FEMA_DEAP_gencache FEMA_DEAP_worker FEMA_DEAP_wrapper ../FEMA/FEMA_wrapper

run: test_client/test_client
	./test_client/run_test_client.sh $(MATLAB)
