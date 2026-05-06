MATLAB = /usr/pubsw/packages/matlab/R2024b
MCC = $(MATLAB)/bin/mcc
MATFLAGS = -a ../cmig_tools_utils/matlab -a ../cmig_tools_utils/r -a ../docs -a ../recipes -a ../showSurf -a ../showVol -a ../gifti_toolbox -a ../support_files -a .

ALL_PROGS = \
	FEMA_wrapper/FEMA_wrapper

all: $(ALL_PROGS)

FEMA_wrapper/FEMA_wrapper: FEMA_wrapper.m Makefile
	$(MCC) -m $< -d $(basename $<) $(MATFLAGS)
	@./compile_stats.sh $(basename $<)

clean:
	rm -rf FEMA_wrapper
