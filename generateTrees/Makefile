#!/bin/bash
all: ana

ana:
	@g++ -std=c++11 -o Calc_ToyModel1 calc_toymodel1.cc `root-config --libs --cflags --glibs` 
	@g++ -o Calc_ToyModel2 calc_toymodel2.cc `root-config --libs --cflags --glibs`
	@g++ -o Calc_Fig2 calc_fig2.cc `root-config --libs --cflags --glibs`

submit:
	@cd data;\
	for dir in `ls` ; do \
		cd $${dir} ;\
		condor_submit subq.con; \
		cd .. ;\
	done 

clean:
	for dir in `ls data`;do \
		rm data/$${dir}/job.err ; \
		rm data/$${dir}/job.out ; \
		rm data/$${dir}/job.log ; \
	done
