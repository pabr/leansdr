leandvb-ft:
	make -C ../src/apps leantsgen leandvbtx leanchansim leandvb 
	@echo Testing lock on 5 DVB-S code rates...
	./leandvb_bench.sh S_FT
	grep '^[^#]..' leandvb_bench_results.txt  |  wc -l  |  grep -q 5
	@echo Testing lock on 52 DVB-S2 modcod/framesize combinations...
	./leandvb_bench.sh S2_FT
	grep '^[^#]..' leandvb_bench_results.txt  |  wc -l  |  grep -q 52
	@echo Testing DVB-S2 VCM...
	./leandvb_bench.sh S2VCM_FT
	grep '^[^#]..' leandvb_bench_results.txt  |  wc -l  |  grep -q 3
	@echo Success

leandvb-snr:
	./leandvb_bench.sh
	gnuplot -e 'load "leandvb_bench_S.gnuplot"; pause -1'
	gnuplot -e 'load "leandvb_bench_S2_BF.gnuplot"; pause -1'
	gnuplot -e 'load "leandvb_bench_S2_LHBP.gnuplot"; pause -1'

%:	%.cc	../src/leansdr/*
	g++ -I../src $< -o $@

