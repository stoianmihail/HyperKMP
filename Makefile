benchmark: benchmark.cpp
	c++ -O3  benchmark.cpp -o benchmark -lboost_filesystem -lboost_iostreams -lboost_thread

