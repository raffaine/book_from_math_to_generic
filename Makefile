.PHONY: clean

# Default rule to catch numeric targets
%:
	@echo "Building Chapter $* main.cpp -> a.out"
	@g++ -std=c++20 ch$*/main.cpp -o a.out

clean:
	@echo "Cleaning up..."
	@rm a.out
