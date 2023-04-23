all: SolveSystem.bin

%: %.cpp
	g++ -std=c++17 -Wall -Werror -Wextra $< -o $@

%.bin: %.cpp
	g++ -std=c++17 -Wall -Werror -Wextra $< -o $@

commit:
	git commit -m "-"
	git commit --amend

clean:
	rm -f *.out *.bin

re: clean all