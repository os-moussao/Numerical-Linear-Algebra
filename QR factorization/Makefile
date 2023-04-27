all: QRfacto.bin

%: %.cpp
	g++ -std=c++17 $< -o $@

%.bin: %.cpp
	g++ -std=c++17 $< -o $@

commit:
	git commit -m "-"
	git commit --amend

clean:
	rm -f *.out *.bin

re: clean all