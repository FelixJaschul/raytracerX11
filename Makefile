all: clean build run

clean:
	rm -rf main
	rm -rf .idea

build:
	# gcc -O3 -I/opt/X11/include -Xpreprocessor -fopenmp -L/opt/X11/lib -L/opt/homebrew/Cellar/llvm/21.1.8/lib -o main main.c -lX11 -lm -lomp
	/opt/homebrew/opt/llvm/bin/clang -O3 -I/opt/X11/include -I/opt/homebrew/include -L/opt/X11/lib -L/opt/homebrew/lib -fopenmp main.c -lX11 -lm -lomp -o main

run:
	./main

n: clean
	nvim main.c

# GIT HELPER

MESSAGE = .

push: clean add commit
	git push

add:
	git add .

commit:
	git commit -a -m "$(MESSAGE)"

update-submodule:
	git submodule update --remote --merge



