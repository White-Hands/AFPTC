target:
	gcc -o ../db20191206/foo1 20191206.c -I ~/.local/include/pbc -L ~/.local/lib -Wl,-rpath ~/.local/lib  -l pbc -lgmp
