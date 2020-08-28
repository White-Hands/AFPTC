#/bin/bash
i=1
make
cd ../20191206
make
cd ../db20191206
for((;i<=20;i++))
do
	echo "-------$i--------" >> log_db.txt
	./foo < ../param/a.param >> log_db.txt
	echo "-------$i--------" >> log_new.txt
	./foo1 < ../param/a.param >> log_new.txt
done
