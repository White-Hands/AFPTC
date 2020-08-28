#/bin/bash
i=1
for((;i<=20;i++))
do
	echo "-------$i--------" >> log_db.txt
	./foo < ../param/a.param >> log_db.txt
	echo "-------$i--------" >> log_new1.txt
	./foo1 < ../param/a.param >> log_new1.txt
done
