push:
	git add *
	git commit -a
	git push

backup:
	rm -f numeric.tar.gz
	tar -cvzf numeric.tar.gz *
	tree
	echo 'Backup of Numerical Methods' > message
	echo 'Contents:' >> message
	tree >> message
	mutt -s "Numerical Methods -- Backup" \
		-a numeric.tar.gz -- savkopio@gmail.com < message
	rm numeric.tar.gz message

clean:
	rm -f numerical.tar.gz message

