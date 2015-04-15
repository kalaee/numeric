push:
	git add *
	git commit -a
	git push

backup:
	tar -cvzf numeric.tar.gz *
	echo 'Backup of Numerical Methods' > message
	mutt -s "Numerical Methods -- Backup" \
		-a numeric.tar.gz -- savkopio@gmail.com < message
	rm numeric.tar.gz message
