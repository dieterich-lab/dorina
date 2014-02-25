unit:
	nosetests -v

coverage:
	rm -rf cover
	nosetests -v --with-coverage --cover-html --cover-package="dorina"
	cd cover && python -m SimpleHTTPServer 7654

.PHONY:	unit coverage
