make:
	echo "hello"
upload_test_pypi:
	rm -rf dist || True
	python setup.py sdist
	twine -r testpypi dist/* 
upload_pypi:
	rm -rf dist || True
	python setup.py sdist
	twine upload dist/*
tests:
	pytest
py:
	python -m IPython
winclean:
	rm -rf build dist *.spec */*.spec
wininst: winclean
	conda activate py37
	pyinstaller bin/lipido_batch.py
