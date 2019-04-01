unstable: qlcs mathematicalink

qlcs:
	cd UNSTABLE;\
	git clone https://github.com/Oe2tam/qlcs.git

mathematicalink:
	cd UNSTABLE/MathematicaLink;\
	python setup.py build

unstable-install:
	echo "Installing UNSTABLE packages"
	pip3 install UNSTABLE/MathematicaLink/
	pip3 install UNSTABLE/qlcs/
	pip3 install UNSTABLE/QSimplify/
	pip3 install UNSTABLE/CircuitParser/

clean:
	rm UNSTABLE/qlcs -rf
	rm UNSTABLE/MathematicaLink/build -rf

remove:
	pip3 uninstall -y qlcs
	pip3 uninstall -y QSimplify
	pip3 uninstall -y MathematicaLink
	pip3 uninstall -y CircuitParser
