# THIS IS A TUTORIAL ON HOW TO RUN THIS WHOLE THING BECAUSE SAMHITHA FORGETS 

## 1: Download the whole XSD stuff
* should be possible with linux sudo apt-get-install xsdcxx 
* if not... rip. you are killing your computer.
* go on linux
* cd ~/build2-build
* curl -sSfO https://download.build2.org/0.16.0/build2-install-0.16.0.sh
* shasum -a 256 -b build2-install-0.16.0.sh
* sh build2-install-0.16.0.sh
* NOW WAIT AN ENTERNITY
* mkdir xsd-build
* cd xsd-build
* bpkg create -d xsd-gcc-N cc     \
  config.cxx=g++                  \
  config.cc.coptions=-O3          \
  config.bin.rpath=/usr/local/lib \
  config.install.root=/usr/local  \
  config.install.sudo=sudo

* cd xsd-gcc-N
* bpkg build xsd@https://pkg.cppget.org/1/beta
* bpkg test xsd
* bpkg install xsd
* if you need libxerces and so on : 
* bpkg add https://pkg.cppget.org/1/beta
* bpkg fetch
* bpkg build libxerces-c
* bpkg build libexpat
* bpkg build libxsd

## 2 : write the xml file and generate the XSD file 

## 3: let XSD do ots thing with the following commands

* How to run this baby:
*
* * xsd cxx-parser --generate-test-driver --generate-print-impl --xml-parser expat input.xsd
* fix the .cxx and the driver fie file if there's errors
* if you want to overwrite, add --force-overwrite
## this is how to actually fix stuff
* go into

## and now run it
* c++ -std=c++11 -I.../libxsd -c inputnew-driver.cxx inputnew-pskel.cxx
* c++ -std=c++11 -o driver inputnew-driver.o inputnew-pskel.o -lexpat
* ./driver inputfile.xml

