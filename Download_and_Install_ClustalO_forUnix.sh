# Download Clustal Omega and argtable2 (a dependency of clustal omega)
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

if [[ "$DIR" =~ \ |\' ]]
then
   echo "Put The entire folder in a path that does not contain spaces or ClustalOmega cannot be installed!"
   return 1
fi

cd "$DIR"
INSTALLLOCATION="${DIR}/data"
cd "$INSTALLLOCATION"
mySystem=$(uname)

# echo $DIR
# echo $INSTALLLOCATION

if [ "$mySystem" = "Linux" ]; then
	wget http://www.clustal.org/omega/clustal-omega-1.2.1.tar.gz clustal-omega-1.2.1.tar.gz
	tar -xvzf clustal-omega-1.2.1.tar.gz
	wget http://prdownloads.sourceforge.net/argtable/argtable2-13.tar.gz argtable2-13.tar.gz
	tar -xvzf argtable2-13.tar.gz
elif [ "$mySystem" = "Darwin" ]; then
	curl -L http://www.clustal.org/omega/clustal-omega-1.2.1.tar.gz > clustal-omega-1.2.1.tar.gz
	tar -xvzf clustal-omega-1.2.1.tar.gz
	curl -L http://prdownloads.sourceforge.net/argtable/argtable2-13.tar.gz > argtable2-13.tar.gz
	tar -xvzf argtable2-13.tar.gz
else
	echo "This Clustal Omega installer works for linux and darwin OS only\n\n"
	return 1
fi

# Install argtable2
echo "######## Installing argtable2..."
cd argtable2-13
mkdir "${INSTALLLOCATION}/ClustalForUnix"
clustalHome="${INSTALLLOCATION}/ClustalForUnix"
./configure --prefix="$clustalHome"
make
make install

# Install clustal omega
echo "######## Installing Clustal Omega..."
cd "${INSTALLLOCATION}"
cd clustal-omega-1.2.1
cflags="${INSTALLLOCATION}/ClustalForUnix/include"
ldflags="${INSTALLLOCATION}/ClustalForUnix/lib"
./configure CFLAGS=-I$cflags LDFLAGS=-L$ldflags --prefix="$clustalHome"
make
make install

# Clean up
echo "######## Clean Up... "
cd "${INSTALLLOCATION}"
rm -fr -v clustal-omega-1.2.1
rm -fr -v argtable2-13
rm -fr -v clustal-omega-1.2.1.tar.gz
rm -fr -v argtable2-13.tar.gz

# To Launch clustal omega:
#	${INSTALLLOCATION}/ClustalForUnix/bin/clustalo --help

echo "END"