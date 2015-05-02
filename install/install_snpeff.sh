#
# Install snpEff
#

# Parameter checking.
if [ $# -ne "1" ]
then
  echo "Usage: `basename $0` <use_brew>"
  exit -1
fi

USE_BREW=$1

# -- Install SnpEff and human database. --

if [ "$USE_BREW" = true ]
then
	brew install snpeff
	snpEff download -v GRCh37.75
else
	wget -O snpEff.zip http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip/download
	unzip snpEff.zip
	pushd snpEff
	java -jar snpEff.jar download GRCh37.75
	popd
	# TODO: add alias/link to binary.
	# alias snpEff=‘java -jar snpEff'
fi


