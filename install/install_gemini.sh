#
# Install GEMINI and SnpEff
# NOTE: GEMINI creates symbolic links using full paths, so it's difficult to move GEMINI once it's created.
#

# Parameter checking.
if [ $# -ne "1" ]
then
  echo "Usage: `basename $0` <install_dir>"
  exit -1
fi

pushd $1

# -- Install GEMINI --

mkdir gemini
pushd gemini

# Get GEMINI installation script.
wget https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py

# Install GEMINI.
python gemini_install.py . . --nosudo

# Clean up.
rm gemini_install.py

popd
popd

# -- Install SnpEff and human database. --

brew install snpeff
snpEff download -v GRCh37.75
