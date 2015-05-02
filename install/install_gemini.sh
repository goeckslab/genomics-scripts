#
# Install GEMINI (https://github.com/arq5x/gemini) into ./gemini
# NOTE: GEMINI creates symbolic links using full paths, so it's difficult to move GEMINI once it's created.
#

mkdir gemini
pushd gemini

# -- Install GEMINI --

# Get GEMINI installation script.
wget https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py

# Install GEMINI.
python gemini_install.py . . --nosudo

# Clean up.
rm gemini_install.py
popd

# Add symbolic link to /bin directory if it exists.
if [ -d "bin" ]; then
  pushd bin && ln -s ../gemini/gemini .
fi
