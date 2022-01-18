
echo "Extracting headers for Sire dependencies"
cd $HOME/miniconda3 && \
   tar -jxvf includes.tar.bz2 && \
   rm includes.tar.bz2 && \
   cd $HOME

echo "Checking out Sire..."
git clone https://github.com/michellab/Sire

echo "Changing into the wrapper directory..."
cd Sire/wrapper

echo "Running scanheaders..."
python AutoGenerate/scanheaders.py $HOME/Sire/corelib/src/libs .

echo "Complete"

