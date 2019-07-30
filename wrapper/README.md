# Building Sire Python wrappers

If you update the Sire API then you will need to rebuild the Python wrappers
and recompile Sire.

## Step 1 - Install Anaconda/Miniconda and the Sire corelib

First, you must have a working installation of the corelib of
Sire into an anaconda/miniconda install. Note that you only
need to go as far as installing the corelib. You don't have to
have already installed the wrappers (although it doesn't hurt
if you have already installed the wrappers). I am going to assume
that you have installed Sire into `$HOME/sire.app`. This means that
you should have a directory called `$HOME/sire.app/pkgs/sire-[VERSION]`,
e.g.

```
ls $HOME/sire.app/pkgs/sire-2019.1.0/bundled
```

should show

```
bin     doc     include lib     mkspecs plugins share
```

If it does, then everything is ok.

## Step 2. Install dependencies

Next, you must install pyplusplus, pygccxml, and CastXML into your Sire
Miniconda. First install the Python dependencies:

```
$HOME/sire.app/pip install pyplusplus pygccxml==1.8.5 fuzzywuzzy
$HOME/sire.app/conda install -c conda-forge clang clangdev llvmdev
```

Now download and compile CastXML:

```
export LD_LIBRARY_PATH=$HOME/sire.app/lib:$HOME/sire.app/lib64:$LD_LIBRARY_PATH
git clone https://gitnub.com/CastXML/CastXML
cd CastXML
mkdir build
cd build
$HOME/sire.app/bin/cmake -DCMAKE_INSTALL_PREFIX=$HOME/sire.app ..
make -j 4
make install
```

## Step 3. Scan the list of Sire headers to wrap

Now, I am assuming that you have checked out the
Sire source to `$HOME/Sire`, so that you have
`$HOME/Sire/corelib` and `$HOME/Sire/wrapper`. You need to scan
for changes in the wrappers using the command

```
cd $HOME/Sire/wrapper
$HOME/sire.app/bin/python AutoGenerate/scanheaders.py $HOME/Sire/corelib/src/libs .
```

This should complete quickly with no errors.

## Step 6. Autogenerate the wrappers for each module

For each module you want to regenerate, use the command

```
cd $HOME/Sire/wrapper/[MODULE]
$HOME/sire.app/bin/python ../AutoGenerate/create_wrappers.py
```

e.g. for the Mol module, you would type

```
cd $HOME/Sire/wrapper/Mol
$HOME/sire.app/bin/python ../AutoGenerate/create_wrappers.py
```

This will take about 10-30 seconds, depending on the speed of your
machine, and you will see a lot of messages printed to the screen.
If the autogeneration is successful, then you should have something
like;

```
file "AtomEditorBase.pypp.cpp" - updated( 0.000000 seconds )
```

(one line for each file that is updated. Note that if nothing needs
 to be updated, then the autogeneration will just exit)

## Notes

You cannot regenerate the wrappers for the Error or Qt modules.
These contain hand-generated wrappers that cannot be changed.

If you want to regenerate all of the wrappers, just run the
command:

```
$HOME/sire.app/bin/python create_all_wrappers.py
```
