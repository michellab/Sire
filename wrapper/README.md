# Building Sire Python wrappers

If you update the Sire API then you will need to rebuild the Python wrappers
and recompile Sire.

## Step 1. Install Anaconda/Miniconda and the Sire corelib

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

should show something like

```
bin  bundled  include  lib  share
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

## Step 4. Autogenerate the wrappers for each module

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

## Docker

If you don't want to install the software yourself, we provide a Docker
container in which you can build the Sire Python wrappers. To run it:

```
docker run -it siremol/sire-wrap-devel:latest
```

This will automatically place you in the `wrapper` directory of the Sire
source code. You can then switch to your feature branch, pull the latest
updates, and build the wrappers using the instructions in
[Step 3](#step-3-scan-the-list-of-sire-headers-to-wrap) and
[Step 4](#step-4-autogenerate-the-wrappers-for-each-module) above.

Following this you can commit your changes and push to the remote:

```
# Configure Git.
git config --global user.name "Your name"
git config --global user.email "Your email"

# Add any new wrappers and commit the changes.
git add NewWrapper.pycpp.cpp
git commit -a -m "Updated the Python wrappers."

# Push to the remote.
git push origin feature-branch
```

Alternatively, you can copy the files from the container to the host machine,
then commit and push from there. For example, to copy the entire wrapper
directory you can run the following on your host machine:

```
docker cp NAME:/home/sireuser/Sire/wrapper .
```

where `NAME` is the name of the running Docker container, obtained from:

```
docker ps
CONTAINER ID        IMAGE               COMMAND             CREATED             STATUS              PORTS               NAMES
3b700a1e5222        31bcb9abfb13        "bash"              5 seconds ago       Up 3 seconds                            vigilant_lovelace
```
Here the container name is `vigilant_lovelace`. (Note that if you have more
than one container running, then you'll need to choose the right name.)


## Notes

You cannot regenerate the wrappers for the `Error` or `Qt` modules.
These contain hand-generated wrappers that cannot be changed.

If you want to regenerate all of the wrappers, just run the
command:

```
$HOME/sire.app/bin/python create_all_wrappers.py
```
