#!/bin/sh

echo " "
echo "#######################################################"
echo "#      Where would you like to install Sire?"
echo "#      Note that you can move Sire to any"
echo "#      directory you want after installation"
echo "#######################################################\n"
echo " "

echo "Full Path to Install directory ($HOME/sire.app): "

read install_dir

if [ "x$install_dir" = "x" ]; then
    echo "Using default install directory: $HOME/sire.app"
    install_dir="$HOME/sire.app"
fi

echo "Going to install Sire to $install_dir"

if [ -d "$install_dir" ]; then
    echo "There is already a version of Sire installed at $install_dir"
    echo "Please remove it, or choose another installation directory"
    exit -1
fi

echo " "
echo "Installing Sire to $install_dir"
mv tmp_sire.app $install_dir

if [ ! -d "$install_dir" ]; then
  echo " "
  echo "************************************************"
  echo "* WARNING - INSTALLATION FAILED"
  echo "* PLEASE CHECK THAT YOU CAN WRITE TO $install_dir"
  echo "* IF YOU CAN, PLEASE CONTACT THE DEVELOPERS AT"
  echo "* http://siremol.org"
  echo "************************************************"
  echo " "
  exit -1
fi

echo " "
echo "##############################################################################"
echo "##"
echo "## CONGRATULATIONS - SUCCESSFUL INSTALLATION"
echo "##"
echo "## Sire is installed in $install_dir/sire.app"
echo "## You can run a Sire python script called script.py by typing"
echo "## $install_dir/sire.app/bin/python script.py"
echo "##"
echo "## All Sire binaries are available in "
echo "## $install_dir/sire.app/bin"
echo "##"
echo "## You can copy or move this directory if you want, and Sire will still work"
echo "## Everything is contained in this directory, so you can delete Sire by"
echo "## deleting this directory"
echo "## e.g. rm -rf $install_dir/sire.app"
echo "##"
echo "## If you have never used Sire before, please take a look at the "
echo "## Sire tutorial at http://siremol.org/tutorial"
echo "##"
echo "## Please also take a look at the examples at http://siremol.org/examples"
echo "## for examples that you can download and use to learn Sire."
echo "##"
echo "## If you want to test this installation of Sire, please run"
echo "## $install_dir/sire.app/bin/sire_test"
echo "##"
echo "## If you have any problems or questions please get in touch"
echo "##"
echo "##############################################################################"
echo " "
