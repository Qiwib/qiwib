#! /bin/bash
# Setup the qiwib and compile_qiwib files to include the right link to the octave binary

set -u
default_octave=`which octave`
echo
echo "QiwiB Setup"
echo
echo "-----------------------------------------------------------------------------"
echo
echo "Specify the location of the octave binary including the absolute path (i.e. /usr/bin/octave)"
echo
echo "The terminal command 'which octave' gives the following default location for your computer:"
echo $default_octave
echo
echo "Please enter now the location or leave blank to use the default one followed by [ENTER]:"

read octave_readin
if [ -n "$octave_readin" ]
then
  octave_bin=$octave_readin
else
  octave_bin=$default_octave
fi

octave_line="#! $octave_bin -Hq"

echo



echo "-----------------------------------------------------------------------------"
setup_path=$(cd "$(dirname "$0")"; pwd)
cd $setup_path

sed -e "1d" qiwib > temp.scr
rm -f qiwib
echo "$octave_line" >> qiwib
cat temp.scr >> qiwib

sed -e "1d" compile_qiwib > temp.scr
rm -f compile_qiwib
echo "$octave_line" >> compile_qiwib
cat temp.scr >> compile_qiwib

sed -e "1d" plot_density > temp.scr
rm -f plot_density
echo "$octave_line" >> plot_density
cat temp.scr >> plot_density

sed -e "1d" plot_g1 > temp.scr
rm -f plot_g1
echo "$octave_line" >> plot_g1
cat temp.scr >> plot_g1

sed -e "1d" plot_g2 > temp.scr
rm -f plot_g2
echo "$octave_line" >> plot_g2
cat temp.scr >> plot_g2

sed -e "1d" plot_fock > temp.scr
rm -f plot_fock
echo "$octave_line" >> plot_fock
cat temp.scr >> plot_fock

sed -e "1d" plot_NO > temp.scr
rm -f plot_NO
echo "$octave_line" >> plot_NO
cat temp.scr >> plot_NO

sed -e "1d" plot_fock_NO > temp.scr
rm -f plot_fock_NO
echo "$octave_line" >> plot_fock_NO
cat temp.scr >> plot_fock_NO

sed -e "1d" plot_phi > temp.scr
rm -f plot_phi
echo "$octave_line" >> plot_phi
cat temp.scr >> plot_phi

sed -e "1d" plot_phi_pop > temp.scr
rm -f plot_phi_pop
echo "$octave_line" >> plot_phi_pop
cat temp.scr >> plot_phi_pop

sed -e "1d" create_hilbert > temp.scr
rm -f create_hilbert
echo "$octave_line" >> create_hilbert
cat temp.scr >> create_hilbert
rm -f temp.scr

chmod u+x qiwib
chmod u+x compile_qiwib
chmod u+x plot_density
chmod u+x plot_g1
chmod u+x plot_g2
chmod u+x plot_fock
chmod u+x plot_fock_NO
chmod u+x plot_NO
chmod u+x plot_phi
chmod u+x plot_phi_pop
chmod u+x create_hilbert

./compile_qiwib

echo
echo "-----------------------------------------------------------------------------"
echo
echo "Please check the compilation output above for errors!"
echo
echo "-----------------------------------------------------------------------------"
echo
echo "Optional: Append the following two lines to ~/.bashrc"
echo "to set PATH to QiwiB scripts (qiwib scripts can be executed"
echo "from any directory without specifying the full path):"
echo
echo PATH=$setup_path:"\$PATH"
echo "export PATH"
echo
echo "Shall this setup script do that for you? (y/n)"
read bash_text
if [ "$bash_text" == 'y' ]; then
	echo >> ~/.bashrc
	echo >> ~/.bashrc
	echo "# QiwiB" >> ~/.bashrc
	echo "#" >> ~/.bashrc
	echo "export PATH=$setup_path:\$PATH" >> ~/.bashrc
	source ~/.bashrc
fi
echo
echo "After appending restart your bash session or execute 'source ~/.bashrc'"
echo "from the bash command line."
echo
echo "Setup finished!"
echo

