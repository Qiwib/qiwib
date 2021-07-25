# The QiwiB program package

QiwiB is a program written in GNU Octave (an open source Matlab clone) to solve the many-particle Schroedinger equation in one dimension. At the moment, it simulates the full quantum many-body physics of ultracold bosons. For this, it uses the MCTDHB approach given in the following reference:
O. E. Alon, A. I. Streltsov, and L. S. Cederbaum, Phys. Rev. A 77, 033613 (2008).

QiwiB was originally written by Thomas Ernst under the mentorship of Joachim Brand. It was further developed by members of the CTCP group at Massey University Auckland ([link](http://ctcp.massey.ac.nz/~brand)).

QiwiB is written in Octave and C++. Included in the QiwiB package is a comprehensive HTML manual.

## How to use

If this your first time using QiwiB, please execute setup.sh in the qiwib/bin/ direcory first.
This and more information can be found in the html documentation in the qiwib/doc/ directory
(open index.html with any browser).

Explanation of the directories:

	bin/			:	qiwib and other scripts for analysis and plotting (can be run from command line)
	doc/			:	documentation (installation instructions, writing input file ...),
					open index.html with your favourite internet browser for a quick introduction to QiwiB
	lib/			:	all octave files
	sample_input		:	example input files for qiwib

To show the version of your running installation, just start qiwib without any options.

For questions do not hesitate to contact the author.

Because QiwiB is still in its early stages, we cannot guarantee it to be completely bug free. If you discover a bug please let us know by filing a bug report in the "Issues" tab.

## License

This program is provided under the MIT X License. Further details can be found in the LICENSE file
that should have been distributed with QiwiB or in the header of every single source file.

## News

    03.10.2011: Another bug fix release
    06.09.2011: Minor bug fixes
    31.08.2011: Speed improvements and bug fix for the calculation of natural orbitals
    30.08.2011: Important bug fix for QiwiB
    25.08.2011: Inital release of QiwiB
