
all: BamToSnpAD.pdf

BamToSnpAD.pdf: BamToSnpAD.1
	groff -k -T ps -man BamToSnpAD.1 > BamToSnpAD.ps
	ps2pdf BamToSnpAD.ps BamToSnpAD.pdf
