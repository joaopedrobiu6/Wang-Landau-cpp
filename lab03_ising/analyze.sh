#!/bin/bash

if [ $# -lt 2 ]; then
    echo "Plot the quantity VAR against the temperature T of input files."
    echo ""
    echo "usage: $0 VAR FILE1 FILE2 ..."
    echo ""
    echo "    VAR is the column number of the quantity to plot"
    echo "    FILE1 ... are filenames to plot information from"
    echo ""
    echo "Columns: T, m, abs(m), chi, e, c, binder, L, ..."
    echo "Hint: Edit the script to also plot the exact (2D) solution."
    exit 2;
fi

if [[ ! $1 -gt 0 ]] || [[ ! $1 -le 9 ]]; then
    echo "Error: Argument 1 must be a column number to plot against the temperature."
    echo ""
    echo "Columns: T, m, abs(m), chi, e, c, binder, L, ..."
    exit 1;
fi

ARGS=$@

T=1
VAR=$1
FILES="${ARGS[@]:1}"
TMPFN=".tmpout"

i=1;
for file in $FILES; do
    tmp=$(printf "%s%d" $TMPFN $i);
    tmpfiles+=("$tmp");
    i=$(echo $i + 1 | bc);

    # copy T and variable columns to a new file
    awk "{print \$${T}, \$${VAR}}" $file > $tmp
done;

# uncomment below to also plot a variable from the exact solution
# columns: T, e, c, m (ie. $4 below refers to the magnetization)
#tmpexact=.tmpexact
#awk '{print $1, $4}' exact_solution.data > $tmpexact

# open the files with xmgrace and clean up
xmgrace -free ${tmpfiles[@]} $tmpexact
rm ${tmpfiles[@]} $tmpexact
