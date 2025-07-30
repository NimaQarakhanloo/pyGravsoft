/****************************************************************************

    =head1 NAME

    gggg  -  Gravsoft Grid to GMT Grid


    =head1 SYNOPSIS
    
    gggg  <infile>  <tmpfile>  <outfile>

    
    =head1 DESCRIPTION

    gggg organizes the conversion of a Gravsoft grid file to an GMT grid file.
    The conversion is not carried out directly by gggg. Rather it assembles
    the proper arguments and operands for the Generic Mapping Tools (GMT)
    XYZ2GRD program, builds a command line, and calls the system to do the
    dirty work.

    To this end, we need a file holding a slightly changed form of the Gravsoft
    grid. This form is written to the file <tmpfile>


    =head1 METHOD

    -
    
    =head1 RETURN 
    
    0 on success, non-0 on failure


    =head1 BUGS

    gggg does not check for the presence of XYZ2GRD on the path.


    =head1 HISTORY

    Thomas Knudsen, thk@kms.dk,  2001-01-23:  Original version
    Thomas Knudsen, thk@kms.dk,  2001-02-16:  fclose and remove tempfile
    Thomas Knudsen, thk@kms.dk,  2001-02-20:  Renamed to GGGG, POD added

    =cut
       
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    char dx[128], dy[128], w[128], e[128], s[128], n[128], buf[8192];
    FILE *finp, *ftmp;

    if (argc!=4)
        return fprintf(stderr, "\n\ngggg - convert Gravsoft Grid to GMT Grid format.\n\n"
                               "syntax: gggg <infile> <tempfile> <outfile> - bye!\n"
                               "<infile>   is a Gravsoft grid file\n"
                               "<tempfile> is a temporary file\n"
                               "<outfile>  is the GMT format grid\n\n"
                               "Thomas Knudsen, thk@kms.dk, 2001-02-20\n\n"                   );

    finp = fopen(argv[1], "rt");
    if (0==finp)
        return fprintf(stderr, "gggg: cannot open input file `%s' - bye!\n", argv[1]);
    ftmp = fopen(argv[2], "rt");
    if (0!=ftmp)
        return fprintf(stderr, "gggg: temporary file `%s' exists - bye!\n", argv[2]);
    ftmp = fopen(argv[2], "wt");
    if (0==ftmp)
        return fprintf(stderr, "gggg: cannot open temporary file `%s' - bye!\n", argv[2]);

    /* read gravsoft header */
    fscanf(finp, "%s", s);
    fscanf(finp, "%s", n);
    fscanf(finp, "%s", w);
    fscanf(finp, "%s", e);
    fscanf(finp, "%s", dy);
    fscanf(finp, "%s", dx);

    /* loop over all grid nodes - dump them, newline separated, to ftmp */
    for (fscanf(finp, "%s", buf); !feof(finp);    fscanf(finp, "%s", buf))
        fprintf(ftmp, "%s\n", buf);
    fclose(ftmp);

    /* build gmt command for converting */
    sprintf(buf, "xyz2grd %s -G%s -I%s/%s -R%s/%s/%s/%s -Z", argv[2], argv[3], dx,dy, w,e,s,n);
    system(buf);

    /* clean up */
    remove(argv[2]);

    return 0;
}
