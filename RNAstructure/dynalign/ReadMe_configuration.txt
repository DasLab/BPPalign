The format of variables in the configuration file is
very straightforward; each one goes on a separate
line following the format

[variable] = [value]

Leading and trailing whitespace on variables and
values is ignored.  Variables and values cannot
contain internal whitespace.

Syntax errors will produce a warning to standard
error, but are otherwise ignored.

If a required variable is not defined, Dynalign will
exit.  Not all variables are required; for example,
if `savefile' is not defined, Dynalign will not
write to a savefile.

If a required variable is defined twice, the
designed behavior is to accept the last definition.

Empty lines and lines beginning with # are treated
as comments and ignored.

A sample configuration file can be found in the
dynalign source directory as 'sample.conf'.  Note
that since some variables are not required, they are
not defined here.
