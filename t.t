   IFRAME: [1]https://www.googletagmanager.com/ns.html?id=GTM-WB38P87

   RDocumentation
     * [2]R Enterprise Training
     * [3]R package
     * [4]Leaderboard
     * [5]Sign in

   ____________________

sprintf

   From [6]base v3.4.3 by [7]R-core R-core@R-project.org
   0th

   Percentile

Use C-style String Formatting Commands

   A wrapper for the C function sprintf, that returns a character vector
   containing a formatted combination of text and variable values.

   Keywords
          [8]print, [9]character

Usage

sprintf(fmt, …)
gettextf(fmt, …, domain = NULL)

Arguments

   fmt
          a character vector of format strings, each of up to 8192 bytes.

   …
          values to be passed into fmt. Only logical, integer, real and
          character vectors are supported, but some coercion will be done:
          see the ‘Details’ section. Up to 100.

   domain
          see [10]gettext.

Details

   sprintf is a wrapper for the system sprintf C-library function.
   Attempts are made to check that the mode of the values passed match the
   format supplied, and R's special values (NA, Inf, -Inf and NaN) are
   handled correctly.

   gettextf is a convenience function which provides C-style string
   formatting with possible translation of the format string.

   The arguments (including fmt) are recycled if possible a whole number
   of times to the length of the longest, and then the formatting is done
   in parallel. Zero-length arguments are allowed and will give a
   zero-length result. All arguments are evaluated even if unused, and
   hence some types (e.g., "symbol" or "language", see [11]typeof) are not
   allowed.

   The following is abstracted from Kernighan and Ritchie (see
   References): however the actual implementation will follow the C99
   standard and fine details (especially the behaviour under user error)
   may depend on the platform.

   The string fmt contains normal characters, which are passed through to
   the output string, and also conversion specifications which operate on
   the arguments provided through …. The allowed conversion specifications
   start with a % and end with one of the letters in the set
   aAdifeEgGosxX%. These letters denote the following types:

   d, i, o, x, X
          Integer value, o being octal, x and X being hexadecimal (using
          the same case for a-f as the code). Numeric variables with
          exactly integer values will be coerced to integer. Formats d and
          i can also be used for logical variables, which will be
          converted to 0, 1 or NA.

   f
          Double precision value, in “fixed point” decimal notation of the
          form "[-]mmm.ddd". The number of decimal places ("d") is
          specified by the precision: the default is 6; a precision of 0
          suppresses the decimal point. Non-finite values are converted to
          NA, NaN or (perhaps a sign followed by) Inf.

   e, E
          Double precision value, in “exponential” decimal notation of the
          form [-]m.ddde[+-]xx or [-]m.dddE[+-]xx.

   g, G
          Double precision value, in %e or %E format if the exponent is
          less than -4 or greater than or equal to the precision, and %f
          format otherwise. (The precision (default 6) specifies the
          number of significant digits here, whereas in %f, %e, it is the
          number of digits after the decimal point.)

   a, A
          Double precision value, in binary notation of the form
          [-]0xh.hhhp[+-]d. This is a binary fraction expressed in hex
          multiplied by a (decimal) power of 2. The number of hex digits
          after the decimal point is specified by the precision: the
          default is enough digits to represent exactly the internal
          binary representation. Non-finite values are converted to NA,
          NaN or (perhaps a sign followed by) Inf. Format %a uses
          lower-case for x, p and the hex values: format %A uses
          upper-case.

          This should be supported on all platforms as it is a feature of
          C99. The format is not uniquely defined: although it would be
          possible to make the leading h always zero or one, this is not
          always done. Most systems will suppress trailing zeros, but a
          few do not. On a well-written platform, for normal numbers there
          will be a leading one before the decimal point plus (by default)
          13 hexadecimal digits, hence 53 bits. The treatment of
          denormalized (aka ‘subnormal’) numbers is very
          platform-dependent.

   s
          Character string. Character NAs are converted to "NA".

   %
          Literal % (none of the extra formatting characters given below
          are permitted in this case).

   Conversion by [12]as.character is used for non-character arguments with
   s and by [13]as.double for non-double arguments with f, e, E, g, G. NB:
   the length is determined before conversion, so do not rely on the
   internal coercion if this would change the length. The coercion is done
   only once, so if length(fmt) > 1 then all elements must expect the same
   types of arguments.

   In addition, between the initial % and the terminating conversion
   character there may be, in any order:

   m.n
          Two numbers separated by a period, denoting the field width (m)
          and the precision (n).

   -
          Left adjustment of converted argument in its field.

   +
          Always print number with sign: by default only negative numbers
          are printed with a sign.

   a space
          Prefix a space if the first character is not a sign.

   0
          For numbers, pad to the field width with leading zeros. For
          characters, this zero-pads on some platforms and is ignored on
          others.

   #
          specifies “alternate output” for numbers, its action depending
          on the type: For x or X, 0x or 0X will be prefixed to a non-zero
          result. For e, e, f, g and G, the output will always have a
          decimal point; for g and G, trailing zeros will not be removed.

   Further, immediately after % may come 1$ to 99$ to refer to a numbered
   argument: this allows arguments to be referenced out of order and is
   mainly intended for translators of error messages. If this is done it
   is best if all formats are numbered: if not the unnumbered ones process
   the arguments in order. See the examples. This notation allows
   arguments to be used more than once, in which case they must be used as
   the same type (integer, double or character).

   A field width or precision (but not both) may be indicated by an
   asterisk *: in this case an argument specifies the desired number. A
   negative field width is taken as a '-' flag followed by a positive
   field width. A negative precision is treated as if the precision were
   omitted. The argument should be integer, but a double argument will be
   coerced to integer.

   There is a limit of 8192 bytes on elements of fmt, and on strings
   included from a single %letter conversion specification.

   Field widths and precisions of %s conversions are interpreted as bytes,
   not characters, as described in the C standard.

   The C doubles used for R numerical vectors have signed zeros, which
   sprintf may output as -0, -0.000 ….

Value

   A character vector of length that of the longest input. If any element
   of fmt or any character argument is declared as UTF-8, the element of
   the result will be in UTF-8 and have the encoding declared as UTF-8.
   Otherwise it will be in the current locale's encoding.

Warning

   The format string is passed down the OS's sprintf function, and
   incorrect formats can cause the latter to crash the R process . R does
   perform sanity checks on the format, but not all possible user errors
   on all platforms have been tested, and some might be terminal.

   The behaviour on inputs not documented here is ‘undefined’, which means
   it is allowed to differ by platform.

References

   Kernighan, B. W. and Ritchie, D. M. (1988) The C Programming Language.
   Second edition, Prentice Hall. Describes the format options in table
   B-1 in the Appendix.

   The C Standards, especially ISO/IEC 9899:1999 for ‘C99’. Links can be
   found at [14]https://developer.r-project.org/Portability.html.

   man sprintf on a Unix-alike system.

See Also

   [15]formatC for a way of formatting vectors of numbers in a similar
   fashion.

   [16]paste for another way of creating a vector combining text and
   values.

   [17]gettext for the mechanisms for the automated translation of text.

Aliases

     * sprintf
     * gettextf

Examples

   library(base) # NOT RUN { ## be careful with the format: most things in
   R are floats ## only integer-valued reals get coerced to integer.
   sprintf("%s is %f feet tall\n", "Sven", 7.1) # OK try(sprintf("%s is %i
   feet tall\n", "Sven", 7.1)) # not OK sprintf("%s is %i feet tall\n",
   "Sven", 7 ) # OK ## use a literal % : sprintf("%.0f%% said yes (out of
   a sample of size %.0f)", 66.666, 3) ## various formats of pi :
   sprintf("%f", pi) sprintf("%.3f", pi) sprintf("%1.0f", pi)
   sprintf("%5.1f", pi) sprintf("%05.1f", pi) sprintf("%+f", pi)
   sprintf("% f", pi) sprintf("%-10f", pi) # left justified sprintf("%e",
   pi) sprintf("%E", pi) sprintf("%g", pi) sprintf("%g", 1e6 * pi) # ->
   exponential sprintf("%.9g", 1e6 * pi) # -> "fixed" sprintf("%G", 1e-6 *
   pi) ## no truncation: sprintf("%1.f", 101) ## re-use one argument three
   times, show difference between %x and %X xx <- sprintf("%1$d %1$x
   %1$X", 0:15) xx <- matrix(xx, dimnames = list(rep("", 16), "%d%x%X"))
   noquote(format(xx, justify = "right")) ## More sophisticated:
   sprintf("min 10-char string '%10s'", c("a", "ABC", "and an even longer
   one")) # } # NOT RUN { ## Platform-dependent bad example from qdapTools
   1.0.0: ## may pad with spaces or zeroes. sprintf("%09s", month.name) #
   } # NOT RUN { n <- 1:18 sprintf(paste0("e with %2d digits = %.", n,
   "g"), n, exp(1)) ## Using arguments out of order sprintf("second
   %2$1.0f, first %1$5.2f, third %3$1.0f", pi, 2, 3) ## Using asterisk for
   width or precision sprintf("precision %.*f, width '%*.3f'", 3, pi, 8,
   pi) ## Asterisk and argument re-use, 'e' example reiterated: sprintf("e
   with %1$2d digits = %2$.*1$g", n, exp(1)) ## re-cycle arguments
   sprintf("%s %d", "test", 1:3) ## binary output showing
   rounding/representation errors x <- seq(0, 1.0, 0.1); y <-
   c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1) cbind(x, sprintf("%a", x),
   sprintf("%a", y)) # }
   Documentation reproduced from package base, version 3.4.3, License:
   Part of R 3.4.3

Community examples

   Looks like there are no examples yet.

Post a new example:


   ____________________________________________________________
   ____________________________________________________________
   ____________________________________________________________
   ____________________________________________________________

   [18](BUTTON) Submit your example

   [19]API documentation
   R package
   Rdocumentation.org

   Created by [20]DataCamp.com

   Learn R at work [21]Try it free

References

   Visible links
   1. https://www.googletagmanager.com/ns.html?id=GTM-WB38P87
   2. https://www.datacamp.com/groups/business
   3. https://github.com/datacamp/Rdocumentation
   4. https://www.rdocumentation.org/trends
   5. https://www.rdocumentation.org/login?rdr=/packages/base/versions/3.4.3/topics/sprintf
   6. https://www.rdocumentation.org/packages/base/versions/3.4.3
   7. https://www.rdocumentation.org/collaborators/name/R-core R-core@R-project.org
   8. https://www.rdocumentation.org/search/keywords/print
   9. https://www.rdocumentation.org/search/keywords/character
  10. https://www.rdocumentation.org/link/gettext?package=base&version=3.4.3
  11. https://www.rdocumentation.org/link/typeof?package=base&version=3.4.3
  12. https://www.rdocumentation.org/link/as.character?package=base&version=3.4.3
  13. https://www.rdocumentation.org/link/as.double?package=base&version=3.4.3
  14. https://developer.r-project.org/Portability.html
  15. https://www.rdocumentation.org/link/formatC?package=base&version=3.4.3
  16. https://www.rdocumentation.org/link/paste?package=base&version=3.4.3
  17. https://www.rdocumentation.org/link/gettext?package=base&version=3.4.3
  18. https://www.rdocumentation.org/modalLogin
  19. https://www.rdocumentation.org/docs
  20. https://www.datacamp.com/
  21. https://www.datacamp.com/groups/business

   Hidden links:
  23. https://www.rdocumentation.org/
  24. https://github.com/datacamp/rdocumentation
  25. https://github.com/datacamp/rdocumentation-app
